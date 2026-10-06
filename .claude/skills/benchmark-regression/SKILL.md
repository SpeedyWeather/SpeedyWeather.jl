---
name: benchmark-regression
description: Check whether the current branch or PR of SpeedyWeather.jl has a performance regression. Runs the quick debug mode of the benchmark suite (PrimitiveWet resolution sweep, 8 layers, truncation ≤ 128) on the branch's HEAD and as references on origin/main, the latest release and the latest benchmarked revision; if the branch is significantly slower, bisects to the culprit commit and writes a markdown summary with a table, the culprit and a plain-language hypothesis, ready to paste as a GitHub comment. Use when asked to check a branch or PR for performance regressions, compare benchmark performance across versions, or find the commit that made the model slower. Optional argument - an architecture (cpu, gpu, amdgpu; default cpu).
---

# Benchmark regression check for a branch

The measuring, comparing and bisecting is done by
`SpeedyWeather/benchmark/regression/regression.jl` — the same script the scheduled CI job
(`.github/workflows/benchmark_regression.yml`) runs on main. This skill runs it on the current
branch, sanity-checks the result, adds a hypothesis for a culprit and writes a markdown file the
user can paste into a GitHub comment. It never posts anything to GitHub itself.

## Ground rules

- `regression.jl` benchmarks every revision in a detached worktree under the work directory,
  with the benchmark harness of the current checkout. **Never check out other revisions in the
  user's working tree** yourself.
- It benchmarks the **committed** `HEAD`. If `git status --porcelain` shows uncommitted changes
  to source files, tell the user they are not included and ask whether to continue or commit first.
- **One benchmark at a time**, and no other heavy work (tests, builds, other agents) in parallel —
  contention makes the timings meaningless.
- On Apple silicon (M3) each revision takes ~10 min (~7.5 min benchmark, 2–4 min precompile):
  ~40 min for the four revisions without a regression; a confirmed regression adds ~30 min and
  a bisection ~10 min per step (≈ log₂ of the commits in range). Run it with `run_in_background`
  and wait for the completion notification; don't poll.
- Every Bash call starts a fresh shell: repeat the variable block at the top of each call.

```bash
REPO=$(git rev-parse --show-toplevel)
W=<scratchpad dir if the session has one, else ${TMPDIR:-/tmp}>/benchmark-regression
OUT=$REPO/SpeedyWeather/benchmark/regression/results     # gitignored, survives the session
ARCH=cpu        # or gpu / amdgpu from the skill arguments
R="julia --startup-file=no $REPO/SpeedyWeather/benchmark/regression/regression.jl"
```

## 1. Check

```bash
rm -rf $OUT && $R check --candidate HEAD --arch $ARCH --workdir $W --output-dir $OUT --bisect
```

This benchmarks the branch (`HEAD`) and as references `origin/main`, the latest release tag and
the latest benchmarked revision of this architecture (the first-parent commit that introduced its
stored results in `SpeedyWeather/benchmark/assets/benchmark_results.json`). A **regression** is a
geometric-mean SYPD ratio branch / reference below 0.85 for all configurations or for one of the
two transforms (LT+FFT, MT). It is confirmed by benchmarking both a second time (best of both runs
counts) and then bisected (first parent) from the most recent regressed reference, or from where
the branch forked off it. On `main` itself the candidate is main and main is no reference.

Results: `$OUT/report.md` (table, verdicts, bisection) and `$OUT/summary.json`
(`regression`, `comparisons`, `bisect.culprit`, `notes`), plus one JSON + `.log` per run.
Revisions that cannot run with the current harness are listed under notes.

## 2. Sanity-check a culprit

Only if there is one. In the bisection table of `report.md`, the culprit's ratio must be clearly
below the cutoff and its parent's clearly above. If either is within ~5% of the cutoff, benchmark
both again and compare:

```bash
$R benchmark <sha> $OUT/recheck-<sha7>.json --arch $ARCH --workdir $W
$R table --candidate culprit --result parent=<files,...> --result culprit=<files,...>
```

## 3. Hypothesis

Only if there is a culprit. Read what it changed — `git -C "$REPO" show --stat <sha>`, its PR
(`gh pr view <N>` with the `#N` from the commit title, if any) and the diff of code that runs every
time step: `dynamics/`, `parameterizations/`, time stepping, `SpeedyTransforms`, `RingGrids`,
`LowerTriangularArrays`. Typical causes: new work or allocations inside the time loop, type
instabilities, changed default parameters or components (e.g. a new parameterization enabled by
default, shorter time step), changed loop order or memory layout. Whether only MT or only LT+FFT
regressed points at the transforms; both regressed points at dynamics or physics. Write 2–3
short sentences in plain, simple language and phrase it as a hypothesis, not a finding.

## 4. Write the comment

Write `$OUT/comment.md`: the content of `report.md`, and if there is a culprit, append

```markdown
**Hypothesis:** <2–3 plain sentences>
```

after the bisection section. Do not add anything else — it is meant to be pasted as is into a
GitHub comment. Then give the user the path to `comment.md`, the one-line verdict (the headline of
the report) and the hypothesis if any. Keep `$OUT` for later inspection; `regression.jl` removes its
worktrees itself.
