---
name: benchmark-regression
description: Check whether SpeedyWeather.jl's main branch has a performance regression. Runs the quick debug mode of the benchmark suite (PrimitiveWet resolution sweep, truncation ≤ 128) on origin/main, the latest release and the latest benchmarked revision; if main is significantly slower, bisects main to the culprit commit and files a short GitHub issue with a table, the culprit and a plain-language hypothesis. Use when asked to check for performance regressions, compare benchmark performance across versions, or find the commit that made the model slower. Optional arguments - an architecture (cpu, gpu, amdgpu; default cpu) and "dry-run" to report without filing an issue.
---

# Benchmark regression check

The measuring, comparing and bisecting is done by
`SpeedyWeather/benchmark/regression/regression.jl` — the same script the scheduled CI job
(`.github/workflows/benchmark_regression.yml`) runs. This skill runs it, sanity-checks the result,
adds a hypothesis for the culprit and files the GitHub issue. Run `regression.jl` without arguments
for its full usage.

## Ground rules

- `regression.jl` benchmarks every revision in a detached worktree under the work directory, with
  the benchmark harness of the current checkout. **Never check out other revisions in the user's
  working tree** yourself.
- **One benchmark at a time**, and no other heavy work (tests, builds, other agents) in parallel —
  contention makes the timings meaningless.
- Expect 30–60 min without and 1–3 h with a bisection, mostly precompiling every revision. Run it
  with `run_in_background` and wait for the completion notification; don't poll.
- Every Bash call starts a fresh shell: repeat the variable block at the top of each call.

```bash
REPO=$(git rev-parse --show-toplevel)
W=<scratchpad dir if the session has one, else ${TMPDIR:-/tmp}>/benchmark-regression
ARCH=cpu        # or gpu / amdgpu from the skill arguments
R="julia --startup-file=no $REPO/SpeedyWeather/benchmark/regression/regression.jl"
```

## 1. Check

```bash
$R check --arch $ARCH --workdir $W --output-dir $W/results --bisect
```

This benchmarks `origin/main`, the latest release tag and the latest benchmarked revision of this
architecture (the first-parent commit that introduced its stored results in
`SpeedyWeather/benchmark/assets/benchmark_results.json`). A **regression** is a geometric-mean
SYPD ratio main / reference below 0.85 for all configurations or for one of the two transforms
(LT+FFT, MT). It is confirmed by benchmarking main and that reference a second time (best of both
runs counts) and then bisected (first parent) from the most recent regressed reference.

Results: `$W/results/report.md` (table, verdicts, bisection) and `summary.json`
(`regression`, `comparisons`, `bisect.culprit`, `notes`), plus one JSON + `.log` per run.
Revisions that cannot run with the current harness are listed under `notes`; mention them.

**No regression** → tell the user in a few lines with the table from `report.md`. Done.

## 2. Sanity-check the culprit

In the bisection table of `report.md`, the culprit's ratio must be clearly below the cutoff and
its parent's clearly above. If either is within ~5% of the cutoff, benchmark both again and compare:

```bash
$R benchmark <sha> $W/results/recheck-<sha7>.json --arch $ARCH --workdir $W
$R table --candidate culprit --result parent=<files,...> --result culprit=<files,...>
```

If the bisection ended with only skipped commits, report that list instead of a single culprit.

## 3. Hypothesis

Read what the culprit changed — `git -C "$REPO" show --stat <sha>`, the PR (`gh pr view <N>` with
the `#N` from the commit title) and the diff of code that runs every time step: `dynamics/`,
`parameterizations/`, time stepping, `SpeedyTransforms`, `RingGrids`, `LowerTriangularArrays`.
Typical causes: new work or allocations inside the time loop, type instabilities, changed default
parameters or components (e.g. a new parameterization enabled by default, shorter time step),
changed loop order or memory layout. Whether only MT or only LT+FFT regressed points at the
transforms; both regressed points at dynamics or physics. Write 2–3 short sentences in plain,
simple language and phrase it as a hypothesis, not a finding.

## 4. Issue

Skip this step if the user asked for a dry run. First look for an open issue:
`gh issue list --repo SpeedyWeather/SpeedyWeather.jl --state open --search "Performance regression in:title"`.
If one already names the same culprit, comment on it with the new numbers instead. Otherwise:

```bash
gh issue create --repo SpeedyWeather/SpeedyWeather.jl --label "performance :rocket:" \
    --title "Performance regression: main ~X% slower than <reference> since #N" --body-file $W/issue.md
```

Keep the body short — this template, nothing more:

```markdown
`main` (<sha7>) is **~X% slower** than <reference> (<sha7>) in the quick benchmark
(`manual_benchmarking.jl --debug`: PrimitiveWet, T ≤ 128, <arch label>). Higher SYPD is better.

<table from report.md; drop the columns of a reference that did not regress if it gets wide>

**Culprit (git bisect):** <sha7> <commit title> (#N)

**Hypothesis:** <2–3 plain sentences>

<details><summary>How this was measured</summary>

<the machine line of report.md>; `regression.jl check --bisect`, best of N runs per configuration,
bisect cutoff <cutoff> on the <group> geometric mean.
</details>
```

Give the user the issue link and the one-line verdict. Keep `$W/results` for later inspection;
`regression.jl` removes its worktrees itself.
