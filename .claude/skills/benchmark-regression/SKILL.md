---
name: benchmark-regression
description: Check whether SpeedyWeather.jl's main branch has a performance regression. Runs the quick debug mode of the benchmark suite (PrimitiveWet resolution sweep, truncation ≤ 128) on origin/main, the latest release and the latest benchmarked revision; if main is significantly slower, bisects main to the culprit commit and files a short GitHub issue with a table, the culprit and a plain-language hypothesis. Use when asked to check for performance regressions, compare benchmark performance across versions, or find the commit that made the model slower. Optional arguments - an architecture (cpu, gpu, amdgpu; default cpu) and "dry-run" to report without filing an issue.
---

# Benchmark regression check

Compares the speed of `origin/main` against the latest release and the latest benchmarked
revision with `manual_benchmarking.jl --debug`, bisects a regression to its culprit commit and
files a GitHub issue. Expect roughly 15–30 min without and 1–2 h with a bisection; most of it is
precompiling every revision. Run long commands with `run_in_background` and wait for completion.

## Ground rules

- **Never check out another revision in the user's working tree.** Every revision is benchmarked
  in a detached worktree under `$W/trees` (bisection in `$W/bisect`); `run_debug_benchmark.sh`
  handles this and always uses the benchmark harness of the current checkout, so all revisions are
  measured identically.
- **One benchmark at a time**, and no other heavy work (tests, builds, other agents) in parallel —
  contention makes the timings meaningless. Same machine, arch and Julia for every revision.
- Every Bash call starts a fresh shell: repeat the variable block below at the top of each call.
- Benchmarks are noisy. Never call a regression from a single run; confirm it (step 4).

```bash
REPO=$(git rev-parse --show-toplevel)
S=$REPO/.claude/skills/benchmark-regression/scripts
export SPEEDY_BENCH_WORKDIR=<scratchpad dir if the session has one, else ${TMPDIR:-/tmp}>/benchmark-regression
W=$SPEEDY_BENCH_WORKDIR
ARCH=cpu        # or gpu / amdgpu from the skill arguments
RJ="julia --startup-file=no $S/regression.jl"   # own env in scripts/, instantiated on first use
```

## 1. Revisions

```bash
git -C "$REPO" fetch origin --tags --quiet
MAIN=$(git -C "$REPO" rev-parse origin/main)
RELEASE=$(git -C "$REPO" tag --list 'v[0-9]*' --sort=-v:refname | head -1)
BENCHMARKED=$($RJ benchmarked-commit --arch $ARCH --repo "$REPO")
```

`benchmarked-commit` finds the first-parent commit on main that introduced the currently stored
results of this architecture in `SpeedyWeather/benchmark/assets/benchmark_results.json`.
Drop a revision that resolves to the same commit as another one. If a revision fails to run with
the current harness (e.g. it predates an API the harness uses), say so in the report and continue
with the remaining ones.

## 2. Benchmark

Sequentially, in one background command (each run prints its table, logs go to `<json>.log`):

```bash
$S/run_debug_benchmark.sh "$BENCHMARKED" $W/results/benchmarked-1.json $ARCH
$S/run_debug_benchmark.sh "$RELEASE"     $W/results/release-1.json     $ARCH
$S/run_debug_benchmark.sh "$MAIN"        $W/results/main-1.json        $ARCH
```

## 3. Compare

```bash
$RJ table --candidate main \
    --result benchmarked=$W/results/benchmarked-1.json \
    --result release=$W/results/release-1.json \
    --result main=$W/results/main-1.json
```

This prints a markdown table (SYPD per configuration, change of main vs each reference, geometric
means for all configurations, LT+FFT and MT separately) and a verdict per reference. A
**regression** is a geometric-mean SYPD ratio below 0.85 (main ≥ 15% slower) in any group.

No regression → skip to step 7.

## 4. Confirm

Re-run main and every reference it regressed against once more (`main-2.json`, `release-2.json`, …),
then re-run `table` with both files per revision, e.g. `--result main=$W/results/main-1.json,$W/results/main-2.json`
(the best SYPD per configuration is used). Only a regression that persists counts. If it
disappears, report it as noise and stop.

## 5. Bisect

Good revision = the most recent reference that main regressed against (`git merge-base
--is-ancestor` tells which is newer; if it is not an ancestor of main use `git merge-base` with
main). Take `GROUP` and `CUTOFF` from the verdict line of that reference. The cutoff is the
geometric midpoint between the reference and main, so single-run noise rarely flips a verdict.

```bash
git -C "$REPO" worktree add --detach $W/bisect $MAIN
git -C $W/bisect bisect start --first-parent $MAIN $GOOD
git -C $W/bisect bisect run $S/bisect_step.sh $W/results/<ref>-1.json,$W/results/<ref>-2.json $CUTOFF $GROUP $ARCH
git -C $W/bisect bisect log > $W/results/bisect.log
git -C $W/bisect bisect reset
```

`bisect_step.sh` benchmarks each step into `$W/results/bisect-<sha>.json` and skips (exit 125)
revisions that fail to run. Sanity-check the culprit: its ratio (from the bisect output) must be
clearly below the cutoff and its parent's clearly above. If either is within ~5% of the cutoff,
re-run both once more with `run_debug_benchmark.sh` and judge on the best of the runs.
If bisect ends with only skipped commits, report the remaining range instead of a single culprit.

## 6. Hypothesis

Read what the culprit changed — `git -C "$REPO" show --stat <sha>`, the PR (`gh pr view <N>` with
the `#N` from the commit title) and the diff of code that runs every time step: `dynamics/`,
`parameterizations/`, time stepping, `SpeedyTransforms`, `RingGrids`, `LowerTriangularArrays`.
Typical causes: new work or allocations inside the time loop, type instabilities, changed default
parameters or components (e.g. a new parameterization enabled by default, shorter time step),
changed loop order or memory layout. Whether only MT or only LT+FFT regressed points at the
transforms; both regressed points at dynamics or physics. Write 2–3 short sentences in plain,
simple language and phrase it as a hypothesis, not a finding.

## 7. Report

**No regression:** tell the user in a few lines with the table. No issue.

**Regression** (unless the user asked for a dry run): first check for an existing open issue,
`gh issue list --repo SpeedyWeather/SpeedyWeather.jl --state open --search "Performance regression in:title"`.
If one already names the same culprit, add a comment with the new numbers instead of a new issue.
Otherwise file one with label `performance :rocket:`:

```bash
gh issue create --repo SpeedyWeather/SpeedyWeather.jl --label "performance :rocket:" \
    --title "Performance regression: main ~X% slower than <reference> since #N" --body-file $W/issue.md
```

Keep the body short — this template, nothing more:

```markdown
`main` (<sha7>) is **~X% slower** than <reference> (<sha7>) in the quick benchmark
(`manual_benchmarking.jl --debug`: PrimitiveWet, T ≤ 128, <arch label>). Higher SYPD is better.

<table from regression.jl, the main and reference columns are enough if it gets wide>

**Culprit (git bisect):** <sha7> <commit title> (#N)

**Hypothesis:** <2–3 plain sentences>

<details><summary>How this was measured</summary>

<machine, Julia version, threads>; each revision benchmarked with the same harness in its own
worktree, best of N runs per configuration; bisect cutoff <CUTOFF> on the <GROUP> geometric mean.
</details>
```

Give the user the issue link and the one-line verdict.

## 8. Clean up

```bash
git -C $W/bisect bisect reset 2>/dev/null
for tree in $W/trees/* $W/bisect; do [ -d "$tree" ] && git -C "$REPO" worktree remove --force "$tree"; done
git -C "$REPO" worktree prune
```

Keep `$W/results` (JSON + logs) so the numbers can be inspected later.
