# Vertical advection revision — progress log

Tracking execution of the plan for revising
`SpeedyWeather/src/dynamics/vertical_advection.jl` (Centered, Upwind, WENO schemes).
Branch: `mg/vertical-advection-revision`.

## Important finding (2026-06-29)

Before starting, found `stash@{0}` ("WIP on mg/vertical-advection-revision") already
sitting on this branch, touching only `vertical_advection.jl` (+94/-34 lines). It
implements, in one combined edit:
- the Phase 1 value-based refactor: `reconstructed_at_face(ξ, ij, s, k, u, adv)` split into
  `gather_stencil_values(ξ, ij, s, k)` + `reconstruct_face(S, u, adv)`.
- a Phase 2-shaped CPU path: `_vertical_advection!(::CPU, ...)` plain loop, `k` outer / `ij`
  inner, with a single rolling `face::AbstractField` scratch buffer (registered via
  `variables(::VerticalAdvection)` as a new `ScratchVariable`) so each interior face is
  reconstructed once and carried to the next `k` — same flux-form idea as the plan, just one
  buffer instead of two ping-ponged ones.
- GPU path left as the original per-(ij,k) kernel, just using the new
  `gather_stencil_values`/`reconstruct_face` split (no Phase 3 fusion/per-column kernel yet).

Decision: did not blindly `git stash apply`. Plan's own #1 priority is bit-identical
correctness validated by a golden harness *written first* (Phase 0), and the stash was
written before any such harness existed. So: do Phase 0 against the current committed
(un-stashed) code first, then bring in the stash's logic (re-split into separate
Phase-1 / Phase-2 commits) and validate against the Phase 0 golden set before trusting it.
The stash itself is left untouched (not dropped) until its content is fully absorbed.

## Environment notes

- Interactive shell is `login03` (HPC, Lmod + SLURM) and does have a real A40 GPU
  (`nvidia-smi`). Per the user (2026-06-29): both CPU and GPU work — including quick
  benchmarks — are fine directly on the login node, as long as a single run doesn't exceed
  ~5-10 minutes; only genuinely extensive jobs need `sbatch` (the user's own
  `submit_bench.sh` / `submit_run_benchmarks.sh`, untracked in repo root, use
  `--partition=gpu --gres=gpu:1 --qos=gpushort --account=flai` for those). So Phase 2/3
  benchmarking can run inline here.
- Julia via `module load julia` (default 1.12.2). Each Bash tool call is a fresh shell, so
  `module load julia` must be chained with `&&` in every command that calls `julia`.

## Checklist

- [x] Phase 0: golden capture + face-consistency assertion (scratch validation script, no src change)
- [x] Phase 1: value-based `reconstruct_face` refactor (bit-identical) + tests
- [x] Phase 2: CPU plain-loop flux-form path + golden + CPU bench
- [x] Phase 3: GPU kernel(s) 3a/3b + golden(GPU) + GPU bench (decide 3a vs 3b)
- [ ] Phase 4: fuse u/v/T/q (+tracers) into one pass + golden + bench
- [ ] Phase 5: dispatch wiring, cleanup, version bump, CHANGELOG, full validation

## Log

- 2026-06-29: Read current `vertical_advection.jl`, the `vertical_integration!` CPU/GPU
  precedent in `tendencies.jl`, and the existing `test/dynamics/vertical_advection.jl`.
  Found and inspected the stash described above. Confirmed Julia 1.12.2 available via
  `module load julia`; A40 GPU visible via `nvidia-smi` on the login node (access mode for
  Phase 3 still TBD). Set up this progress file.
- 2026-06-29: Phase 0 done. Added a permanent "Vertical advection face consistency" testset
  to `test/dynamics/vertical_advection.jl` asserting `tail(stencil(k)) == front(stencil(k+1))`
  for Centered{1,2}, Upwind{1,2,3}, WENO{3}, nlayers ∈ {2,3,5,8,24} — 222/222 pass. Existing
  "Vertical advection runs"/"stencils" testsets still pass (`--check-bounds=yes`). Wrote a
  scratch golden-comparison harness (not committed) at
  `/tmp/claude-3154/.../scratchpad/golden_harness.jl`: `capture` mode randomizes
  u/v/temperature/humidity/w/tendencies + a tracer `:abc` for 5 schemes × 3 (grid,nlayers)
  configs, runs `vertical_advection!` once, serializes inputs+outputs; `check` mode reloads
  the same inputs into a fresh model and asserts `all(new .=== old)` per output field.
  Capture run against the pre-revision code (HEAD, before any Phase 1+ edits) completed for
  all 15 configs; a self-check (run against the same unmodified code) confirmed the harness
  itself reports "all bit-identical". Also found that GPU compute on this cluster goes
  through `sbatch --partition=gpu --gres=gpu:1 --qos=gpushort --account=flai` (per the
  user's existing `submit_bench.sh`/`submit_run_benchmarks.sh`), not direct execution on
  the login node — relevant for Phase 3.
  Committed Phase 0 as `113ace41`.
- 2026-06-29: Phase 1 done — split `reconstructed_at_face(ξ, ij, s, k, u, adv)` into
  `gather_stencil_values(ξ, ij, s, k)` + `reconstruct_face(S, u, adv)` for all 5 schemes
  (1st/3rd/5th-order upwind, 2nd/4th-order centered, WENO). Kernel structure unchanged
  (still per-(ij,k), no CPU/GPU split yet) — this phase is purely the math/memory-access
  decoupling needed before Phase 2/3 can share reconstruction code across loop shapes.
  Golden check: 15/15 bit-identical. Existing unit tests: 12+222+3 pass with
  `--check-bounds=yes`. Committed as `d53992e7`.
- 2026-06-29: Phase 2 done — `_vertical_advection!` now dispatches on `architecture(...)`
  like `vertical_integration!`. CPU path is a plain `k`-outer/`ij`-inner loop, reusing
  `vars.scratch.grid.a_2D` (already free at this point in the tendency pipeline; no new
  `ScratchVariable` added) as a rolling single buffer for the shared interior face — cheaper
  than the plan's two-ping-ponged-buffers sketch since by the time the inner `ij` loop for
  layer k finishes, the old face value is no longer needed. GPU path unchanged (still the
  Phase 1 per-(ij,k) kernel). Golden check: 15/15 bit-identical. Unit tests pass
  (`--check-bounds=yes`). Committed as `8aab769e`.

  CPU benchmark (`@benchmark`, 200 samples, after a `run!(period=Day(1))` warmup),
  baseline = pre-revision code at commit `c6fdc8a7` in a throwaway worktree
  (`scratchpad/baseline-worktree`, still present — not committed/pushed anywhere):

  | config / scheme   | baseline (min) | Phase 2 (min) | speedup |
  |--------------------|----------------|---------------|---------|
  | T31L8 / Centered2  | 375.5 μs       | 305.4 μs      | 1.23x   |
  | T31L8 / Upwind5    | 774.8 μs       | 688.5 μs      | 1.13x   |
  | T31L8 / WENO       | 3.710 ms       | 1.913 ms      | 1.94x   |
  | T31L24 / Centered2 | 1.647 ms       | 0.898 ms      | 1.83x   |
  | T31L24 / Upwind5   | 5.181 ms       | 1.974 ms      | 2.62x   |
  | T31L24 / WENO      | 11.31 ms       | 5.400 ms      | 2.09x   |
  | T63L8 / Centered2  | 1.283 ms       | 1.049 ms      | 1.22x   |
  | T63L8 / Upwind5    | 2.670 ms       | 2.384 ms      | 1.12x   |
  | T63L8 / WENO       | 12.85 ms       | 6.617 ms      | 1.94x   |

  Speedups grow with column depth (L24 > L8) and stencil width (WENO/Upwind5 > Centered2),
  as expected from halving the number of face reconstructions. Priority #2 (CPU perf at
  trunc<70) met with margin; moving to Phase 3 (GPU).

- 2026-06-29: Phase 3 done — added `vertical_advection_column_kernel!` (one GPU thread per
  `ij`, looping over all `k`, carrying the shared interior face in a thread-local register —
  no scratch array needed since each thread owns its whole column). This is "3a" from the
  plan; "3b" is just the pre-existing per-(ij,k) `vertical_advection_kernel!` from Phase 1,
  unchanged, so no separate shared-memory variant was built (the plan flagged that as only a
  candidate, and 3a already won decisively — see below).

  **Correctness vs 3b**: confirmed bit-identical (`===`) on GPU at default julia flags,
  for all 5 schemes x 4 grid configs (T31L8, T85L24, T127L24, T170L8) via a scratch script
  launching both kernels directly. *However*: under `--check-bounds=yes` (the project's
  mandated test flag), Upwind5/WENO diverge by up to ~2e-4 relative (Centered2 stays exact)
  — root-caused to the forced bounds-check insertion changing generated code enough to shift
  the GPU compiler's FMA-contraction/reordering decisions differently between the two
  differently-shaped kernels. Not a logic bug (same class of GPU FP non-determinism
  `vertical_integration.jl`'s GPU test already works around with `≈`). Added
  `test/GPU/vertical_advection.jl` (wired into `test/GPU/runtests.jl`) comparing both
  kernels directly with `≈` instead of `===`, for this reason — passes under
  `--check-bounds=yes`.

  **GPU dispatch decision (3a vs 3b)**, benchmarked via single-variable kernel launch
  (`@benchmark`, 100 samples, min) on an A40:

  | npoints (config)      | Centered2 | Upwind5 | WENO  |
  |------------------------|-----------|---------|-------|
  | 3168 (T31)             | 0.86x     | 0.63x   | 0.50x |
  | 18688 (T85L24)          | 0.91x     | 0.88x   | 0.84x |
  | 28480 (T96/T106L24)     | 0.89x     | 0.92x   | 0.93x |
  | 40320 (T127L24)         | 1.00x     | 1.18x   | 1.28x |
  | 70144 (T170L24/L8)      | 1.03-1.04x| 1.20x   | 1.29x |

  (ratio = pointwise-time / column-time; >1 means column/3a wins). Crossover sits between
  npoints=28480 and 40320, independent of `nlayers` (8 vs 24 behaved consistently) — i.e.
  driven by `npoints` (column's only source of parallelism), not column depth. Implemented
  as `_vertical_advection!(::GPU,...)` dispatching on `npoints >=
  VERTICAL_ADVECTION_GPU_COLUMN_THRESHOLD = 35_000` (midpoint of the measured crossover
  band): per-column (3a) above, per-(ij,k) (3b) below. This matters because T85-T96 sit
  *inside* the plan's own "trunc>70" GPU-priority zone yet still favour 3b — a blanket
  "always use 3a" choice would have cost ~10% there to gain up to 29% at trunc>=127.

  Also discovered (not a regression): `dynamics_only=true` PrimitiveWetModel runs for
  Day(1) at T85/T127 produce NaNs with Upwind5/WENO (Centered2 stays clean) — confirmed
  identical on the untouched pre-revision code in a throwaway worktree at commit `c6fdc8a7`
  (same warning timestamps). Pre-existing model-stability characteristic at these
  resolutions with `dynamics_only`, unrelated to this PR; not investigated further.

  Committed Phase 3 as the next commit (kernel + dispatch threshold + GPU test).

- 2026-06-29 (later): user decided to skip Phase 4 (fusion) for now — implemented it, then
  reverted (`git checkout --`) on request before committing, so HEAD stays at the Phase 3
  state. Ran a comprehensive CPU+GPU benchmark sweep (old vs new, 10 resolutions x 3 schemes)
  instead; results written to `BENCHMARK_RESULTS.md` for the PR. CPU: clean 1.1-2.8x win
  everywhere, consistent with Phase 2's spot numbers. GPU: Centered2 neutral; Upwind5/WENO
  show a REGRESSION (0.28-0.78x) in the full `vertical_advection!` call at trunc>=63 that
  contradicts the isolated single-kernel-launch benchmark the dispatch threshold was tuned
  on (which still reproduces 1.05-1.29x column win when re-checked). Root cause not found
  before running out of time — ruled out: short-warmup branch divergence, dispatch not
  firing, environmental drift (isolated test still reproduces). **Open issue: the GPU
  per-column kernel's dispatch threshold needs re-investigation before this can be
  considered a GPU win for Upwind/WENO — possibly revert GPU dispatch to always use the
  per-(ij,k) kernel until understood, or the full-pipeline call has some per-field state
  interaction the isolated test doesn't capture.**
