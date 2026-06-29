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

- Interactive shell is `login03` (HPC, Lmod + SLURM). `nvidia-smi` shows an A40 GPU visible
  from the login node, but the user's own `submit_bench.sh` / `submit_run_benchmarks.sh`
  (untracked, already in the repo root) submit GPU work via
  `sbatch --partition=gpu --gres=gpu:1 --qos=gpushort --account=flai` — so Phase 3 GPU runs
  will go through `sbatch`, not direct execution on the login node.
- Julia via `module load julia` (default 1.12.2). Each Bash tool call is a fresh shell, so
  `module load julia` must be chained with `&&` in every command that calls `julia`.

## Checklist

- [x] Phase 0: golden capture + face-consistency assertion (scratch validation script, no src change)
- [ ] Phase 1: value-based `reconstruct_face` refactor (bit-identical) + tests
- [ ] Phase 2: CPU plain-loop flux-form path + golden + CPU bench
- [ ] Phase 3: GPU kernel(s) 3a/3b + golden(GPU) + GPU bench (decide 3a vs 3b)
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
  Capture run against the pre-revision code (HEAD, before any Phase 1+ edits) kicked off
  in background; confirming completion before starting Phase 1 edits.
