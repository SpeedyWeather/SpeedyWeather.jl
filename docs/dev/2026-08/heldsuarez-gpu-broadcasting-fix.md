# Fix GPU scalar indexing when a `Field` broadcasts against a bare GPU array

> Status: **completed**. `RingGrids/src/field.jl` gained explicit `BroadcastStyle` combine
> rules for `FieldGPUStyle` vs bare GPU arrays and vs itself at a different dimensionality;
> covered by `SpeedyWeather/test/GPU/broadcasting.jl`; merged into
> [PR #1219](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1219).

Date of initial draft: 2026-08-22 (written retroactively, after the fix; see revision log)

Base revision: `557b38d5` (`main`, branch `mc/heldgpu`)

## Originating prompt

> Can you run the Held Saurez example from the docs locally but on the GPU?

followed, after the run failed, by:

> Can you create a branch called mc/heldgpu to fix this? First test whether leaving the field
> wrapper with .data works which is a reasonable work around but have a look into broadcasting
> too. I don't want you to rewrite a lot of it as it affects many other parts of the code but
> if you find a fix there without many consequences then implement that one.

## Revision log

- **2026-08-22, initial draft.** Written after the fix was already implemented, tested, and
  opened as PR #1219, at the user's request to check whether a dev plan existed for this
  change (`CLAUDE.md` requires one for substantial changes) and to add one retroactively since
  it did not. This document reconstructs the investigation and decisions made during that
  session; no new code changes accompany it.

## Problem description

Running the docs' Held-Suarez example (`docs/src/examples_3D.md`) with `architecture =
SpeedyWeather.GPU()` crashed one timestep in in `SpeedyWeather/src/dynamics/drag.jl:62`:

```julia
@. Fu -= drag.drag_coefs' .* u
```

with

```
ERROR: Scalar indexing is disallowed.
Invocation of getindex resulted in scalar indexing of a GPU array.
```

`drag.drag_coefs'` is an `Adjoint` of a plain `CuVector`; `u`/`Fu` are RingGrids `Field`s
(GPU-backed, via a `SubArray` of a `CuArray`). Both operands are, individually, entirely
GPU-broadcast-compatible; the model runs equivalent CPU code without issue.

## Background

RingGrids' `Field` defines its own `Base.Broadcast.BroadcastStyle`, `FieldStyle{N,Grid}` for
CPU and `FieldGPUStyle{N,Grid} <: GPUArraysCore.AbstractGPUArrayStyle{N}` for GPU-backed data
(`RingGrids/src/field.jl:636-698`, before this change). Both styles carry `Grid` as a type
parameter, so that broadcasting refuses to mix fields on incompatible grids, and both define
the handful of `(::Val{M})`-style dimension-upgrade constructors that
`Base.Broadcast.BroadcastStyle` combine logic calls when it needs to resolve two array styles
of different dimensionality to a common one.

On CPU, `drag.drag_coefs'` gets `Broadcast.DefaultArrayStyle{2}()` (Julia's default style for
any plain `Array`), and Base's broadcast machinery has a dedicated, general rule for combining
*any* `AbstractArrayStyle` with `DefaultArrayStyle`:

```julia
# base/broadcast.jl
BroadcastStyle(a::AbstractArrayStyle{M}, ::DefaultArrayStyle{N}) where {M,N} =
    typeof(a)(Val(max(M, N)))
```

so `FieldStyle` always wins over a bare `Array`/`Adjoint` without RingGrids having to do
anything. On GPU, however, `drag.drag_coefs'` gets `CUDA.CuArrayStyle{1}` — not
`DefaultArrayStyle` — so that special case does not apply. The only other generic rule Base
provides for combining two `AbstractArrayStyle`s is

```julia
function BroadcastStyle(a::A, b::B) where {A<:AbstractArrayStyle{M},B<:AbstractArrayStyle{N}} where {M,N}
    if Base.typename(A) === Base.typename(B)
        return A(Val(max(M, N)))
    end
    return Unknown()
end
```

which requires the two styles to share a type *name* — `FieldGPUStyle !== CuArrayStyle`, so
this returns `Broadcast.Unknown()`. `Unknown()` broadcasts fall back to a generic,
non-GPU-aware materialization path, which is what performs scalar `getindex` calls on the
`CuArray` and trips GPUArrays' scalar-indexing guard. This is a broadcasting-dispatch gap, not
a bug in `LinearDrag` or `drag.jl` — any code broadcasting a `Field` directly against a bare
GPU array (rather than another `Field`) would hit the same failure. It went unnoticed because
the existing GPU test suite does not exercise this combination.

## Summary of changes

Investigated in order, per the user's request:

1. **Confirmed the `.data` workaround works.** `@. Fu.data -= drag.drag_coefs' .* u.data`
   broadcasts plain `CuArray`s on both sides and succeeds — but would need to be applied at
   every call site with this pattern, and silently reintroduces the same failure mode anywhere
   future code broadcasts a `Field` against a bare GPU array.

2. **Traced the root cause** using a small reproduction script that constructed the actual
   `Variables`/`LinearDrag` state from an initialized GPU model and printed
   `Base.BroadcastStyle(...)` for each operand (see Problem description / Background above).

3. **Implemented the general fix** in `RingGrids/src/field.jl` (~10 lines, no other files
   touched): added

   ```julia
   Base.BroadcastStyle(::FieldGPUStyle{N, Grid}, ::GPUArrays.AbstractGPUArrayStyle{M}) where {N, Grid, M} =
       FieldGPUStyle{N, Grid}(Val(max(N, M)))
   Base.BroadcastStyle(::GPUArrays.AbstractGPUArrayStyle{M}, ::FieldGPUStyle{N, Grid}) where {N, Grid, M} =
       FieldGPUStyle{N, Grid}(Val(max(N, M)))
   ```

   mirroring Base's `DefaultArrayStyle` special case, so a `Field` always wins over a bare GPU
   array on GPU too, at whichever dimensionality is larger (matching the existing CPU
   behavior). This alone introduced a new *ambiguity*: `FieldGPUStyle{N,Grid}` is itself an
   `AbstractGPUArrayStyle{N}`, so combining two `FieldGPUStyle`s of different `N` (e.g. a 3D
   wind field against a 2D field on the same grid, hit inside `pressure_gradient_flux!`) now
   matched both the new rule and Base's generic same-typename rule. Julia's own `MethodError`
   for the ambiguity named the exact fix needed, which was added as a third method:

   ```julia
   Base.BroadcastStyle(::FieldGPUStyle{M, Grid}, ::FieldGPUStyle{N, Grid}) where {M, N, Grid} =
       FieldGPUStyle{M, Grid}(Val(max(M, N)))
   ```

   `drag.jl` itself is unchanged — the original, unmodified `@. Fu -= drag.drag_coefs' .* u`
   now works.

4. **Added a regression test**, `SpeedyWeather/test/GPU/broadcasting.jl`, covering three cases
   on GPU: `Field .* Adjoint(vector)`, `Field .* reshaped vector`, and `Field .* Field` with
   mismatched `N` on the same grid. Wired into `SpeedyWeather/test/GPU/runtests.jl`.

5. **Bumped `RingGrids`** from `0.1.9` to `0.1.9+DEV` (internal fix, non-breaking) and added the
   `CHANGELOG.md` entry, per repository convention.

## Testing and verification

- The unmodified Held-Suarez docs example (`PrimitiveDryModel` + `HeldSuarez` forcing +
  `LinearDrag` + `EarthOrography`, truncation=32, nlayers=8) run to completion on GPU for
  `period = Day(20)`.
- New test file:

  ```bash
  CUDA_VISIBLE_DEVICES=<idle GPU> julia --project=SpeedyWeather/test/GPU/CUDA \
    SpeedyWeather/test/GPU/runtests.jl
  ```

  (run standalone against just `broadcasting.jl` during development: 3/3 passed).
- `RingGrids/test/field_types.jl` (CPU) — no regressions (11/11 passed).
- The full GPU and RingGrids test suites were **not** run end-to-end for this change (the user
  asked not to, citing runtime); this is a known gap, see below.

## Documentation changes

None beyond this plan and the `CHANGELOG.md` entry. The docs' Held-Suarez example
(`docs/src/examples_3D.md`) is unchanged — it already showed correct, idiomatic code; the bug
was in the library, not the example.

## Known limitations

- The full `RingGrids` and GPU (`SpeedyWeather/test/GPU`) test suites were not run against this
  change; only the new targeted test and the pre-existing `field_types.jl` CPU test were
  verified. It is possible, though unlikely given the previous behavior was an outright crash,
  that some other code path relied on the old `Unknown()` fallback behavior.
- The fix only covers `FieldGPUStyle` combined with `AbstractGPUArrayStyle`. If another custom
  `AbstractArrayStyle` is introduced elsewhere in the codebase that needs to broadcast against
  a `Field` on GPU, it will need the same treatment.

## Future work

- Consider whether `Test.detect_ambiguities` (or an equivalent CI check restricted to
  `RingGrids`) should be run to catch this class of `BroadcastStyle` ambiguity earlier than a
  runtime `MethodError`.
- If more GPU broadcast styles are added in the future, consider generalizing the two-argument
  combine rules above into a single macro or helper rather than hand-writing each pair, to keep
  the growth in cases manageable.
