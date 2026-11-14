# Implicit solver arrays as variables, implicit solvers reduced to options

> Status: **planned**. Plan only, implementation to follow (depends on #1402 being merged).

Date of initial draft: 2026-11-14

Base revision: caed82d6 (main), to be implemented on top of
[#1402](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1402)

## Originating prompt

> Idea: Currently the implicit arrays live in model.time_stepping.implicit what about model them
> to the variables? That way implicit reduces to its options and arrays can be allocated at
> initialisation time depending on the model.

> Yes please write a follow-up PR plan and push this in a separate draft PR with just the plan in
> it and a summary description. That way another session can easily pick that up

## Revision log

- Initial draft.

## Problem description

With #1402 the implicit solver is part of the time stepper, `model.time_stepping.implicit`
(`DefaultImplicit()` is replaced by the model's default on construction). `ImplicitShallowWater`
only holds an option (`centering`), but `ImplicitPrimitiveEquation` mixes options and state:

| Kind | Fields of `ImplicitPrimitiveEquation` |
|---|---|
| Options | `centering`, `reference_temperature`, `schedule` |
| Dimensions | `truncation`, `nlayers`, `nlayers_val::Val{NLAYERS}` |
| State | `Δt::Ref` (time step the operators were built with), `temp_profile` (nlayers), `R`, `L`, `L1`, `L3`, `S` (nlayers × nlayers), `U`, `W`, `L0`, `L2`, `L4` (nlayers), `S⁻¹` ((truncation+1) × nlayers × nlayers) |

Consequences:

1. The implicit solver needs a `SpectralGrid` to be constructed and carries 9 type parameters
   (array types, `NLAYERS`, the `Ref` type for Reactant).
2. Differentiability: Enzyme tests pass `Const(model)`. With `GlobalMeanReferenceTemperature` (or
   `MaximumReferenceTemperature`) `temp_profile` and all operators depend on the model state, so
   treating them as constant is incorrect for gradients through `reinitialize!`.
3. Reactant: model fields are baked in as constants; `implicit.Δt` therefore needs `scalar_ref`
   to be traceable, and `reinitialize!` has a separate Reactant method
   (`ext/SpeedyWeatherReactantExt.jl`) that writes into model arrays.
4. Inconsistent with the design direction of the `Variables` refactor: the model holds
   parameters, `Variables` hold all mutable state.

## Background

Relevant code (on the #1402 branch):

- `SpeedyWeather/src/time_stepping/implicit/implicit_primitive_equations.jl`: struct (~l. 60-135),
  generator `ImplicitPrimitiveEquation(spectral_grid; ...)` (~l. 136), `default_implicit`,
  `initialize!(implicit, model)` (schedule + `warn_uniform_spacing`), `reinitialize!(implicit,
  model, vars)` (decides on `update` from `implicit.Δt[]` and the schedule) and
  `reinitialize!(implicit, update::Bool, Δt, model, vars)` (sets `temp_profile` from
  `reference_temperature(...)`), `initialize!(implicit, Δt, temp_profile, geometry, geopotential,
  atmosphere, adiabatic_conversion)` (builds `R, U, L, W, L0-L4`, then `invert_S!`),
  `invert_S_lu!` (CPU, via `on_architecture`), `invert_S_gauss_jordan!` + kernel
  `implicit_invert_S_kernel!` (device/Reactant), `implicit_correction!` + kernel
  `implicit_primitive_leapfrog_kernel!` (reads `S⁻¹, R, U, L, W, nlayers_val`).
- Readers of `temp_profile` in `SpeedyWeather/src/dynamics/tendencies.jl`: temperature anomaly
  (~l. 148, 169), `linear_pressure_gradient!` (~l. 486), `vordiv_tendencies!` (~l. 541),
  temperature tendency (~l. 690). These receive `implicit` as an argument (function barriers).
- `ext/SpeedyWeatherReactantExt.jl`: `reinitialize!(implicit, update::TracedRNumber{Bool}, ...)`.
- Variables infrastructure: `SpeedyWeather/src/variables/dimensions.jl` has `ScalarDim`
  (allocates a `Ref`), `VectorDim`, `MatrixDim(m, n)`; no 3D dimension yet.
  `variables(time_stepping, model)` (in `steppers/general.jl`) already forwards to
  `variables(implicit_solver(time_stepping), model)`.
- Tests touching these arrays: `test/dynamics/implicit.jl` (S⁻¹ via LU vs Gauss-Jordan, operator
  `L` vs linearised tendency), `test/dynamics/time_stepping.jl` (reference temperature profiles,
  ~l. 370-405), `test/dynamics/vertical_advection.jl` (~l. 77, modifies `temp_profile`),
  `test/reactant/test_implicit.jl` (sets `temp_profile`, `implicit.Δt[]`),
  `test/differentiability/primitivewet.jl` (`Const(model.time_stepping.implicit)`).

## Summary of changes

1. **New variable dimension** in `variables/dimensions.jl` for a matrix per spherical harmonic
   degree, e.g. `DegreeMatrixDim(n)` allocating `zeros(NF, truncation + 1, n, n)` on the model's
   architecture (name to be decided; alternatively a general `TensorDim(dims...)`). Plus tests in
   the variables tests.
2. **`variables(implicit::ImplicitPrimitiveEquation, model)`** declares, in a namespace
   `:implicit` of `vars.dynamics` (`DynamicsVariable`, or a dedicated group if preferred):
   - `Δt` as `ScalarDim()` (time step the operators were built with),
   - `temp_profile`, `U`, `W`, `L0`, `L2`, `L4` as `VectorDim(nlayers)`,
   - `R`, `L`, `L1`, `L3`, `S` as `MatrixDim(nlayers, nlayers)` (`S` possibly a `ScratchVariable`),
   - `S⁻¹` with the new dimension.
   `nlayers` from `model.geometry` / `model.spectral_grid`.
3. **`ImplicitPrimitiveEquation` reduced to options**: `centering`, `reference_temperature`,
   `schedule`; type parameters `NF` (or none, converting at use), schedule and reference
   temperature types. Generator `ImplicitPrimitiveEquation(spectral_grid; kwargs...)` kept for
   the number format (and backwards compatibility of the call), `ImplicitPrimitiveEquation(; kwargs...)`
   added.
4. **Functions take the arrays from `vars`**: `initialize!(implicit, Δt, ...)` →
   `initialize!(vars.dynamics.implicit, implicit, Δt, ...)` (or pass the NamedTuple of arrays as
   `implicit_vars`), likewise `invert_S_lu!`, `invert_S_gauss_jordan!`, `reinitialize!` (and its
   Reactant method), `implicit_correction!` and its kernel. Keep function barriers: unpack
   `vars.dynamics.implicit` once and pass arrays to kernels.
5. **`temp_profile` readers** in `tendencies.jl` read `vars.dynamics.implicit.temp_profile` (they
   already receive `vars`).
6. **`nlayers_val`**: kernels currently use `Val(nlayers)` from a type parameter for type-stable
   vertical loops. Options: derive `Val(nlayers)` at the call site from the size of `S⁻¹`
   (`Val(size(S⁻¹, 2))`, a dynamic `Val` but behind a function barrier, so one dynamic dispatch
   per call) or from `model.geometry` if that carries `nlayers` as a type parameter. Benchmark
   both against main (see Testing).
7. **Geometry-only operators**: `R` and `L0`-`L4` depend only on geometry/atmosphere, not on Δt or
   the reference temperature. Optionally compute them once in `initialize!` (still variables), not
   on every `reinitialize!`. Not required for this PR.
8. **ShallowWater**: `ImplicitShallowWater` already holds options only; no variables needed. Check
   that `default_implicit` and `DefaultImplicit` (from #1402) no longer need the `SpectralGrid` for
   any implicit solver and simplify if so.
9. **Reactant**: `implicit.Δt` becomes the scalar variable (no `scalar_ref`), the Reactant method
   of `reinitialize!` writes into variables; remove what becomes unnecessary in the extension.
10. **Breaking**: `model.time_stepping.implicit.temp_profile` (etc.) → `vars.dynamics.implicit.temp_profile`.
    CHANGELOG line, docs (`docs/src/primitiveequation.md` implicit section, `time_integration.md`
    implicit solver section from #1402), SpeedyWeather version: already a `-DEV` minor bump.

## Testing and verification

- Existing implicit tests adapted to read arrays from `vars`: `test/dynamics/implicit.jl`,
  `test/dynamics/time_stepping.jl` (reference temperature profiles),
  `test/dynamics/vertical_advection.jl`, `test/reactant/test_implicit.jl`.
- Bit-identical results: run PrimitiveWet/Dry (e.g. T32L8, 10 days, `GlobalMeanReferenceTemperature`
  and the default isothermal) on main+#1402 and on this branch, compare prognostic variables
  (should be identical, same operations in the same order).
- Differentiability: `test/differentiability/primitivewet.jl` with `vars` now carrying the implicit
  arrays (`Duplicated`); add a check that the gradient w.r.t. temperature through `reinitialize!`
  with `GlobalMeanReferenceTemperature` is non-zero (it was implicitly zero before).
- GPU (CUDA) and Reactant test suites, as the arrays move between model and variables.
- Benchmark `manual_benchmarking.jl --debug` vs main (kernel type stability via `nlayers_val`).

## Documentation changes

- Primitive equation implicit section: where the operators live and that they are (re)computed on
  `reinitialize!`.
- `time_integration.md` "Implicit solver" section: implicit solvers are options, their state is
  in the variables.
- Variables docs: new dimension.

## Known limitations

- Memory of `S⁻¹` ((truncation+1) × nlayers²) moves into `Variables` and is therefore duplicated
  in shadow variables for Enzyme (`dvars`), as it should be for correct gradients.

## Future work

- Implicit solver for `NCycleLorenz` with the primitive equations (not implemented, see #917).
- Compute geometry-only operators once (item 7) if not done here.
