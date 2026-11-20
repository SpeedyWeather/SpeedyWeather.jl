# Precompile workload for the default model

> Status: **completed**. A PrecompileTools workload for the default `PrimitiveWetModel` and `nlayers` removed from the type of `ImplicitPrimitiveEquation`. In a fresh session `initialize!` takes 0.5 s instead of 54 s.

Date of initial draft: 2026-11-19

Base revision: 884416c7 (main)

## Originating prompt

> Have a look at PR1212. The model initialisation of the default model is very slow. Merge main and reevaluate what can be done to make it faster

> This PR is now massively big, can you turn this PR into just doing the precompile targeting main? Then we can split the compile time issue into 2, offloading to precompile and generally getting compile time faster regardless when its done

## Revision log

- 2026-11-19: developed on top of #1212 (`set!` type erasure, variable collection as a vector),
  see `docs/dev/2026-08/reduce-initialize-compile-time.md` on that branch.
- 2026-11-20: split out of #1212 so that precompilation and general compile-time reductions
  are reviewed separately. Re-measured against plain main.

## Problem description

In a fresh session, constructing, initializing and running the default `PrimitiveWetModel` is
almost entirely compilation: ~9 s for the model constructor, ~54 s for `initialize!(model)`
and ~40 s for the first `run!`, while the actual work is well below a second. SpeedyWeather
has no precompile workload, so all of this is paid in every new session.

## Background

#1212 profiled `initialize!` with SnoopCompile. After its two structural fixes the remaining
cost is flat: NetCDF reading and interpolation of boundary conditions, the spectral transform,
KernelAbstractions launch code, `Variables` construction, initial conditions. Whatever step
first touches one of these pays for it. This is the case for precompilation: compile once at
`Pkg.precompile` time, cached in the package image.

While testing the workload, changing `truncation` reused the cache but **any `nlayers` other
than 8 recompiled everything**. `ImplicitPrimitiveEquation` stored `nlayers_val::Val{NLAYERS}`
(from #1163, so that the implicit kernel can hold the vertical geopotential in an `ntuple` of
static length), which made `nlayers` a type parameter of the model.

## Summary of changes

1. `SpeedyWeather/src/precompile.jl`: a `@setup_workload`/`@compile_workload` that builds the
   default `PrimitiveWetModel` with `NetCDFOutput` to a temporary directory, `initialize!`s it
   and runs two time steps with output. It is wrapped in `try`/`catch`: the default model
   downloads boundary-condition assets on first use, and a failed download (no network at
   precompile time) skips the workload with a warning instead of failing precompilation.
   It can be disabled with the PrecompileTools preference `precompile_workload = false`.
2. `ImplicitPrimitiveEquation` loses its `NLAYERS` type parameter and `nlayers_val` field.
   `Val(nlayers)` is created at the kernel launch in `implicit_correction!`, so only the
   implicit kernel specializes on `nlayers`. Cost: one dynamic dispatch per time step.
3. `PrecompileTools` added as a dependency of SpeedyWeather.

## Testing and verification

Fresh Julia session each, default `SpectralGrid()` (T31, Float32, CPU), Linux x86_64,
Julia 1.13.4:

| | main | this PR |
|---|---|---|
| `using SpeedyWeather` | 3.6 s | 3.3 s |
| `PrimitiveWetModel(spectral_grid)` | 9.0 s | 0.1 s |
| `initialize!(model)` | 52–55 s | **0.5 s** |
| first `run!(simulation, period = Day(1))` | 39.7 s | **0.6 s** |
| first `run!(..., output = true)` | 59.9 s | **0.8 s** |
| `nlayers = 16`: `initialize!` / first `run!` | 52.1 s / 39.5 s | 0.6 s / 1.1 s |
| precompiling SpeedyWeather | 11.7 s | **88.7 s** |
| package image | 12 MB | 101 MB |

Without change 2, `nlayers = 10` cost 10.2 s / 19.0 s even with the workload.

Existing tests cover the implicit solver for several `nlayers` (`dynamics/implicit.jl`,
`time_stepping.jl`, `vertical_coordinates.jl` with 8, 16, 24 layers). Enzyme and Reactant
compatibility of the runtime `Val(nlayers)` is checked by their CI.

## Documentation changes

None, the workload is transparent to users. The preference to disable it is documented in
`precompile.jl`.

## Known limitations

- Only the default `PrimitiveWetModel` with `NF = Float32` on CPU is covered. Measured on
  the #1212 branch with the workload, still compiling on first use: `PrimitiveDryModel` ~28 s,
  `ShallowWaterModel` ~15 s, `BarotropicModel` ~8 s, `NF = Float64` ~57 s. Each could be added
  to the workload at roughly that cost in precompile time.
- Precompiling SpeedyWeather takes ~89 s instead of ~12 s, image ~101 MB instead of 12 MB.
- If the asset download fails at precompile time, the workload is skipped and the first
  `initialize!` is slow again until the next precompilation.

## Future work

- General compile-time reductions independent of precompilation (#1212): type-erased `set!`,
  variable definitions collected in a vector. These also shrink the workload's precompile time.
- Decide whether to cover further model types or `Float64` in the workload.
