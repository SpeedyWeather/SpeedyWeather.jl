# `Parameterization` function stub for parameterizations defined in other packages

> Status: **in progress**. Implemented, tested locally and opened as
> [#1261](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1261); awaiting review.

Date of initial draft: 2026-09-19

Base revision: `35d764b6` (`main`)

## Originating prompt

> Could you have a look at PR1089 and review it? Max made a comment about introducing a
> Parameterization stub, can you evaluate this and then open a new pull request for it, if
> deemed a good idea?

Max's comment on #1089 (`docs/src/co2_forcing.md`):

> What do you think about just defining a function stub `Parameterization` in SpeedyWeather.
> And then all external libraries that define parameterisations can use that to extend it and
> we avoid this annoying no-exports-from-extensions problem like that? Basically:
>
> ```julia
> external_param = CustomParam(....) # defined in main code of package like AnalyticRadiation without SpeedyWeather dependency
> param = SpeedyWeather.Parameterization(spectral_grid, external_param; kwargs...) # load into SpeedyWeather, extended in extension.
> ```

## Revision log

- **2026-09-19, implementation.** Implemented directly after evaluating the proposal. This
  plan was written retroactively and was not signed off before implementation.
- **2026-09-19, "I think a docs/dev plan is missing, can you add it?"** Added this document.

## Problem description

Julia package extensions cannot export new names. A package that defines a parameterization
without depending on SpeedyWeather has to put its SpeedyWeather wrapper (a subtype of
`AbstractParameterization`) in an extension. Users can then only reach that wrapper through
`Base.get_extension`. #1089 documents exactly this for AnalyticBandRadiation.jl (now
NumericalRadiation.jl):

```julia
const SpeedyExt = Base.get_extension(AnalyticBandRadiation, :AnalyticBandRadiationSpeedyWeatherExt)
longwave = SpeedyExt.SpeedyAnalyticBandLongwave(spectral_grid)
```

This is verbose. It also exposes the extension's module name, which has already changed with the
package rename, so the documented snippet is out of date.

## Background

- SpeedyWeather already defines exported function stubs whose methods live in extensions:
  `animate` (GeoMakie extension) and `TerrariumOutput` (Terrarium extension).
- An extension that adds a method to a SpeedyWeather function, dispatching on a type owned by
  its parent package, is not type piracy.
- Parameterizations are conventionally built with a generator taking the `SpectralGrid` as the
  first positional argument and options as keyword arguments (see
  `docs/src/parameterizations.md`, "Define the generator function").

## Summary of changes

- `SpeedyWeather/src/parameterizations/general.jl`: `function Parameterization end` with a
  docstring that explains the pattern, plus one fallback
  `Parameterization(::SpectralGrid, p::AbstractParameterization) = p`. Generic code can
  therefore call it on SpeedyWeather parameterizations too. Anything else has no method and
  throws a `MethodError`.
- `SpeedyWeather/src/SpeedyWeather.jl`: `export Parameterization`, next to `TerrariumOutput`.

Intended use in an external package's extension:

```julia
SpeedyWeather.Parameterization(spectral_grid::SpectralGrid, scheme::AnalyticBandLongwave; kwargs...) =
    SpeedyAnalyticBandLongwave(spectral_grid, scheme; kwargs...)
```

and by users:

```julia
longwave_radiation = Parameterization(spectral_grid, AnalyticBandLongwave())
model = PrimitiveWetModel(spectral_grid; longwave_radiation)
```

## Testing and verification

`SpeedyWeather/test/parameterizations/custom_parameterization.jl` gains the testset
"Parameterization from external packages". It mocks an external type `ExternalCooling` (not
a subtype of anything in SpeedyWeather) and adds a `Parameterization` method that maps it
onto `UniformCooling`. The testset checks:

- the returned type, number format and fields, with keyword arguments forwarded
- the pass-through for an `AbstractParameterization`
- `MethodError` for unsupported inputs
- that a `PrimitiveDryModel` with this parameterization initializes and runs 3 time steps

Run locally on Julia 1.13: all 6 new tests and the 6337 existing custom-parameterization tests
pass.

## Documentation changes

- New subsection "Parameterizations from other packages" in `docs/src/parameterizations.md`,
  under "Use your parameterization".
- `CHANGELOG.md` entry.

## Known limitations

- Users still have to call `Parameterization(spectral_grid, scheme)` explicitly. Passing
  `scheme` directly to a model constructor doesn't work yet.
- Nothing enforces that extension methods return an `AbstractParameterization`.
- Adoption requires a one-line change in NumericalRadiation.jl's SpeedyWeather extension, and
  #1089's docs need updating to use it.

## Future work

- Have model constructors call `Parameterization(spectral_grid, x)` on each parameterization
  keyword. Then `PrimitiveWetModel(spectral_grid; longwave_radiation = AnalyticBandLongwave())`
  would work directly.
- Add a `SpeedyWeather.Parameterization` method to NumericalRadiation.jl's SpeedyWeather
  extension, and switch #1089's examples from `Base.get_extension` to `Parameterization`.
