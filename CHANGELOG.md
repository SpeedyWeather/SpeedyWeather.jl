# Changelog

## Unreleased

- set! for time step [#650](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/650)
- Bug: ZonalWind had sqrt of negative for uncommon resolutions [#649](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/649)
- Spectral Gradients are now differentiable, with correctness check in extended CI [#638](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/638)
- Tracer advection [#579](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/579)
- Buildkite CI with dummy pipeline [#646](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/646)
- ConvectiveHeating implemented [#639](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/639)
- Number format flexibility with set! [#634](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/634)
- Forcing/drag for primitive models [#635](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/635)
- RingGrids indexing with leading Colon should now always return another RingGrid instance [#637](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/637)
- Roll back GPUArrays upgrade to ensure CUDA compatibility [#636](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/636)
- Change default timestep to 40min at T31 [#623](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/623)

## v0.13

- AbstractSurfacePerturbation introduced [#631](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/631)
- AbstractStochasticPhysics and SPPT implemented [#629](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/629)
- Introduced seperate extended tests that are not run on every commit [#628](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/628)
- One-band longwave radiation [#624](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/624)
- compat entry for FiniteDifferences.jl [#620](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/620)
- Slab ocean and net surface fluxes [#613](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/613)
- added compat entries for FiniteDifferences.jl and Enzyme.jl [#620](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/620) [#622](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/622)
- update to GPUArrays v11, JLArrays v0.2 and remove CUDA from tests [#590](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/590)
- Power spectrum for n-dim LowerTriangularArrays [#618](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/618)
- Added custom EnzymeRules for the SpectralTransform and an extension for compatibility with FiniteDifferences.jl [#589](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/589) [#625](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/625) [#627](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/627)
- Optical depth introduced and array-agnostic ColumnVariables [#606](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/606)
- Include large-scale condensation tests [#615](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/615)
- bugfix: large-scale condensation also at <100% [#609](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/609)

## v0.12.1

- ConstantLandTemperature implemented [#612](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/612)
- set! for more boundary conditions [#611](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/611)
- SpectralFilter for horizontal diffusion [#601](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/601)
- GeoMakie weak dependency, globe function for 3D data visualisation [#600](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/600)
- Zonal mean for AbstractGridArray [#603](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/603)
- Rossby-Haurwitz wave with initial conditions for interface displacement for shallow water models[#604](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/604)

## v0.12.0

- OctaminimalGaussianArray/Grid to start with 4 points around the poles [#595](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/595)
- Output both accumulated and precipitation rate as netCDF [#596](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/596)
- Random processes for random pattern generation [#592](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/592)
- Also allow SpectralGrid as positional argument to model constructors [#593](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/593)
- De-interweave SpectralTransform [#587](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/587)
- Rossby-Haurwitz wave initial conditions [#591](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/591)
- Haversine formula and AbstractSphericalDistance [#588](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/588)
- Array-agnostic SpectralTransform [#583](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/583)
- Move CUDA dependency into extension [#586](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/586)
- Stop supporting Julia v1.9 [#585](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/585)
- `feedback.verbose` (de/activate the progressbar) is now set to `isinteractive()` to disable automatically for documentation [#582](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/582)

## v0.11.0

- Extend `set!` with `orography` keyword argument to set the orography with `set!` [#578](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/578)
- Added a new benchmark suite for the dynamics functions [#577](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/577)
- Introduced a new `set!` function that allows to set `PrognosticVariables` to new values with keyword arguments [#563](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/563)
- Restructured dynamical core with prognostic/diagnostic variables array-agnostic and 3-dimensional [#525](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/525)
- Modularised NetCDF output [#573](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/573)
- Fixed a bug in RingGrids, now broadcasts are defined even when the dimensions mismatch in some cases [#568](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/568)
- RingGrids: To wrap an Array with the horizontal dimension in matrix shape into a full grid, one has to use e.g. `FullGaussianGrid(map, input_as=Matrix)` now. [#572](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/572)
- CompatHelper: Allow for JLD2.jl v0.5
