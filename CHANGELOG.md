# Changelog

## Unreleased

- Added documentation to output surface variables [#708](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/708)
- `get_gridcell_polygons` in RingGrids [#706](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/706)
- Revised extended CI tests for differentiability can be launched on comment in PRs [#695](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/695), [#696](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/696), [#697](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/697), [#698](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/698)
- Reverse for AbstractGridArrays [#691](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/691)
- LandGeometry, LandThermodynamics as component of LandModel, DryLandModel [#677](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/677)
- Diagnostic albedo separate for ocean/land [#677](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/677)
- Increasing precision for accumulated precipitation output [#685](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/685)
- Allow steps to be specified for for run!(simulation, steps=n) [#684](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/684)
- Restrict extended CI tests to Julia v1.10 (due to Enzyme) [#681](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/681)
- SigmaCoordinates constructors [#679](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/679)
- Full timestepping differentiable by Enzyme with Julia 1.10 [#656](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/656)
- Type instabilities etc flagged by Reactant [#674](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/674)
- Added JLD2 Output to directly output PrognosticVariables and DiagnosticVariables [#675](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/675)
- Compat for GPUArrays v11, JLArrays 0.2, KernelAbstractions 0.9 [#658](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/658)
- MITgcm's 2-layer land bucket model [#672](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/672)
- Fix grid_cell_average! bugs for some grids and resolutions [#673](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/673)
- Interfaces for interpolation of AbstractGridArray [#671](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/671)
- Test folder sorted into subfolders [#671](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/671)
- Land model modularised + land netCDF output [#671](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/671)

## v0.14

- CompatHelper: Allow for Makie.jl v0.22 [#663](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/663)
- Prescribed ocean/land heat fluxes modular [#667](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/667)
- NetCDF output for tracers defined, precipitation rate output at initialize! again [#657](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/657)
- Precipitation rate from large-scale condensation/convection for coupling [#669](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/669)
- Fixed a bug, so that sum now works for LowerTraingularArrays [#668](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/668)
- Require only AbstractVectors in the interpolation [#665](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/665)
- Intitialize! simulation restructured [#666](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/666)
- First time steps in main time loop [#659](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/659)
- Generalise OneBandRadiation to NBandRadiation [#633](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/633)
- New output variables definition simplified [#653](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/653)
- Bug: get_vertices for full grids had southpole at 90˚N [#654](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/654)
- Group initialize!, timestep!, finalize! of main time loop [#652](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/652)
- Coordinates defined in degrees lond, latd, not colat, lon [#647](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/647)
- Spectral Gradients are now differentiable, with correctness check in extended CI [#638](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/638)
- set! for time step [#650](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/650)
- Bug: ZonalWind had sqrt of negative for uncommon resolutions [#649](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/649)
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
