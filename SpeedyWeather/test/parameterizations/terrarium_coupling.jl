using Terrarium
using Test
using Dates

# The SpeedyWeatherTerrariumExt extension activates once both packages are
# loaded. Pull the public names out of the extension module.
const SWTerrarium = Base.get_extension(SpeedyWeather, :SpeedyWeatherTerrariumExt)
@assert SWTerrarium !== nothing "SpeedyWeatherTerrariumExt failed to load"
const TerrariumLand = SWTerrarium.TerrariumLand
const AbstractTerrariumLandModel = SWTerrarium.AbstractTerrariumLandModel

# helper: build the coupled wet model from a given Terrarium land mask (boolean Field)
function build_terrarium_wet_model(ring_grid, land_sea_mask, terrarium_mask; Δt = 300.0, kwargs...)
    spectral_grid = SpectralGrid(ring_grid)
    Nz = 4
    Δz_min = 0.05
    column_grid = Terrarium.ColumnRingGrid(
        Terrarium.CPU(), Float32,
        Terrarium.ExponentialSpacing(; N = Nz, Δz_min),
        ring_grid,
        terrarium_mask,
    )
    soil_initializer = Terrarium.SoilInitializer(eltype(column_grid))
    soil = Terrarium.SoilEnergyWaterCarbon(
        eltype(column_grid);
        hydrology = Terrarium.SoilHydrology(eltype(column_grid)),
    )
    terrarium_model = Terrarium.LandModel(
        column_grid; initializer = soil_initializer, vegetation = nothing, soil,
    )
    land = SpeedyWeather.LandModel(spectral_grid, terrarium_model; Δt, kwargs...)
    model = PrimitiveWetModel(
        spectral_grid;
        land,
        land_sea_mask,
        surface_heat_flux = SurfaceHeatFlux(spectral_grid, land = PrescribedLandHeatFlux()),
        surface_humidity_flux = SurfaceHumidityFlux(spectral_grid, land = PrescribedLandHumidityFlux()),
        time_stepping = Leapfrog(spectral_grid, Δt_at_T32 = Minute(15)),
    )
    return model
end

@testset "Terrarium coupling: type hierarchy" begin
    @test AbstractTerrariumLandModel <: SpeedyWeather.AbstractLand
    @test TerrariumLand <: AbstractTerrariumLandModel
end

@testset "Terrarium coupling: initialize + run + NaN check" begin
    # Small ring grid + matching Terrarium column grid; keep Nz tiny so the
    # test runs in seconds even on CI. Derive the boolean Terrarium mask from the
    # fractional SpeedyWeather land-sea mask (a column wherever there is any land).
    ring_grid = SpeedyWeather.RingGrids.FullGaussianGrid(12)
    spectral_grid = SpectralGrid(ring_grid)
    land_sea_mask = EarthLandSeaMask(spectral_grid)
    SpeedyWeather.load_mask!(land_sea_mask)

    # a maximum sub-step that divides neither the SpeedyWeather step nor the Leapfrog half steps
    model = build_terrarium_wet_model(ring_grid, land_sea_mask, land_sea_mask.land_fraction .> 0; Δt = Second(421))
    @test model.land isa AbstractTerrariumLandModel
    @test model.land isa SpeedyWeather.AbstractLand
    @test SpeedyWeather.get_nlayers(model.land) == 1

    sim = SpeedyWeather.initialize!(model)

    # The Terrarium state and the SpeedyWeather-side mirrors must all be
    # allocated in the Variables tree under the :land namespace.
    @test haskey(sim.variables.prognostic.land, :terrarium)
    @test haskey(sim.variables.prognostic.land, :soil_temperature)
    @test haskey(sim.variables.prognostic.land, :soil_moisture)

    # Run a handful of steps and confirm nothing has gone NaN.
    SpeedyWeather.run!(sim, steps = 3)

    # Terrarium sub-steps each SpeedyWeather step and both clocks stay in sync
    @test sim.variables.prognostic.land.terrarium.clock.iteration > 3
    @test sim.variables.prognostic.land.terrarium.clock.time == sim.variables.prognostic.clock.time

    @test all(isfinite, sim.variables.prognostic.land.soil_temperature)
    @test all(isfinite, sim.variables.prognostic.land.soil_moisture)
    @test all(isfinite, sim.variables.grid.temperature)
    @test all(isfinite, sim.variables.grid.humidity)
    @test all(isfinite, sim.variables.grid.pressure)

    # Spot-check the Terrarium state itself for NaNs in the soil column.
    state = sim.variables.prognostic.land.terrarium
    @test all(isfinite, Terrarium.interior(state.temperature))
    @test all(isfinite, Terrarium.interior(state.saturation_water_ice))
    @test all(isfinite, Terrarium.interior(state.skin_temperature))

    # The Terrarium state can be moved to the host, e.g. for restart files
    host_state = SpeedyWeather.on_architecture(SpeedyWeather.CPU(), state)
    @test host_state isa Terrarium.StateVariables
    @test Terrarium.interior(host_state.temperature) == Terrarium.interior(state.temperature)
    @test SpeedyWeather.on_architecture(SpeedyWeather.CPU(), sim.variables.prognostic) isa NamedTuple

    # Check soil moisture capacity works correctly
    land_vars = sim.variables.parameterizations.land
    @test haskey(land_vars, :soil_moisture_capacity)
    mask = SWTerrarium.land_mask(model.land)
    @test all(0 .< land_vars.soil_moisture_capacity[mask] .<= 1)   # porosity is a volume fraction
    @test all(isfinite, land_vars.albedo)                           # LearnedLandAlbedo now sees real soil moisture
end

@testset "Terrarium coupling: synchronized time steps and Adapt" begin
    ring_grid = SpeedyWeather.RingGrids.FullGaussianGrid(12)
    spectral_grid = SpectralGrid(ring_grid)
    land_sea_mask = EarthLandSeaMask(spectral_grid)
    SpeedyWeather.load_mask!(land_sea_mask)

    # by default (`Δt = nothing`), Terrarium takes one step per SpeedyWeather step of the same length
    model = build_terrarium_wet_model(ring_grid, land_sea_mask, land_sea_mask.land_fraction .> 0; Δt = nothing)
    @test isnothing(model.land.Δt)
    sim = SpeedyWeather.initialize!(model)
    nsteps = 3
    SpeedyWeather.run!(sim, steps = nsteps)
    terrarium_clock = sim.variables.prognostic.land.terrarium.clock
    @test terrarium_clock.iteration == nsteps
    @test terrarium_clock.time == sim.variables.prognostic.clock.time

    # adapting the variables (e.g. into kernel arguments) adapts the mutable Terrarium state to a NamedTuple
    # (Terrarium's own `adapt_structure`) and keeps everything else
    adapted = SpeedyWeather.Adapt.adapt(Array, sim.variables)
    @test adapted.prognostic.land.terrarium isa NamedTuple
    @test adapted.prognostic.land.soil_temperature == sim.variables.prognostic.land.soil_temperature
    @test adapted.grid.temperature == sim.variables.grid.temperature
    @test sim.variables.prognostic.land.terrarium isa Terrarium.StateVariables

    # a maximum sub-step `Δt` splits each SpeedyWeather step into the fewest equal sub-steps not longer than `Δt`
    @test SWTerrarium.terrarium_substeps(nothing, Millisecond(Minute(20))) == 1
    @test SWTerrarium.terrarium_substeps(Minute(10), Millisecond(Minute(20))) == 2
    @test SWTerrarium.terrarium_substeps(Minute(7), Millisecond(Minute(20))) == 3
    @test SWTerrarium.terrarium_substeps(600.0f0, Millisecond(Minute(15))) == 2
    @test_throws ArgumentError build_terrarium_wet_model(ring_grid, land_sea_mask, land_sea_mask.land_fraction .> 0; Δt = 0.0)
end

@testset "Terrarium coupling (dry): initialize + run + NaN check" begin
    ring_grid = SpeedyWeather.RingGrids.FullGaussianGrid(12)
    spectral_grid = SpectralGrid(ring_grid)

    Nz = 4
    Δz_min = 0.05
    column_grid = Terrarium.ColumnRingGrid(
        Terrarium.CPU(), Float32,
        Terrarium.ExponentialSpacing(; N = Nz, Δz_min),
        ring_grid,
    )

    soil_initializer = Terrarium.SoilInitializer(eltype(column_grid))
    # the value boundary condition below (`PrescribedSurfaceTemperature`) requires an explicit timestepper
    timestepper = Terrarium.ForwardEuler(eltype(column_grid))
    soil_model = Terrarium.SoilModel(column_grid; initializer = soil_initializer, timestepper)

    # InputSource for air temperature — consumed by the dry land boundary condition
    air_temperature_field = Terrarium.Field(column_grid, Terrarium.XY())
    Tair_input = Terrarium.InputSource(column_grid, air_temperature_field; name = :air_temperature)
    bcs = Terrarium.PrescribedSurfaceTemperature(:air_temperature)

    land = SpeedyWeather.LandModel(
        spectral_grid, soil_model;
        inputs = Tair_input,
        boundary_conditions = bcs,
        Δt = 421.0,     # maximum sub-step (seconds) that divides neither the step nor the Leapfrog half steps
    )
    @test land isa AbstractTerrariumLandModel
    @test land isa SpeedyWeather.AbstractLand
    @test SpeedyWeather.get_nlayers(land) == 1

    land_sea_mask = RockyPlanetMask(land.spectral_grid)
    time_stepping = Leapfrog(land.spectral_grid, Δt_at_T32 = Minute(15))

    model = PrimitiveDryModel(
        land.spectral_grid;
        land,
        land_sea_mask,
        time_stepping,
    )

    sim = SpeedyWeather.initialize!(model)

    @test haskey(sim.variables.prognostic.land, :terrarium)
    @test haskey(sim.variables.prognostic.land, :soil_temperature)

    SpeedyWeather.run!(sim, steps = 3)

    @test sim.variables.prognostic.land.terrarium.clock.iteration > 3
    @test sim.variables.prognostic.land.terrarium.clock.time == sim.variables.prognostic.clock.time
    @test all(isfinite, sim.variables.prognostic.land.soil_temperature)
    @test all(isfinite, sim.variables.grid.temperature)
    @test all(isfinite, sim.variables.grid.pressure)

    state = sim.variables.prognostic.land.terrarium
    @test all(isfinite, Terrarium.interior(state.temperature))
end

@testset "Terrarium coupling: land mask helper + convenience constructor" begin
    ring_grid = SpeedyWeather.RingGrids.FullGaussianGrid(12)
    spectral_grid = SpectralGrid(ring_grid)
    land_sea_mask = EarthLandSeaMask(spectral_grid)
    SpeedyWeather.load_mask!(land_sea_mask)

    # land_mask(land_sea_mask) returns a boolean Field equal to land_fraction .> 0
    bool_mask = SWTerrarium.land_mask(land_sea_mask)
    @test eltype(bool_mask) == Bool
    @test all(bool_mask.data .== (land_sea_mask.land_fraction.data .> 0))

    # a stricter threshold selects fewer points
    strict = SWTerrarium.land_mask(land_sea_mask; threshold = 0.5)
    @test all(strict.data .== (land_sea_mask.land_fraction.data .> 0.5))
    @test count(strict.data) < count(bool_mask.data)

    # the convenience ColumnRingGrid constructor derives the same mask from land_sea_mask
    Nz = 4
    Δz_min = 0.05
    spacing = Terrarium.ExponentialSpacing(; N = Nz, Δz_min)
    column_grid = Terrarium.ColumnRingGrid(Terrarium.CPU(), Float32, spacing, ring_grid, land_sea_mask)
    @test all(column_grid.mask.data .== (land_sea_mask.land_fraction.data .> 0))
    @test count(column_grid.mask.data) == count(bool_mask.data)
end

@testset "Terrarium coupling: ocean fallback fill + mask-consistency warning" begin
    ring_grid = SpeedyWeather.RingGrids.FullGaussianGrid(12)
    spectral_grid = SpectralGrid(ring_grid)
    land_sea_mask = EarthLandSeaMask(spectral_grid)
    SpeedyWeather.load_mask!(land_sea_mask)

    # Consistent mask (superset of SpeedyWeather land): no warning expected.
    consistent_mask = SWTerrarium.land_mask(land_sea_mask)
    model = @test_nowarn build_terrarium_wet_model(
        ring_grid, land_sea_mask, consistent_mask;
        ocean_temperature = 290, ocean_moisture = 0,
    )
    @test model.land.ocean_temperature == 290

    sim = SpeedyWeather.initialize!(model)
    lf = land_sea_mask.land_fraction.data
    ocean_pts = lf .== 0
    st = sim.variables.prognostic.land.soil_temperature.data
    sm = sim.variables.prognostic.land.soil_moisture.data

    # ocean-only points hold the fallback values, not the 0 K / 0 allocation default
    @test all(st[ocean_pts] .== 290)
    @test all(sm[ocean_pts] .== 0)
    # land points were seeded from Terrarium (not the fallback), so they differ
    @test all(st[.!ocean_pts] .!= 290)
    @test all(isfinite, st)

    SpeedyWeather.run!(sim, steps = 3)
    # ocean points are never touched by the land step, still the fallback value
    @test all(st[ocean_pts] .== 290)
    @test all(isfinite, st)

    # Stricter Terrarium mask leaves land cells uncovered -> warning must fire.
    land_sea_mask2 = EarthLandSeaMask(spectral_grid)
    SpeedyWeather.load_mask!(land_sea_mask2)
    strict_mask = SWTerrarium.land_mask(land_sea_mask2; threshold = 0.5)
    model2 = build_terrarium_wet_model(ring_grid, land_sea_mask2, strict_mask)
    @test_logs (:warn,) match_mode = :any SpeedyWeather.initialize!(model2)

    # orphaned points (land_fraction > 0 without a column) get the fallback, not 0 K
    sim2 = SpeedyWeather.initialize!(model2)
    orphaned = (land_sea_mask2.land_fraction.data .> 0) .& .!strict_mask.data
    @test any(orphaned)     # sanity: the stricter mask really does orphan some land
    @test all(sim2.variables.prognostic.land.soil_temperature.data[orphaned] .== model2.land.ocean_temperature)
end
