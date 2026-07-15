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
function build_terrarium_wet_model(ring_grid, land_sea_mask, terrarium_mask; kwargs...)
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
    land = SpeedyWeather.LandModel(spectral_grid, terrarium_model; Δt = 300.0, kwargs...)
    model = PrimitiveWetModel(
        spectral_grid;
        land,
        land_sea_mask,
        surface_heat_flux = SurfaceHeatFlux(spectral_grid, land = PrescribedLandHeatFlux()),
        surface_humidity_flux = SurfaceHumidityFlux(spectral_grid, land = PrescribedLandHumidityFlux()),
        time_stepping = Leapfrog(spectral_grid, Δt_at_T31 = Minute(15)),
    )
    return model
end

@testset "Terrarium coupling: type hierarchy" begin
    @test AbstractTerrariumLandModel <: SpeedyWeather.AbstractLand
    @test TerrariumLand <: AbstractTerrariumLandModel
end

@testset "Terrarium coupling: initialize + run + NaN check" begin
    # Small ring grid + matching Terrarium column grid; keep Nz tiny so the
    # test runs in seconds even on CI.
    ring_grid = SpeedyWeather.RingGrids.FullGaussianGrid(12)
    spectral_grid = SpectralGrid(ring_grid)
    land_sea_mask = EarthLandSeaMask(spectral_grid)
    SpeedyWeather.load_mask!(land_sea_mask)

    Nz = 4
    Δz_min = 0.05
    # Terrarium needs a boolean land mask: a column is allocated wherever there is
    # any land. Derive it from the fractional SpeedyWeather land-sea mask.
    land_mask = land_sea_mask.land_fraction .> 0
    column_grid = Terrarium.ColumnRingGrid(
        Terrarium.CPU(), Float32,
        Terrarium.ExponentialSpacing(; N = Nz, Δz_min),
        ring_grid,
        land_mask,
    )

    soil_initializer = Terrarium.SoilInitializer(eltype(column_grid))
    soil = Terrarium.SoilEnergyWaterCarbon(
        eltype(column_grid);
        hydrology = Terrarium.SoilHydrology(eltype(column_grid)),
    )
    terrarium_model = Terrarium.LandModel(
        column_grid;
        initializer = soil_initializer,
        vegetation = nothing,
        soil,
    )

    land = SpeedyWeather.LandModel(spectral_grid, terrarium_model; Δt = 300.0)
    @test land isa AbstractTerrariumLandModel
    @test land isa SpeedyWeather.AbstractLand
    @test SpeedyWeather.get_nlayers(land) == 1

    surface_heat_flux = SurfaceHeatFlux(land.spectral_grid, land = PrescribedLandHeatFlux())
    surface_humidity_flux = SurfaceHumidityFlux(land.spectral_grid, land = PrescribedLandHumidityFlux())
    time_stepping = Leapfrog(land.spectral_grid, Δt_at_T31 = Minute(15))

    model = PrimitiveWetModel(
        land.spectral_grid;
        land,
        surface_heat_flux,
        surface_humidity_flux,
        land_sea_mask,
        time_stepping,
    )

    sim = SpeedyWeather.initialize!(model)

    # The Terrarium state and the SpeedyWeather-side mirrors must all be
    # allocated in the Variables tree under the :land namespace.
    @test haskey(sim.variables.prognostic.land, :terrarium)
    @test haskey(sim.variables.prognostic.land, :soil_temperature)
    @test haskey(sim.variables.prognostic.land, :soil_moisture)

    # Run a handful of steps and confirm nothing has gone NaN.
    SpeedyWeather.run!(sim, steps = 3)

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
    soil_model = Terrarium.SoilModel(column_grid; initializer = soil_initializer)

    # InputSource for air temperature — consumed by the dry land boundary condition
    air_temperature_field = Terrarium.Field(column_grid, Terrarium.XY())
    Tair_input = Terrarium.InputSource(column_grid, air_temperature_field; name = :air_temperature)
    bcs = Terrarium.PrescribedSurfaceTemperature(:air_temperature)

    land = SpeedyWeather.LandModel(
        spectral_grid, soil_model;
        boundary_conditions = bcs,
        input_variables = Terrarium.variables(Tair_input),
        Δt = 300.0,
    )
    @test land isa AbstractTerrariumLandModel
    @test land isa SpeedyWeather.AbstractLand
    @test SpeedyWeather.get_nlayers(land) == 1

    land_sea_mask = RockyPlanetMask(land.spectral_grid)
    time_stepping = Leapfrog(land.spectral_grid, Δt_at_T31 = Minute(15))

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
    @test model.land.mask == true
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
