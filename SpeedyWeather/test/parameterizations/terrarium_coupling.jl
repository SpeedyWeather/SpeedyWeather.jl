using Terrarium
using Test
using Dates

# The SpeedyWeatherTerrariumExt extension activates once both packages are
# loaded. Pull the public names out of the extension module.
const SWTerrarium = Base.get_extension(SpeedyWeather, :SpeedyWeatherTerrariumExt)
@assert SWTerrarium !== nothing "SpeedyWeatherTerrariumExt failed to load"
const TerrariumLand = SWTerrarium.TerrariumLand
const AbstractTerrariumLandModel = SWTerrarium.AbstractTerrariumLandModel

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
