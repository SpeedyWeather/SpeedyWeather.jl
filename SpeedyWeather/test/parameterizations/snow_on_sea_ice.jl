@testset "Snow on sea ice: flag off is inert" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 5)

    # default is off, so no snow depth variable in the ocean namespace
    model = PrimitiveWetModel(spectral_grid)
    vars = Variables(model)
    @test model.land.snow.sea_ice_snow == false
    @test !haskey(vars.prognostic.ocean, :snow_depth)

    # explicitly on: variable is allocated and initialized to zero
    land = LandModel(spectral_grid, snow = SnowModel(spectral_grid, sea_ice_snow = true))
    model = PrimitiveWetModel(spectral_grid; land)
    simulation = initialize!(model)
    @test haskey(simulation.variables.prognostic.ocean, :snow_depth)
    @test all(iszero, simulation.variables.prognostic.ocean.snow_depth)
end

@testset "Snow on sea ice: accumulation, melt and cap" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 5)
    snow = SnowModel(spectral_grid, sea_ice_snow = true, melting_threshold = 275, snow_depth_cap = 0.5)
    land = LandModel(spectral_grid; snow)
    model = PrimitiveWetModel(spectral_grid; land)
    simulation = initialize!(model)
    vars = simulation.variables

    (; Δt) = model.time_stepping
    S = vars.prognostic.ocean.snow_depth
    ℵ = vars.prognostic.ocean.sea_ice_concentration
    snow_fall_rate = vars.parameterizations.snow_rate

    # fully ice covered everywhere, so the land-sea mask is the only mask left
    fill!(ℵ, 1)
    ocean = model.land_sea_mask.land_fraction .< 1

    # ACCUMULATION: constant snow fall, no melt (surface air temperature below threshold)
    fill!(snow_fall_rate, 1.0e-5)
    fill!(vars.parameterizations.surface_air_temperature, 250)
    fill!(S, 0)

    SpeedyWeather.sea_ice_snow_timestep!(vars, snow, snow_fall_rate, model)
    @test all(S[ocean] .≈ 1.0e-5 * Δt)       # depth per ice area, not weighted by concentration

    SpeedyWeather.sea_ice_snow_timestep!(vars, snow, snow_fall_rate, model)
    @test all(S[ocean] .≈ 2.0e-5 * Δt)       # accumulates further

    # MELT: no snow fall, surface air temperature above the melting threshold
    fill!(snow_fall_rate, 0)
    fill!(vars.parameterizations.surface_air_temperature, 285)
    S_before = copy(S)
    SpeedyWeather.sea_ice_snow_timestep!(vars, snow, snow_fall_rate, model)
    @test all(S[ocean] .< S_before[ocean])   # monotonically decreasing
    @test all(S .>= 0)                       # floored at zero

    # melt to exhaustion never goes negative
    for _ in 1:50
        SpeedyWeather.sea_ice_snow_timestep!(vars, snow, snow_fall_rate, model)
    end
    @test all(S .>= 0)

    # CAP: heavy snow fall, no melt, capped at snow_depth_cap
    fill!(snow_fall_rate, 1)
    fill!(vars.parameterizations.surface_air_temperature, 250)
    for _ in 1:5
        SpeedyWeather.sea_ice_snow_timestep!(vars, snow, snow_fall_rate, model)
    end
    @test maximum(S) <= snow.snow_depth_cap
end

@testset "Snow on sea ice: disappears with the sea ice" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 5)
    snow = SnowModel(spectral_grid, sea_ice_snow = true)
    land = LandModel(spectral_grid; snow)
    model = PrimitiveWetModel(spectral_grid; land)
    simulation = initialize!(model)
    vars = simulation.variables

    S = vars.prognostic.ocean.snow_depth
    ℵ = vars.prognostic.ocean.sea_ice_concentration
    snow_fall_rate = vars.parameterizations.snow_rate

    fill!(S, 0.2)                            # start from a nonzero snow cover on ice
    fill!(snow_fall_rate, 0)
    fill!(vars.parameterizations.surface_air_temperature, 250)   # no melt

    fill!(ℵ, 1)                              # ice present: snow survives
    SpeedyWeather.sea_ice_snow_timestep!(vars, snow, snow_fall_rate, model)
    @test all(S[model.land_sea_mask.land_fraction .< 1] .≈ 0.2)

    fill!(ℵ, 0)                              # ice gone: snow goes with it
    SpeedyWeather.sea_ice_snow_timestep!(vars, snow, snow_fall_rate, model)
    @test all(iszero, S[model.land_sea_mask.land_fraction .< 1])
end

@testset "Snow on sea ice: albedo response" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 5)
    albedo = OceanLandAlbedo(ocean = OceanSeaIceAlbedo(spectral_grid), land = LandSnowAlbedo(spectral_grid))
    snow = SnowModel(spectral_grid, sea_ice_snow = true)
    land = LandModel(spectral_grid; snow)
    model = PrimitiveWetModel(spectral_grid; land, albedo)
    simulation = initialize!(model)
    vars = simulation.variables

    scheme = model.albedo.ocean
    a = vars.parameterizations.ocean.albedo
    S = vars.prognostic.ocean.snow_depth
    ℵ = vars.prognostic.ocean.sea_ice_concentration

    ij = 1
    fill!(ℵ, 1)                              # full ice cover
    fill!(S, 0)
    SpeedyWeather.albedo!(ij, a, vars, scheme, model)
    albedo_no_snow = a[ij]

    fill!(S, 0.2)                            # snow raises the albedo over ice
    SpeedyWeather.albedo!(ij, a, vars, scheme, model)
    @test a[ij] > albedo_no_snow
    @test 0 <= a[ij] <= 1

    fill!(ℵ, 0)                              # no ice: snow has no effect
    SpeedyWeather.albedo!(ij, a, vars, scheme, model)
    albedo_snow_no_ice = a[ij]
    fill!(S, 0)
    SpeedyWeather.albedo!(ij, a, vars, scheme, model)
    @test a[ij] ≈ albedo_snow_no_ice
    @test a[ij] ≈ scheme.albedo_ocean

    # clamped to 1 even when albedo_ice + albedo_snow > 1
    scheme2 = OceanSeaIceAlbedo(spectral_grid, albedo_ice = 0.9, albedo_snow = 0.9, snow_depth_scale = 0.01)
    fill!(ℵ, 1)
    fill!(S, 10)
    SpeedyWeather.albedo!(ij, a, vars, scheme2, model)
    @test a[ij] <= 1
end

@testset "Snow on sea ice: flux insulation" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 5)
    snow = SnowModel(spectral_grid, sea_ice_snow = true)
    land = LandModel(spectral_grid; snow)
    model = PrimitiveWetModel(spectral_grid; land)
    simulation = initialize!(model)
    vars = simulation.variables

    # set the surface state the fluxes read from, not defined before the first time step
    fill!(vars.parameterizations.surface_air_temperature, 280)
    fill!(vars.parameterizations.surface_air_density, 1.2)
    fill!(vars.parameterizations.surface_wind_speed, 10)
    fill!(vars.parameterizations.surface_pressure, 1.0e5)
    fill!(vars.parameterizations.boundary_layer_drag, 1.0e-3)
    fill!(vars.grid.humidity, 0)

    S = vars.prognostic.ocean.snow_depth
    ℵ = vars.prognostic.ocean.sea_ice_concentration

    # pick an ocean point with a defined sea surface temperature, over land the SST is NaN
    SST = SpeedyWeather.get_prognostic_step(vars.prognostic.ocean.sea_surface_temperature, model.time_stepping, model.ocean)
    ij = findfirst(ij -> model.land_sea_mask.land_fraction[ij] < 1 && isfinite(SST[ij]), eachindex(S))
    @test ij !== nothing

    @testset for (flux, field) in (
            (model.surface_heat_flux.ocean, :sensible_heat_flux),
            (model.surface_humidity_flux.ocean, :surface_humidity_flux),
        )
        f! = field == :sensible_heat_flux ? SpeedyWeather.surface_heat_flux! : SpeedyWeather.surface_humidity_flux!
        out = getproperty(vars.parameterizations.ocean, field)

        fill!(ℵ, 1)                          # full ice cover
        fill!(S, 0)
        f!(ij, vars, flux, model)
        flux_no_snow = abs(out[ij])

        fill!(S, 0.5)                        # snow insulates: smaller flux magnitude
        f!(ij, vars, flux, model)
        @test abs(out[ij]) < flux_no_snow

        fill!(ℵ, 0)                          # no ice: snow has no effect
        fill!(S, 0)
        f!(ij, vars, flux, model)
        flux_open_ocean = abs(out[ij])
        fill!(S, 0.5)
        f!(ij, vars, flux, model)
        @test abs(out[ij]) ≈ flux_open_ocean
    end
end
