@testset "Surface skin temperature" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, period = Hour(1))

    vars = simulation.variables
    (; time_stepping) = model
    drag = model.boundary_layer.drag
    land_fraction = Array(model.land_sea_mask.land_fraction.data)
    SST = SpeedyWeather.get_prognostic_step(vars.prognostic.ocean.sea_surface_temperature, time_stepping, drag)
    T_land = vars.prognostic.land.soil_temperature

    # pure ocean point: SST, pure land point: uppermost soil layer
    ij_ocean = findfirst(ij -> land_fraction[ij] == 0 && isfinite(SST[ij]), eachindex(land_fraction))
    ij_land = findfirst(ij -> land_fraction[ij] == 1 && isfinite(T_land[ij, 1]), eachindex(land_fraction))
    @test SpeedyWeather.surface_skin_temperature(ij_ocean, vars, model.land_sea_mask, time_stepping, drag, 0) == SST[ij_ocean]
    @test SpeedyWeather.surface_skin_temperature(ij_land, vars, model.land_sea_mask, time_stepping, drag, 0) == T_land[ij_land, 1]

    # coastal point: weighted by land fraction, bounded by SST and land temperature
    ij_coast = findfirst(ij -> 0 < land_fraction[ij] < 1 && isfinite(SST[ij]) && isfinite(T_land[ij, 1]), eachindex(land_fraction))
    Tₛ = SpeedyWeather.surface_skin_temperature(ij_coast, vars, model.land_sea_mask, time_stepping, drag, 0)
    lf = land_fraction[ij_coast]
    @test Tₛ ≈ (1 - lf) * SST[ij_coast] + lf * T_land[ij_coast, 1]
end

@testset "Bulk Richardson number uses surface temperature" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, period = Hour(1))

    vars = simulation.variables
    (; time_stepping) = model
    drag = model.boundary_layer.drag
    diffusion = model.vertical_diffusion
    nlayers = spectral_grid.nlayers
    land_fraction = Array(model.land_sea_mask.land_fraction.data)
    SST = vars.prognostic.ocean.sea_surface_temperature
    ij = findfirst(ij -> land_fraction[ij] == 0 && isfinite(SST[ij]), eachindex(land_fraction))

    T_air = SpeedyWeather.get_prognostic_step(vars.grid.temperature, time_stepping, drag)[ij, nlayers]
    z = vars.dynamics.geopotential[ij, nlayers] / model.planet.gravity - model.orography.orography[ij]
    ΔΦ₀ = model.planet.gravity * z

    # calm conditions: no resolved wind at the lowermost layer
    for var in (vars.grid.u, vars.grid.v)
        var.data[ij, nlayers, :] .= 0
    end
    vars.parameterizations.surface_wind_speed[ij] = 1     # gusts only

    Ri(ΔSST) = begin
        SST.data[ij, :] .= T_air + ΔSST
        SpeedyWeather.bulk_richardson_surface(ij, ΔΦ₀, vars, model.atmosphere, model.land_sea_mask, drag, time_stepping)
    end

    drag_coefficient(ΔSST) = begin
        SST.data[ij, :] .= T_air + ΔSST
        SpeedyWeather.parameterization!(ij, vars, drag, model)
        vars.parameterizations.boundary_layer_drag[ij]
    end

    boundary_layer_top(ΔSST) = begin
        SST.data[ij, :] .= T_air + ΔSST
        SpeedyWeather.parameterization!(ij, vars, diffusion, model)
        vars.parameterizations.boundary_layer_height[ij]
    end

    # surface warmer than the air: unstable, colder: stable
    @test Ri(10) < 0
    @test Ri(-10) > 0

    # unstable: maximum drag, not the stable minimum drag even in calm conditions
    z₀ = vars.parameterizations.surface_roughness[ij]
    drag_max = (drag.von_Karman / log(z / z₀))^2
    @test drag_coefficient(10) ≈ drag_max
    @test drag_coefficient(-10) < drag_coefficient(10)
    @test drag_coefficient(-10) == drag.drag_min

    # warm surface in calm conditions: there is a boundary layer to diffuse in (top index ≤ nlayers)
    @test boundary_layer_top(10) <= nlayers
    @test boundary_layer_top(-10) == nlayers + 1
end
