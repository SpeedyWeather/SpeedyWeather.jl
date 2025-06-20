@testset "Prescribed surface sensible heat fluxes" begin
    for Model in (PrimitiveDryModel, PrimitiveWetModel)
        tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits

        # prescribe ocean
        spectral_grid = SpectralGrid(trunc=31)
        output = NetCDFOutput(spectral_grid, path=tmp_output_path)
        ocean_heat_flux = PrescribedOceanHeatFlux(spectral_grid)
        land_heat_flux = SurfaceLandHeatFlux(spectral_grid)
        surface_heat_flux = SurfaceHeatFlux(ocean=ocean_heat_flux, land=land_heat_flux)
        model = Model(spectral_grid; surface_heat_flux, output)
        add!(model, SpeedyWeather.SurfaceFluxesOutput()...)

        simulation = initialize!(model)
        set!(simulation.prognostic_variables.ocean.sensible_heat_flux, (λ, ϕ) -> ϕ > 0 ? 10 : 0, model.geometry)
        run!(simulation, period=Day(1), output=true)

        # prescribe land
        ocean_heat_flux = SurfaceOceanHeatFlux(spectral_grid)
        land_heat_flux = PrescribedLandHeatFlux(spectral_grid)
        surface_heat_flux = SurfaceHeatFlux(ocean=ocean_heat_flux, land=land_heat_flux)
        model = Model(spectral_grid; surface_heat_flux, output)
        add!(model, SpeedyWeather.SurfaceFluxesOutput()...)

        simulation = initialize!(model)
        set!(simulation.prognostic_variables.land.sensible_heat_flux, (λ, ϕ) -> ϕ > 0 ? 100 : 0, model.geometry)
        run!(simulation, period=Day(1), output=true)
    end
end

@testset "Prescribed surface evaporative fluxes" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits

    # prescribe ocean
    spectral_grid = SpectralGrid(trunc=31)
    output = NetCDFOutput(spectral_grid, path=tmp_output_path)
    evaporative_flux_ocean = PrescribedOceanEvaporation(spectral_grid)
    evaporative_flux_land = SurfaceLandEvaporation(spectral_grid)
    surface_evaporation = SurfaceEvaporation(ocean=evaporative_flux_ocean, land=evaporative_flux_land)
    model = PrimitiveWetModel(spectral_grid; surface_evaporation, output)
    
    simulation = initialize!(model)
    set!(simulation.prognostic_variables.ocean.evaporative_flux, (λ, ϕ) -> ϕ > 0 ? 5e-5 : 0, model.geometry)
    run!(simulation, period=Day(1))

    # prescribe land
    evaporative_flux_ocean = SurfaceOceanEvaporation(spectral_grid)
    evaporative_flux_land = PrescribedLandEvaporation(spectral_grid)
    surface_evaporation = SurfaceEvaporation(ocean=evaporative_flux_ocean, land=evaporative_flux_land)
    model = PrimitiveWetModel(spectral_grid; surface_evaporation, output)

    simulation = initialize!(model);
    set!(simulation.prognostic_variables.land.evaporative_flux, (λ, ϕ) -> ϕ > 0 ? 5e-5 : 0, model.geometry)
    run!(simulation, period=Day(1))
end