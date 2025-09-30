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

@testset "Prescribed surface humidity fluxes" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits

    # prescribe ocean
    spectral_grid = SpectralGrid(trunc=31)
    output = NetCDFOutput(spectral_grid, path=tmp_output_path)
    humidity_flux_ocean = PrescribedOceanHumidityFlux(spectral_grid)
    humidity_flux_land = SurfaceLandHumidityFlux(spectral_grid)
    surface_humidity_flux = SurfaceHumidityFlux(ocean=humidity_flux_ocean, land=humidity_flux_land)
    model = PrimitiveWetModel(spectral_grid; surface_humidity_flux, output)
    
    simulation = initialize!(model)
    set!(simulation.prognostic_variables.ocean.surface_humidity_flux, (λ, ϕ) -> ϕ > 0 ? 5e-5 : 0, model.geometry)
    run!(simulation, period=Day(1))

    # prescribe land
    humidity_flux_ocean = SurfaceOceanHumidityFlux(spectral_grid)
    humidity_flux_land = PrescribedLandHumidityFlux(spectral_grid)
    surface_humidity_flux = SurfaceHumidityFlux(ocean=humidity_flux_ocean, land=humidity_flux_land)
    model = PrimitiveWetModel(spectral_grid; surface_humidity_flux, output)

    simulation = initialize!(model);
    set!(simulation.prognostic_variables.land.surface_humidity_flux, (λ, ϕ) -> ϕ > 0 ? 5e-5 : 0, model.geometry)
    run!(simulation, period=Day(1))
end