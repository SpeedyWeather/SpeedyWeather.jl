@testset "Longwave radiation" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    @testset for LW in (Nothing, UniformCooling, JeevanjeeRadiation, OneBandGreyLongwave, OneBandLongwave)
        longwave = LW(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; radiation = Radiation(spectral_grid; longwave))

        initialize!(model.radiation, model)

        vars = Variables(model)

        ij = rand(1:model.spectral_grid.npoints)
        SpeedyWeather.parameterization!(ij, vars, model.radiation.longwave, model)
    end
end

@testset "Longwave Transmissivity" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

    @testset for T in (FriersonLongwaveTransmissivity, ByrneOGormanLongwaveTransmissivity, TransparentLongwaveTransmissivity)
        transmissivity = T(spectral_grid)
        longwave = OneBandLongwave(spectral_grid; transmissivity)
        model = PrimitiveWetModel(spectral_grid; radiation = Radiation(spectral_grid; longwave))
        @test model.radiation.longwave === longwave
        initialize!(model.radiation, model)

        vars = Variables(model)

        # transmissivity depends on pressure thickness and thereofre surface pressure should be nonzero
        vars.parameterizations.surface_pressure .= 1.0e5
        t = SpeedyWeather.transmissivity!(1, vars, model.radiation.longwave.transmissivity, model)
        for ij in 2:model.spectral_grid.npoints
            SpeedyWeather.transmissivity!(ij, vars, model.radiation.longwave.transmissivity, model)
        end

        @test all(0 .< t .<= 1)
    end
end

@testset "Byrne and O'Gorman humidity-dependent longwave transmissivity" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 8)
    transmissivity = ByrneOGormanLongwaveTransmissivity(spectral_grid)
    longwave = OneBandLongwave(spectral_grid; transmissivity)
    model = PrimitiveWetModel(spectral_grid; radiation = Radiation(spectral_grid; longwave))
    initialize!(model.radiation, model)
    vars = Variables(model)
    vars.parameterizations.surface_pressure .= 1.0e5

    # dry: only the dry optical depth, summing to dry_absorption over the column
    vars.grid.humidity .= 0
    t_dry = copy(SpeedyWeather.transmissivity!(1, vars, transmissivity, model)[1, :])
    @test prod(t_dry) ≈ exp(-transmissivity.dry_absorption) rtol = 1.0e-5

    # moist: lower transmissivity in every layer
    vars.grid.humidity .= 0.005
    t_moist = SpeedyWeather.transmissivity!(1, vars, transmissivity, model)[1, :]
    @test all(t_moist .< t_dry)
    @test prod(t_moist) ≈ exp(-(transmissivity.dry_absorption + transmissivity.water_vapor_absorption * 0.005)) rtol = 1.0e-4
end
