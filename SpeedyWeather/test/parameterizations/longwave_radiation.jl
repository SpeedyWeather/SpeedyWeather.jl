@testset "Longwave radiation" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
    @testset for LW in (Nothing, UniformCooling, JeevanjeeRadiation, OneBandGreyLongwave, OneBandLongwave)
        longwave_radiation = LW(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; longwave_radiation)

        initialize!(model.longwave_radiation, model)

        vars = Variables(model)

        ij = rand(1:model.spectral_grid.npoints)
        SpeedyWeather.parameterization!(ij, vars, model.longwave_radiation, model)
    end
end

@testset "Longwave Transmissivity" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)

    @testset for T in (FriersonLongwaveTransmissivity, TransparentLongwaveTransmissivity)
        transmissivity = T(spectral_grid)
        longwave_radiation = OneBandLongwave(spectral_grid; transmissivity)
        model = PrimitiveWetModel(spectral_grid; longwave_radiation)
        initialize!(model.longwave_radiation, model)

        vars = Variables(model)

        # transmissivity depends on pressure thickness and thereofre surface pressure should be nonzero
        vars.parameterizations.surface_pressure .= 1e5
        t = SpeedyWeather.transmissivity!(1, vars, model.longwave_radiation.transmissivity, model)
        for ij in 2:model.spectral_grid.npoints
            SpeedyWeather.transmissivity!(ij, vars, model.longwave_radiation.transmissivity, model)
        end

        @test all(0 .< t .<= 1)
    end
end
