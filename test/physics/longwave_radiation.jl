@testset "Longwave radiation" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)
    @testset for LW in (Nothing, UniformCooling, JeevanjeeRadiation, OneBandGreyLongwave, OneBandLongwave)
        longwave_radiation = LW(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; longwave_radiation)

        initialize!(model.longwave_radiation, model)

        progn = PrognosticVariables(model)
        diagn = DiagnosticVariables(model)

        ij = rand(1:model.spectral_grid.npoints)
        SpeedyWeather.parameterization!(ij, diagn, progn, model.longwave_radiation, model)
    end
end

@testset "Longwave Transmissivity" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    @testset for T in (FriersonLongwaveTransmissivity, TransparentLongwaveTransmissivity)
        transmissivity = T(spectral_grid)
        longwave_radiation = OneBandLongwave(spectral_grid; transmissivity)
        model = PrimitiveWetModel(spectral_grid; longwave_radiation)
        initialize!(model.longwave_radiation, model)

        progn = PrognosticVariables(model)
        diagn = DiagnosticVariables(model)

        ij = rand(1:model.spectral_grid.npoints)
        SpeedyWeather.parameterization!(ij, diagn, progn, model.longwave_radiation, model)
        t = SpeedyWeather.transmissivity!(ij, diagn, progn, model.longwave_radiation.transmissivity, model)

        @test all(0 .<= t[ij, :] .<= 1)
    end
end