@testset "Longwave radiation" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
    @testset for LW in (Nothing, UniformCooling, JeevanjeeRadiation, OneBandGreyLongwave, OneBandLongwave)
        longwave_radiation = LW(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; longwave_radiation)

        initialize!(model.longwave_radiation, model)

        vars = Variables(model)
        # fill the precomputed (state-independent) transmissivity field; no-op unless OneBandLongwave
        SpeedyWeather.initialize!(vars, model.longwave_radiation, model)

        ij = rand(1:model.spectral_grid.npoints)
        SpeedyWeather.parameterization!(ij, vars, model.longwave_radiation, model)
    end
end

@testset "Longwave Transmissivity" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
    (; nlayers, nlat, npoints) = spectral_grid

    @testset for T in (FriersonLongwaveTransmissivity, TransparentLongwaveTransmissivity, ConstantLongwaveTransmissivity)
        transmissivity = T(spectral_grid)
        longwave_radiation = OneBandLongwave(spectral_grid; transmissivity)
        model = PrimitiveWetModel(spectral_grid; longwave_radiation)
        initialize!(model.longwave_radiation, model)

        vars = Variables(model)
        SpeedyWeather.initialize!(vars, model.longwave_radiation, model)   # precompute t once on host

        # precomputed field is nlayers × nlat and a valid transmissivity in [0, 1]
        t_pre = vars.parameterizations.longwave_transmissivity
        @test size(t_pre) == (nlayers, nlat)
        @test all(0 .<= t_pre .<= 1)

        # the on-the-fly transmissivity! reads the actual surface pressure while the precompute
        # uses the reference pressure. For SigmaCoordinates (the default here) the transmissivity
        # is independent of surface pressure, so with any nonzero pₛ the two paths must be
        # bit-identical at every grid point; we use the reference pressure to make that explicit.
        # This guards that the precomputed and on-the-fly code paths can't drift apart.
        vars.parameterizations.surface_pressure .= model.atmosphere.reference_pressure
        whichring = model.geometry.whichring
        for ij in 1:npoints
            SpeedyWeather.parameterization!(ij, vars, model.longwave_radiation, model)
            t = SpeedyWeather.transmissivity!(ij, vars, model.longwave_radiation.transmissivity, model)
            @test t[ij, :] == t_pre[:, whichring[ij]]
        end
    end
end
