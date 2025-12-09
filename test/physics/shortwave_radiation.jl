@testset "Shortwave radiation" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)
    @testset for SW in (Nothing, TransparentShortwave,)# OneBandShortwave, OneBandGreyShortwave)
        sw = SW(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)

        initialize!(model.shortwave_radiation, model)

        progn = PrognosticVariables(model)
        diagn = DiagnosticVariables(model)

        (; time) = progn.clock
        SpeedyWeather.cos_zenith!(diagn, time, model)

        ij = rand(1:model.spectral_grid.npoints)
        diagn.physics.ocean.albedo[ij] = 0.5
        SpeedyWeather.parameterization!(ij, diagn, progn, model.shortwave_radiation, model)

        trd = model.planet.solar_constant * diagn.physics.cos_zenith[ij]

        if !(sw isa Nothing)
            osr = diagn.physics.outgoing_shortwave[ij]
            ssrd = diagn.physics.surface_shortwave_down[ij]
            @test 0 < osr < ssrd <= trd
            @test isfinite(osr)
            @test isfinite(ssrd)
        end
    end
end

# To be adapted when translated to new ij-based structure
# @testset "Shortwave radiation transmissivity" begin
#     spectral_grid = SpectralGrid(trunc=31, nlayers=8)
#     @testset for T in (TransparentShortwaveTransmissivity, BackgroundShortwaveTransmissivity)
#         transmissivity = T(spectral_grid)
#         sw = OneBandShortwave(spectral_grid; transmissivity=transmissivity)
#         model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)

#         initialize!(model.shortwave_radiation, model)

#         progn = PrognosticVariables(model)
#         diagn = DiagnosticVariables(model)

#         ij = rand(1:model.spectral_grid.npoints)
#         SpeedyWeather.parameterization!(ij, diagn, progn, model.shortwave_radiation, model)
#     end
# end

# @testset "Shortwave radiation clouds" begin
#     spectral_grid = SpectralGrid(trunc=31, nlayers=8)
#     @testset for C in (DiagnosticClouds, NoClouds)
#         clouds = C(spectral_grid)
#         sw = OneBandShortwave(spectral_grid; clouds=clouds)
#         model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)

#         initialize!(model.shortwave_radiation, model)

#         progn = PrognosticVariables(model)
#         diagn = DiagnosticVariables(model)

#         ij = rand(1:model.spectral_grid.npoints)
#         SpeedyWeather.parameterization!(ij, diagn, progn, model.shortwave_radiation, model)
#     end
# end