@testset "Longwave radiation" begin

    spectral_grid = SpectralGrid(trunc=31, nlayers=5)

    for Radiation in (NoLongwave,
                        UniformCooling,
                        JeevanjeeRadiation)

        longwave_radiation = Radiation(spectral_grid)
        
        for Model in (PrimitiveDryModel,
                        PrimitiveWetModel)
            model = Model(spectral_grid; longwave_radiation)
            model.feedback.verbose = false
            simulation = initialize!(model)
            run!(simulation, period=Day(3))
            @test model.feedback.nars_detected == false
        end
    end     
end