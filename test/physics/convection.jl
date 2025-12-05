@testset "Convection" begin
    
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    for Convection in ( Nothing,
                        BettsMillerConvection,
                        BettsMillerDryConvection,
                        ConvectiveHeating)
        for Model in (  PrimitiveDryModel,
                        PrimitiveWetModel)

            # that combination is not defined
            if ~(Convection == BettsMillerConvection && Model == PrimitiveDryModel)
                convection = Convection(spectral_grid)
                model = Model(spectral_grid; convection)
                model.feedback.verbose = false
                simulation = initialize!(model)

                run!(simulation, steps=36)
            end
        end
    end
end