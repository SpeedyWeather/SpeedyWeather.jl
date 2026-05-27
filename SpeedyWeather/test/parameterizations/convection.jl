@testset "Convection" begin

    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)

    @testset for Convection in (
            Nothing,
            BettsMillerConvection,
            BettsMillerDryConvection,
            ConvectiveHeating,
        )

        @testset for Model in (
                PrimitiveDryModel,
                PrimitiveWetModel,
            )

            # that combination is not defined
            if ~(Convection == BettsMillerConvection && Model == PrimitiveDryModel)
                convection = Convection(spectral_grid)
                model = Model(spectral_grid; convection)
                simulation = initialize!(model; feedback = Feedback(verbose = false))

                run!(simulation, steps = 36)
            end
        end
    end
end
