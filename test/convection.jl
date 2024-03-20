@testset "Convection" begin
    
    spectral_grid = SpectralGrid(trunc=31, nlev=8)

    for Convection in ( NoConvection,
                        SimplifiedBettsMiller,
                        DryBettsMiller)
        for Model in (  PrimitiveDryModel,
                        PrimitiveWetModel)

            convection = Convection(spectral_grid)
            model = Model(;spectral_grid, convection)
            simulation = initialize!(model)
            run!(simulation, period=Day(5))
        end
    end
end