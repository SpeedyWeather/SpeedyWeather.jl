@testset "large_scale_condensation.jl" begin
    @testset "get_large_scale_condensation_tendencies!" begin
        _, diag, model = SpeedyWeather.initialize_speedy()

        # For now, just check that it runs without errors
        SpeedyWeather.get_large_scale_condensation_tendencies!(diag, model)
    end
end
