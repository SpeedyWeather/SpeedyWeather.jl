@testset "Parametrization: large scale condensation" begin @testset for NF in (Float32,
                                                                               Float64)
    _, diagn, model = SpeedyWeather.initialize_speedy(NF, model = PrimitiveEquation)

    column = ColumnVariables{NF}(nlev = diagn.nlev)

    for ij in SpeedyWeather.eachgridpoint(diagn)
        SpeedyWeather.reset_column!(column)
        SpeedyWeather.get_column!(column, diagn, ij, model.geometry)
        SpeedyWeather.large_scale_condensation!(column, model)
    end
end end
