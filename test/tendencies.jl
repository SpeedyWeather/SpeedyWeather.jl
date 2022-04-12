@testset "Consistent run through of tendencies" begin
    # No tests, just make sure get_tendencies runs without bugs
    # No values calculated yet, need consistent definitions of grads, spectral transforms etc.
    for NF in (Float32,Float64)
        P,D,M = SpeedyWeather.initialize_model(NF)
        l2 = 2
        SpeedyWeather.get_tendencies!(P,D,l2,M)
    end    
end