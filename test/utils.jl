@testset "Increasing/decresing vectors" begin
    @test SpeedyWeather.isincreasing(collect(1:10))
    @test SpeedyWeather.isincreasing(sort(rand(10)))
    @test ~SpeedyWeather.isincreasing(rand(Float32,10))
    @test ~SpeedyWeather.isincreasing(randn(10))
end

@testset "clip negatives" begin
    for T in (Float16,Float32,Float64)
        A = randn(T,30,50)
        SpeedyWeather.clip_negatives!(A)
        @test all(A .>= 0)
    end
end