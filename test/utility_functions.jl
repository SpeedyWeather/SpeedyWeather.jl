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

@testset "roundup nlon for FFT" begin
    for i in 1:10
        @test 2^i == SpeedyWeather.roundup_fft(2^i)
        @test 2^i*3 == SpeedyWeather.roundup_fft(2^i*3)
    end
    for n in 1:10
        i = rand(1:1000)
        @test i <= SpeedyWeather.roundup_fft(i)
    end
end