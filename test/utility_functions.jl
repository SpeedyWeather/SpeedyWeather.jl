using Dates: CompoundPeriod, Day, Hour, Minute, Second, Millisecond

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
        @test 2^i*5 == SpeedyWeather.roundup_fft(2^i*5)
    end
    for n in 1:10
        i = rand(2:1000)
        @test i <= SpeedyWeather.roundup_fft(i)
    end
end

@testset "readable secs feedback" begin
    @test SpeedyWeather.readable_secs(123456) == CompoundPeriod(Day(1), Hour(10))
    @test SpeedyWeather.readable_secs(12345) == CompoundPeriod(Hour(3), Minute(26))
    @test SpeedyWeather.readable_secs(1234) == CompoundPeriod(Minute(20), Second(34))
    @test SpeedyWeather.readable_secs(123) == CompoundPeriod(Minute(2), Second(3))
    @test SpeedyWeather.readable_secs(12.3) == CompoundPeriod(Second(12), Millisecond(300))
    @test SpeedyWeather.readable_secs(1.23) == CompoundPeriod(Second(1), Millisecond(230))
    @test SpeedyWeather.readable_secs(0.123) == CompoundPeriod(Second(0), Millisecond(120))
end
