using Dates: CompoundPeriod, Day, Hour, Minute, Second, Millisecond, coarserperiod

@testset "Increasing/decresing vectors" begin
    @test SpeedyWeather.isincreasing(collect(1:10))
    @test SpeedyWeather.isincreasing(sort(rand(10)))
    @test ~SpeedyWeather.isincreasing(rand(Float32, 10))
    @test ~SpeedyWeather.isincreasing(randn(10))
end

@testset "clip negatives" begin
    for T in (Float16, Float32, Float64)
        A = randn(T, 30, 50)
        SpeedyWeather.clip_negatives!(A)
        @test all(A .>= 0)
    end
end

@testset "flip sign" begin
    for T in (Float16, Float32, Float64)
        A = randn(T, 30, 50)
        A2 = copy(A)
        SpeedyWeather.flipsign!(A)
        SpeedyWeather.flipsign!(A)
        @test all(A .== A2)
    end
end

@testset "roundup nlon for FFT" begin
    for i in 1:10
        @test 2^i == SpeedyTransforms.roundup_fft(2^i)
        @test 2^i*3 == SpeedyTransforms.roundup_fft(2^i*3)
        @test 2^i*5 == SpeedyTransforms.roundup_fft(2^i*5)
    end
    for n in 1:10
        i = rand(2:1000)
        @test i <= SpeedyTransforms.roundup_fft(i)
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

@testset "nans" begin
    for s in ((3,), (3, 4), (3, 4, 5))
        for T in (Float16, Float32, Float64)
            A = SpeedyWeather.nans(T, s...)
            for a in A
                @test isnan(a)
            end
            @test size(A) == s
        end

        A = SpeedyWeather.nans(s...)
        for a in A
            @test isnan(a)
        end
    end
end

@testset "time convsersions" begin
    # Priated conversions from integer/float types
    @test convert(Second, 1) == Second(1)
    @test convert(Second, 1.4) == Second(1)
    @test convert(Second, 1.5) == Second(2)

    # Pirated conversions for month and year
    @test Second(Year(1)) == Second(365 * 24 * 60 * 60)
    @test Day(Year(1)) == Day(365)
    @test Hour(Year(1)) == Hour(Second(Year(1)))
    @test Minute(Year(1)) == Minute(Second(Year(1)))
    @test Second(Month(1)) == Second(30 * 24 * 60 * 60)
    @test Millisecond(Year(1)) == Millisecond(Second(Year(1)))
    @test Day(Month(1)) == Day(30)
    @test Minute(Month(1)) == Minute(Second(Month(1)))
    @test Millisecond(Month(1)) == Millisecond(Second(Month(1)))
    @test_throws MethodError Week(Month(1))

    # Century
    @test convert(Year, Century(1)) == Year(100)
    @test Century(1) == Year(100)
    @test Second(Century(1)) == Second(100 * 365 * 24 * 60 * 60)
    @test Millisecond(Century(1)) == Millisecond(Second(Century(1)))
    @test coarserperiod(Year) == (Century, 100)
    @test_throws MethodError Week(Century(1))

    # Millenium
    @test convert(Year, Millenium(1)) == Year(1000)
    @test convert(Century, Millenium(1)) == Century(10)
    @test Millenium(1) == Year(1000)
    @test Millenium(1) == Century(10)
    @test Second(Millenium(1)) == Second(1000 * 365 * 24 * 60 * 60)
    @test Millisecond(Millenium(1)) == Millisecond(Second(Millenium(1)))
    @test coarserperiod(Century) == (Millenium, 10)
    @test_throws MethodError Week(Millenium(1))

    # Type promotion rules
    @test promote(Year(1), Century(1)) == (Year(1), Year(100))
    @test promote(Year(1), Millenium(1)) == (Year(1), Year(1000))
    @test promote(Century(1), Millenium(1)) == (Century(1), Century(10))
end
