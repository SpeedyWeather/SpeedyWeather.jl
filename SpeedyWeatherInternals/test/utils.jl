using SpeedyWeatherInternals.Utils

@testset "Increasing/decresing vectors" begin
    @test isincreasing(collect(1:10))
    @test isincreasing(sort(rand(10)))
    @test ~isincreasing(rand(Float32, 10))
    @test ~isincreasing(randn(10))
end

@testset "clip negatives" begin
    for T in (Float16, Float32, Float64)
        A = randn(T, 30, 50)
        clip_negatives!(A)
        @test all(A .>= 0)
    end
end

@testset "flip sign" begin
    for T in (Float16, Float32, Float64)
        A = randn(T, 30, 50)
        A2 = copy(A)
        flipsign!(A)
        flipsign!(A)
        @test all(A .== A2)
    end
end

@testset "nans" begin
    for s in ((3,), (3, 4), (3, 4, 5))
        for T in (Float16, Float32, Float64)
            A = nans(T, s...)
            for a in A
                @test isnan(a)
            end
            @test size(A) == s
        end

        A = nans(s...)
        for a in A
            @test isnan(a)
        end
    end
end

@testset "readable secs feedback" begin
    using SpeedyWeatherInternals.Utils.TracableDates:
        CompoundPeriod, Day, Hour, Minute, Second, Millisecond
    @test readable_secs(123456) == CompoundPeriod(Day(1), Hour(10))
    @test readable_secs(12345) == CompoundPeriod(Hour(3), Minute(26))
    @test readable_secs(1234) == CompoundPeriod(Minute(20), Second(34))
    @test readable_secs(123) == CompoundPeriod(Minute(2), Second(3))
    @test readable_secs(12.3) == CompoundPeriod(Second(12), Millisecond(300))
    @test readable_secs(1.23) == CompoundPeriod(Second(1), Millisecond(230))
    @test readable_secs(0.123) == CompoundPeriod(Millisecond(120))
end

@testset "@maybe_jit without Reactant" begin
    # just test that it works without error, even with kwargs
    A = rand(10, 10)
    arch = SpeedyWeather.CPU()

    res = @maybe_jit arch sum(A)
    @test res == sum(A)

    res2 = @maybe_jit arch sum(A; dims = 1)
    @test res2 == sum(A; dims = 1)
end
