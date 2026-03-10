using SpeedyWeatherInternals.Utils

@testset "Increasing/decresing vectors" begin
    @test isincreasing(collect(1:10))
    @test isincreasing(sort(rand(10)))
    @test ~isincreasing(rand(Float32, 10))
    @test ~isincreasing(randn(10))
end

using Dates: CompoundPeriod, Day, Hour, Minute, Second, Millisecond, coarserperiod

@testset "readable secs feedback" begin
    @test readable_secs(123456) == CompoundPeriod(Day(1), Hour(10))
    @test readable_secs(12345) == CompoundPeriod(Hour(3), Minute(26))
    @test readable_secs(1234) == CompoundPeriod(Minute(20), Second(34))
    @test readable_secs(123) == CompoundPeriod(Minute(2), Second(3))
    @test readable_secs(12.3) == CompoundPeriod(Second(12), Millisecond(300))
    @test readable_secs(1.23) == CompoundPeriod(Second(1), Millisecond(230))
    @test readable_secs(0.123) == CompoundPeriod(Second(0), Millisecond(120))
end
