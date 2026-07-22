using SpeedyWeatherInternals.Utils

@testset "Increasing/decresing vectors" begin
    @test isincreasing(collect(1:10))
    @test isincreasing(sort(rand(10)))
    @test ~isincreasing(rand(Float32, 10))
    @test ~isincreasing(randn(10))
end

@testset "readable secs feedback" begin
    using Dates: CompoundPeriod, Day, Hour, Minute, Second, Millisecond
    @test readable_secs(123456) == CompoundPeriod(Day(1), Hour(10))
    @test readable_secs(12345) == CompoundPeriod(Hour(3), Minute(26))
    @test readable_secs(1234) == CompoundPeriod(Minute(20), Second(34))
    @test readable_secs(123) == CompoundPeriod(Minute(2), Second(3))
    @test readable_secs(12.3) == CompoundPeriod(Second(12), Millisecond(300))
    @test readable_secs(1.23) == CompoundPeriod(Second(1), Millisecond(230))
    @test readable_secs(0.123) == CompoundPeriod(Millisecond(120))
end

@testset "set_scalar!" begin
    x = rand(Float32, 5)
    set_scalar!(x, 2, 42.0f0)
    @test x[2] == 42.0f0

    y = rand(ComplexF32, 4)
    set_scalar!(y, 1, zero(ComplexF32))
    @test y[1] == 0
end

@testset "@maybe_jit without Reactant" begin
    # just test that it works without error, even with kwargs
    A = rand(10, 10)
    arch = SpeedyWeatherInternals.Architectures.CPU()

    res = @maybe_jit arch sum(A)
    @test res == sum(A)

    res2 = @maybe_jit arch sum(A; dims = 1)
    @test res2 == sum(A; dims = 1)
end
