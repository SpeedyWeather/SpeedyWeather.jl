@testset "Sec, min, hrs arguments" begin
    SG = SpectralGrid(trunc=42, nlayers=1)
    L1 = Leapfrog(SG, Δt_at_T31=30)
    L2 = Leapfrog(SG, Δt_at_T31=Second(30))
    @test L1.Δt == L2.Δt
    
    L3 = Leapfrog(SG, Δt_at_T31=Second(300))
    L4 = Leapfrog(SG, Δt_at_T31=Minute(5))
    @test L3.Δt == L4.Δt

    L4 = Leapfrog(SG, Δt_at_T31=Minute(60))
    L5 = Leapfrog(SG, Δt_at_T31=Hour(1))
    @test L4.Δt == L5.Δt

    # without adjustment
    L1 = Leapfrog(SG, Δt_at_T31=30, adjust_with_output=false)
    L2 = Leapfrog(SG, Δt_at_T31=Second(30), adjust_with_output=false)
    @test L1.Δt == L2.Δt
    
    L3 = Leapfrog(SG, Δt_at_T31=Second(300), adjust_with_output=false)
    L4 = Leapfrog(SG, Δt_at_T31=Minute(5), adjust_with_output=false)
    @test L3.Δt == L4.Δt

    L4 = Leapfrog(SG, Δt_at_T31=Minute(60), adjust_with_output=false)
    L5 = Leapfrog(SG, Δt_at_T31=Hour(1), adjust_with_output=false)
    @test L4.Δt == L5.Δt

    # clock tests
    c1 = SpeedyWeather.Clock()
    SG2 = SpectralGrid(trunc=31, nlayers=1)
    L6 = Leapfrog(SG2, Δt_at_T31=Hour(1), adjust_with_output=false)

    # set period
    SpeedyWeather.initialize!(c1, L6, Hour(10))
    @test c1.n_timesteps == 10
    SpeedyWeather.initialize!(c1, L6, Day(2))
    @test c1.n_timesteps == 48
    
    # set n_timesteps
    SpeedyWeather.initialize!(c1, L6, 10)
    @test c1.n_timesteps == 10
end

@testset "Set clock" begin
    spectral_grid = SpectralGrid(nlayers=1)
    time_stepping = Leapfrog(spectral_grid)
    Δt = time_stepping.Δt_at_T31

    # set n_timesteps
    clock = Clock()
    n_timesteps = 100
    initialize!(clock, time_stepping, n_timesteps)
    @test clock.n_timesteps == 100
    @test clock.period == Second(100*Δt)

    # set period
    clock = Clock()
    period = Day(10)
    initialize!(clock, time_stepping, period)
    @test clock.period == period
    @test clock.n_timesteps == ceil(Int, Millisecond(period).value/time_stepping.Δt_millisec.value)

    model = BarotropicModel(spectral_grid)
    simulation = initialize!(model)
    run!(simulation, steps=1)
    run!(simulation, period=Hour(1))
    @test_throws AssertionError run!(simulation, steps=1, period=Day(1))
end
    
import Dates

@testset "time conversions" begin
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
    @test Dates.coarserperiod(Year) == (Century, 100)
    @test_throws MethodError Week(Century(1))

    # Millenium
    @test convert(Year, Millenium(1)) == Year(1000)
    @test convert(Century, Millenium(1)) == Century(10)
    @test Millenium(1) == Year(1000)
    @test Millenium(1) == Century(10)
    @test Second(Millenium(1)) == Second(1000 * 365 * 24 * 60 * 60)
    @test Millisecond(Millenium(1)) == Millisecond(Second(Millenium(1)))
    @test Dates.coarserperiod(Century) == (Millenium, 10)
    @test_throws MethodError Week(Millenium(1))

    # Type promotion rules
    @test promote(Year(1), Century(1)) == (Year(1), Year(100))
    @test promote(Year(1), Millenium(1)) == (Year(1), Year(1000))
    @test promote(Century(1), Millenium(1)) == (Century(1), Century(10))
end