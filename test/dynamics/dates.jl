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
    
