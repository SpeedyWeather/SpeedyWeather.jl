@testset "Sec, min, hrs arguments" begin
    SG = SpectralGrid(trunc=42, nlev=1)
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
    SG2 = SpectralGrid(trunc=31, nlev=1)
    L6 = Leapfrog(SG2, Δt_at_T31=Hour(1), adjust_with_output=false)

    SpeedyWeather.set_period!(c1, Hour(10))
    SpeedyWeather.initialize!(c1, L6)
    @test c1.n_timesteps == 10 

    SpeedyWeather.set_period!(c1, 10)   # assumed to be in days
    SpeedyWeather.initialize!(c1, L6)
    @test c1.n_timesteps == 24*10

    SpeedyWeather.set_period!(c1, 10.0) # also assumed to be in days
    SpeedyWeather.initialize!(c1, L6)
    @test c1.n_timesteps == 24*10 

    SpeedyWeather.set_period!(c1, Day(2))
    SpeedyWeather.initialize!(c1, L6)
    @test c1.n_timesteps == 48
end