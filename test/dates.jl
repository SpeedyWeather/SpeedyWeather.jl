@testset "Sec, min, hrs arguments" begin
    SG = SpectralGrid(trunc=42,nlev=1)
    L1 = Leapfrog(SG,Δt_at_T31=30)
    L2 = Leapfrog(SG,Δt_at_T31=Second(30))
    @test L1.Δt == L2.Δt
    
    L3 = Leapfrog(SG,Δt_at_T31=Second(300))
    L4 = Leapfrog(SG,Δt_at_T31=Minute(5))
    @test L3.Δt == L4.Δt

    L4 = Leapfrog(SG,Δt_at_T31=Minute(60))
    L5 = Leapfrog(SG,Δt_at_T31=Hour(1))
    @test L4.Δt == L5.Δt

    # without adjustment
    L1 = Leapfrog(SG,Δt_at_T31=30,adjust_with_output=false)
    L2 = Leapfrog(SG,Δt_at_T31=Second(30),adjust_with_output=false)
    @test L1.Δt == L2.Δt
    
    L3 = Leapfrog(SG,Δt_at_T31=Second(300),adjust_with_output=false)
    L4 = Leapfrog(SG,Δt_at_T31=Minute(5),adjust_with_output=false)
    @test L3.Δt == L4.Δt

    L4 = Leapfrog(SG,Δt_at_T31=Minute(60),adjust_with_output=false)
    L5 = Leapfrog(SG,Δt_at_T31=Hour(1),adjust_with_output=false)
    @test L4.Δt == L5.Δt
end