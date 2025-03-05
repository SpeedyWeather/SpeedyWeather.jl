# Williams (2009), MWR oscillation test case
# dF/dt = iωF 
F(x::Complex{T}, ω::T) where T = im*ω*x
function F(L::LowerTriangularMatrix{Complex{T}}, ω::T) where T
    tend = LowerTriangularMatrix{Complex{T}}(undef, size(L, as=Matrix)...)
    for lm in SpeedyWeather.eachharmonic(L, tend)
        tend[lm] = F(L[lm], ω)
    end
    return tend
end

@testset "Leapfrog oscillation" begin

    ω = 1.0             # frequency
    Δt = 2π/100         # time step 
    n_rotations = 1     # times around the circle
    n_timesteps = round(Int, 2π*n_rotations/(ω*Δt))

    # loop over different precisions
    @testset for NF in (Float16, Float32, Float64)

        spectral_grid = SpectralGrid(; NF)
        L = Leapfrog(spectral_grid)

        # INITIAL CONDITIONS
        lmax, mmax = 3, 3
        X_old = ones(LowerTriangularMatrix{Complex{NF}}, 3, 3)
        X_new = copy(X_old)
        
        # store only 1 of the 3x3 values (all the same) per time step
        X_out = zeros(Complex{NF}, n_timesteps+1)

        # exact 2nd leapfrog step
        X_new .*= exp(im*ω*Δt)
        X_out[1] = X_old[1, 1]      # store initial conditions

        # leapfrog forward
        for i in 2:n_timesteps+1
            # always evaluate F with lf = 2
            lf = 2
            SpeedyWeather.leapfrog!(X_old, X_new, F(X_new, NF(ω)), NF(2Δt), lf, L)
            X_out[i] = X_old[1, 1]
        end

        # absolute error to exact result 1+0i
        error = abs(X_out[end]-1)
        @test error < 1e-2
    end
end

@testset "Leapfrog stability" begin

    # LONG TERM STABILITY
    ω = 1.0             # frequency
    Δt = 2π/100         # time step 
    n_rotations = 10
    n_timesteps = round(Int, 2π*n_rotations/(ω*Δt))

    # loop over different precisions
    @testset for NF in (Float16, Float32, Float64)

        spectral_grid = SpectralGrid(; NF)
        L = Leapfrog(spectral_grid)

        # INITIAL CONDITIONS
        lmax, mmax = 3, 3
        X_old = ones(LowerTriangularMatrix{Complex{NF}}, 3, 3)
        X_new = copy(X_old)
        
        # store only 1 of the 3x3 values (all the same) per time step
        X_out = zeros(Complex{NF}, n_timesteps+1)

        # exact 2nd leapfrog step
        X_new .*= exp(im*ω*Δt)
        X_out[1] = X_old[1, 1]      # store initial conditions

        # leapfrog forward
        for i in 2:n_timesteps+1
            # always evaluate F with lf = 2
            lf = 2
            SpeedyWeather.leapfrog!(X_old, X_new, F(X_new, NF(ω)), NF(2Δt), lf, L)
            X_out[i] = X_old[1, 1]
        end

        # magnitude at last time step < 1 for stability
        M_RAW = abs(X_out[end])
        @test M_RAW < 1

        # CHECK THAT NO WILLIAMS FILTER IS WORSE
        spectral_grid = SpectralGrid(; NF)
        L = Leapfrog(spectral_grid, williams_filter=1)

        # INITIAL CONDITIONS
        lmax, mmax = 3, 3
        X_old = ones(LowerTriangularMatrix{Complex{NF}}, 3, 3)
        X_new = copy(X_old)
        
        # store only 1 of the 3x3 values (all the same) per time step
        X_out = zeros(Complex{NF}, n_timesteps+1)

        # exact 2nd leapfrog step
        X_new .*= exp(im*ω*Δt)
        X_out[1] = X_old[1, 1]      # store initial conditions

        # leapfrog forward
        for i in 2:n_timesteps+1
            # always evaluate F with lf = 2
            lf = 2
            SpeedyWeather.leapfrog!(X_old, X_new, F(X_new, NF(ω)), NF(2Δt), lf, L)
            X_out[i] = X_old[1, 1]
        end

        M_Ronly = abs(X_out[end])
        @test M_Ronly < 1
        @test M_Ronly <= M_RAW
    end
end

@testset "Set timestep manually" begin
    @testset for trunc in (31, 63, 127)
        spectral_grid = SpectralGrid(; trunc)
        time_stepping = Leapfrog(spectral_grid)
        set!(time_stepping, Minute(10))
        @test time_stepping.Δt_sec == 10*60

        set!(time_stepping, Δt=Minute(20))
        @test time_stepping.Δt_sec == 20*60
    end
end

@testset "Set clock" begin
    spectral_grid = SpectralGrid()
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
end
    
