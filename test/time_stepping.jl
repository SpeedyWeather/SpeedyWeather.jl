# William (2009), MWR oscillation test case
# dF/dt = iωF 
F(x::Complex{T},ω::T) where T = im*ω*x
function F(L::LowerTriangularMatrix{Complex{T}},ω::T) where T
    tend = LowerTriangularMatrix{Complex{T}}(undef,size(L)...)
    for lm in SpeedyWeather.eachharmonic(L,tend)
        tend[lm] = F(L[lm],ω)
    end
    return tend
end

@testset "Leapfrog oscillation" begin

    ω = 1.0             # frequency
    Δt = 2π/100         # time step 
    n_rotations = 1     # times around the circle
    n_timesteps = round(Int,2π*n_rotations/(ω*Δt))

    # loop over different precisions
    @testset for NF in (Float16,Float32,Float64)
        P = Parameters{SpeedyWeather.BarotropicModel}(NF=NF)
        C = DynamicsConstants(P)

        # INITIAL CONDITIONS
        lmax,mmax = 3,3
        X_old = ones(LowerTriangularMatrix{Complex{NF}},3,3)
        X_new = copy(X_old)
        
        # store only 1 of the 3x3 values (all the same) per time step
        X_out = zeros(Complex{NF},n_timesteps+1)

        # exact 2nd leapfrog step
        X_new .*= exp(im*ω*Δt)
        X_out[1] = X_old[1,1]      # store initial conditions

        # leapfrog forward
        for i in 2:n_timesteps+1
            # always evaluate F with lf = 2
            lf = 2
            SpeedyWeather.leapfrog!(X_old,X_new,F(X_new,NF(ω)),NF(2Δt),lf,C)
            X_out[i] = X_old[1,1]
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
    n_timesteps = round(Int,2π*n_rotations/(ω*Δt))

    # loop over different precisions
    @testset for NF in (Float16,Float32,Float64)
        P = Parameters{SpeedyWeather.BarotropicModel}(NF=NF)
        C = DynamicsConstants(P)

        # INITIAL CONDITIONS
        lmax,mmax = 3,3
        X_old = ones(LowerTriangularMatrix{Complex{NF}},3,3)
        X_new = copy(X_old)
        
        # store only 1 of the 3x3 values (all the same) per time step
        X_out = zeros(Complex{NF},n_timesteps+1)

        # exact 2nd leapfrog step
        X_new .*= exp(im*ω*Δt)
        X_out[1] = X_old[1,1]      # store initial conditions

        # leapfrog forward
        for i in 2:n_timesteps+1
            # always evaluate F with lf = 2
            lf = 2
            SpeedyWeather.leapfrog!(X_old,X_new,F(X_new,NF(ω)),NF(2Δt),lf,C)
            X_out[i] = X_old[1,1]
        end

        # magnitude at last time step < 1 for stability
        M_RAW = abs(X_out[end])
        @test M_RAW < 1

        # CHECK THAT NO WILLIAM'S FILTER IS WORSE
        P = Parameters{SpeedyWeather.BarotropicModel}(NF=NF,williams_filter=1)     # Robert's filter only
        C = DynamicsConstants(P)

        # INITIAL CONDITIONS
        lmax,mmax = 3,3
        X_old = ones(LowerTriangularMatrix{Complex{NF}},3,3)
        X_new = copy(X_old)
        
        # store only 1 of the 3x3 values (all the same) per time step
        X_out = zeros(Complex{NF},n_timesteps+1)

        # exact 2nd leapfrog step
        X_new .*= exp(im*ω*Δt)
        X_out[1] = X_old[1,1]      # store initial conditions

        # leapfrog forward
        for i in 2:n_timesteps+1
            # always evaluate F with lf = 2
            lf = 2
            SpeedyWeather.leapfrog!(X_old,X_new,F(X_new,NF(ω)),NF(2Δt),lf,C)
            X_out[i] = X_old[1,1]
        end

        M_Ronly = abs(X_out[end])
        @test M_Ronly < 1
        @test M_Ronly <= M_RAW
    end
end