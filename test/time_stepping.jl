@testset "Leapfrog oscillation" begin
    
    # William (2009), MWR oscillation test case
    # dF/dt = iωF 
    F(x::Complex{T},ω::T) where T = im*ω*x
    ω = 1               # frequency
    Δt = 2π/100         # time step 
    n_rotations = 1     # times around the circle
    n_timesteps = round(Int,2π*n_rotations/(ω*Δt))

    # loop over different precisions
    for NF in (Float16,Float32,Float64)
        P = Parameters(NF=NF)
        C = Constants(P)

        # INITIAL CONDITIONS
        # with lmax x mmax x nleapfrog x nlev
        X = ones(Complex{NF},3,3,2,2)
        
        # store only 1 of the 3x3x2 values (all the same) per time step
        Xout = zeros(Complex{NF},n_timesteps+1)

        # exact 2nd leapfrog step 
        X[:,:,2,:] = X[:,:,1,:]*Complex{NF}(exp(im*ω*Δt))
        Xout[1] = X[1,1,1,1]      # store initial conditions

        # leapfrog forward
        lf = 2  # leapfrog index to be used for tendency
        for i in 2:n_timesteps+1
            # always evaluate F with lf=2
            SpeedyWeather.leapfrog!(X,F.(X[:,:,lf,:],NF(ω)),NF(2Δt),C)
            Xout[i] = X[1,1,1,1]
        end

        # absolute error to exact result 1+0i
        error = abs(Xout[end]-1)
        # println(error)
        @test error < 1e-2
    end

    # LONG TERM STABILITY
    n_rotations = 35
    n_timesteps = round(Int,2π*n_rotations/(ω*Δt))

    for NF in (Float16,Float32,Float64)
        P = Parameters(NF=NF)
        C = Constants(P)

        # INITIAL CONDITIONS
        # with lmax x mmax x nleapfrog
        X = ones(Complex{NF},3,3,2)
        
        # store only 1 of the 3x3 values (all the same) per time step
        Xout = zeros(Complex{NF},n_timesteps+1)

        # exact 2nd leapfrog step 
        X[:,:,2] = X[:,:,1]*Complex{NF}(exp(im*ω*Δt))
        Xout[1] = X[1,1,1]      # store initial conditions

        # leapfrog forward
        lf = 2  # leapfrog index to be used for tendency
        for i in 2:n_timesteps+1
            SpeedyWeather.leapfrog!(X,F.(X[:,:,lf],NF(ω)),NF(2Δt),C)
            Xout[i] = X[1,1,1]
        end
        #println(Xout)
        # absolute error to exact result 1+0i
        error = abs(Xout[end]-1)

        # magnitude at last time step < 1 for stability
        M_RAW = abs(Xout[end])
        println("Error after " * string(n_timesteps) * " timesteps for NF " * string(NF) * " is: ", error)

        @test M_RAW < 1

        # CHECK THAT NO WILLIAM'S FILTER IS WORSE
        P = Parameters(NF=NF,williams_filter=1)     # Robert's filter only
        C = Constants(P)

        # INITIAL CONDITIONS
        # with 3x3 some spectral coeffs, x2 for leapfrog dimension
        X = ones(Complex{NF},3,3,2)
        
        # store only 1 of the 3x3 values (all the same) per time step
        Xout = zeros(Complex{NF},n_timesteps+1)

        # exact 2nd leapfrog step 
        X[:,:,2] = X[:,:,1]*Complex{NF}(exp(im*ω*Δt))
        Xout[1] = X[1,1,1]      # store initial conditions

        # leapfrog forward
        for i in 2:n_timesteps+1
            SpeedyWeather.leapfrog!(X,F.(X[:,:,lf],NF(ω)),NF(2Δt),C)
            Xout[i] = X[1,1,1]
        end

        M_Ronly = abs(Xout[end])

        @test M_Ronly <= M_RAW
    end
end