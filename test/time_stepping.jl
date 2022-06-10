using StochasticRounding

@testset "Leapfrog oscillation" begin
    
    # William (2009), MWR oscillation test case
    # dF/dt = iωF 
    F(x::Complex{T},ω::T) where T = im*ω*x
    ω = 1.143               # frequency
    Δt = 2π/101        # time step 
    n_rotations = 500.5     # times around the circle
    n_timesteps = round(Int,2π*n_rotations/(ω*Δt))

    # loop over different precisions
    for NF in (Float16sr,Float16,Float32,Float64)
        P = Parameters(NF=NF)
        C = Constants(P)

        # INITIAL CONDITIONS
        # with lmax x mmax x nleapfrog
        X = ones(Complex{NF},3,3,2) #.* Complex{NF}(2.1+0im) #Modify initial conditions slightly
        
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
            
            if i % n_rotations == 0
               println(Xout[i], ' ', abs(Xout[i]))
            end 

        end
        #println(Xout)
        # Error in the magnitude - ideally should always be = 1 if on unit circle
        error = abs(Xout[end]) - 1.0
        println("Error after " * string(n_timesteps) * " timesteps for NF " * string(NF) * " is: ", error)

    end
end