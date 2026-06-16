# Williams (2009), MWR oscillation test case
# dF/dt = iωF
F(x, ω) = im * ω * x

@testset "Leapfrog oscillation" begin

    ω = 1.0             # frequency
    Δt = 2π / 192       # time step
    n_rotations = 1     # times around the circle
    n_time_steps = round(Int, 2π * n_rotations / (ω * Δt))

    # loop over different precisions
    @testset for NF in (Float32, Float64)
        spectral_grid = SpectralGrid(; NF, trunc=5, nlayers = 1)
        L = Leapfrog(spectral_grid, adjust_with_output=false, robert_filter=0.05, williams_filter=0.51)
        model = BarotropicModel(spectral_grid; time_stepping=L)
        simulation = initialize!(model)
        (; clock) = simulation.variables.prognostic
        L.Δt = Δt
        clock.step_counter = 2
        clock.time_step_counter = 2

        # INITIAL CONDITIONS
        for space in (:grid, :spectrum)
            if space == :grid
                X  = ones( Complex{NF}, spectral_grid.grid, 2)
                dX = zeros(Complex{NF}, spectral_grid.grid, 1)
            elseif space == :spectrum
                X  = ones( Complex{NF}, spectral_grid.spectrum, 2)
                dX = zeros(Complex{NF}, spectral_grid.spectrum, 1)
            end

            X[:, 2] .*= exp(im * ω * Δt)    # exact 2nd leapfrog step

            # leapfrog forward
            for i in 1:n_time_steps-1
                dX.data .= F.(X[:, 2], NF(ω))
                SpeedyWeather.update_prognostic!(X, dX, clock, L, nothing, model)
            end

            # absolute error to exact result 1+0i
            error = abs.(X[:, 2] .- 1)
            @info Leapfrog, error[1], abs(X[1])
            @test all(error .< 1.0e-2)
            @test all(abs.(X) .<= 1)         # stable integration?

            # long term stability
            for i in 1:10*n_time_steps
                dX.data .= F.(X[:, 2], NF(ω))
                SpeedyWeather.update_prognostic!(X, dX, clock, L, nothing, model)
            end

            @test all(abs.(X) .<= 1)         # still stable?
        end
    end
end

@testset "Leapfrog spinup" begin

    spectral_grid = SpectralGrid(trunc=5, nlayers = 1)

    # disable RAW filters
    time_stepping = Leapfrog(spectral_grid, adjust_with_output=false, robert_filter=0, williams_filter=1)
    model = BarotropicModel(spectral_grid; time_stepping)
    simulation = initialize!(model)
    (; clock) = simulation.variables.prognostic
    (; NF) = spectral_grid

    Δt = 1
    time_stepping.Δt = Δt

    # initial conditions in step 1, 0 in step 2
    X  = simulation.variables.prognostic.vorticity
    X1 = get_step(X, 1)
    X2 = get_step(X, 2)

    X0 = rand(Complex{NF}, spectral_grid.spectrum, 1)      # initial conditions
    X1 .= X0        # set the first time step
    X2 .= 0         # 2nd step is always 0 at start
    @test all(X1 .!= X2)

    steps = 10
    initialize!(clock, time_stepping, steps)
    transform!(simulation.variables, model, initialize = true)
    
    # test that initial conditions have been copied to step 2
    @test all(X2 .== X1)

    # random tendencies
    dX = rand(Complex{NF}, spectral_grid.spectrum, 1, 1)
    dX .= real.(dX)             # make them real only to better tell them apart
    dX1 = get_step(dX, 1)
    SpeedyWeather.update_prognostic!(X, dX, clock, time_stepping, model.implicit, model)
    SpeedyWeather.time_step!(clock, time_stepping)

    # this Euler step does not count as time step but as step
    @test clock.step_counter == 1
    @test clock.time_step_counter == 0

    # do Euler step manually and compare
    @test all(X2 .== X0 .+ (Δt // 2)*dX1)
    @test all(X1 .== X0)         # previous time step still initial conditions

    # X2old = deepcopy(X2)

    # new time step, new random tendencies
    dX .= rand(Complex{NF}, spectral_grid.spectrum, 1, 1)
    dX .= im*imag.(dX)          # make them imaginary only to better tell them apart
    SpeedyWeather.update_prognostic!(X, dX, clock, time_stepping, model.implicit, model)
    SpeedyWeather.time_step!(clock, time_stepping)

    # this Leapfrog step does count!
    @test clock.step_counter == 2
    @test clock.time_step_counter == 1

    # with Δt step size
    @test all(X2 .== X1 .+ Δt*dX1)
    @test all(X1 .== X0)       # first step still unchanged

    X2old = deepcopy(X2)

    # new time step, new random tendencies
    dX .= rand(Complex{NF}, spectral_grid.spectrum, 1, 1)
    SpeedyWeather.update_prognostic!(X, dX, clock, time_stepping, model.implicit, model)
    SpeedyWeather.time_step!(clock, time_stepping)

    # this Leapfrog step does count!
    @test clock.step_counter == 3
    @test clock.time_step_counter == 2

    # with 2Δt step size
    @test X1 == X2old
    @test all(X2 .== X0 .+ 2Δt*dX1)

    X2old = deepcopy(X2)
    X1old = deepcopy(X1)

    # new time step, new random tendencies
    dX .= rand(Complex{NF}, spectral_grid.spectrum, 1, 1)
    SpeedyWeather.update_prognostic!(X, dX, clock, time_stepping, model.implicit, model)
    SpeedyWeather.time_step!(clock, time_stepping)

    # this Leapfrog step does count!
    @test clock.step_counter == 4
    @test clock.time_step_counter == 3

    # with 2Δt step size
    @test X1 == X2old
    @test all(X2 .== X1old .+ 2Δt*dX1)
end

@testset "NCycleLorenz oscillation" begin

    ω = 1.0             # frequency
    Δt = 2π / 192        # time step, choose 120 as both 3 and 4 are divisors
    n_rotations = 1     # times around the circle
    n_time_steps = round(Int, 2π * n_rotations / (ω * Δt))

    # loop over different precisions
    @testset for NF in (Float32, Float64)
        @testset for Variant in (SpeedyWeather.NCycleLorenzA,
                                    SpeedyWeather.NCycleLorenzB,
                                    SpeedyWeather.NCycleLorenzAB,
                                    SpeedyWeather.NCycleLorenzABBA,
                                    )
            @testset for steps in (3, 4)
                spectral_grid = SpectralGrid(; NF, trunc=5, nlayers = 1)
                L = NCycleLorenz(spectral_grid; steps=steps, variant=Variant(), adjust_with_output=false)
                model = BarotropicModel(spectral_grid; time_stepping=L)
                simulation = initialize!(model)
                (; clock) = simulation.variables.prognostic
                L.Δt = Δt

                # INITIAL CONDITIONS
                X  = ones( LowerTriangularArray{Complex{NF}}, spectral_grid.spectrum, 1)
                dX = zeros(LowerTriangularArray{Complex{NF}}, spectral_grid.spectrum, 2)

                for i in 1:n_time_steps
                    dX[:, 1] .= F.(X, NF(ω))        # tendency in the first step (second is multistep averaged tendency)
                    SpeedyWeather.update_prognostic!(X, dX, clock, L, nothing, model)
                    SpeedyWeather.time_step!(clock, L)
                end

                # absolute error to exact result 1+0i
                error = abs.(X .- 1)
                @info (steps, Variant, error[1], abs(X[1]))
                if steps == 3
                    @test all(error .< 1.0e-2)
                else                            
                    @test all(error .< 1.0e-3)      # more steps, higher order, lower error
                end
                @test all(abs.(X) .<= 1)
            end
        end
    end
end

@testset "NCycleLorenz: weight coefficients in cycle" begin
    @testset for NF in (Float32, Float64)
        w3 = [SpeedyWeather.weight_coefficient(NF, SpeedyWeather.NCycleLorenzABBA(), i-1, 3) for i in 1:12]
        @test w3 == NF[1.0, 1.5, 3.0, 1.0, 3.0, 1.5, 1.0, 3.0, 1.5, 1.0, 1.5, 3.0]

        w4 = [SpeedyWeather.weight_coefficient(NF, SpeedyWeather.NCycleLorenzABBA(), i-1, 4) for i in 1:16]
        @test w4 == NF[1.0, 1 + 1/3, 2.0, 4.0,
                        1.0, 4.0, 2.0, 1 + 1/3,
                        1.0, 4.0, 2.0, 1 + 1/3,
                        1.0, 1 + 1/3, 2.0, 4.0]
    end
end

@testset "Set timestep manually" begin
    @testset for TS in (Leapfrog, NCycleLorenz)
        @testset for trunc in (31, 63, 127)
            @testset for Δt in (Minute(10), Minute(20))
                spectral_grid = SpectralGrid(; trunc)
                time_stepping = TS(spectral_grid)
                set!(time_stepping, Δt=Δt)
                @test time_stepping.Δt_sec == Minute(Δt).value * 60
                @test time_stepping.Δt ≈ time_stepping.Δt_sec / SpeedyWeather.DEFAULT_RADIUS
                @test time_stepping.Δt_millisec == Millisecond(Second(time_stepping.Δt_sec))
                @test time_stepping.Δt_at_T31 == Second(Second(Δt).value / ((trunc + 1) / (SpeedyWeather.DEFAULT_TRUNC + 1)))
            end
        end
    end
end

@testset "Bit reproducibility with NCycleLorenz" begin
    @testset for steps in (3, 4)
        @testset for Variant in (SpeedyWeather.NCycleLorenzA,
                            SpeedyWeather.NCycleLorenzB,
                            SpeedyWeather.NCycleLorenzAB,
                            SpeedyWeather.NCycleLorenzABBA,
                            )
            s = 4      # run longer? As testing for approximate below, s can't be too large

            spectral_grid = SpectralGrid(nlayers = 1)
            time_stepping = NCycleLorenz(spectral_grid; steps, variant = Variant())
            planet = Earth(spectral_grid, radius = 2^22)  # use radius that is power of 2 to avoid rounding errors in scaling

            ic = RandomVelocity(spectral_grid, seed = 1234)
            model = BarotropicModel(spectral_grid; time_stepping, initial_conditions = ic)

            simulation = initialize!(model)
            run!(simulation, steps = 2*s*8*steps)

            vor_restarted = deepcopy(simulation.variables.prognostic.vorticity)
            time_restarted = simulation.variables.prognostic.clock.time

            # do a new simulation from same model
            simulation = initialize!(model)
            run!(simulation, steps = 2*s*8*steps)
            @test vor_restarted == simulation.variables.prognostic.vorticity
            @test time_restarted == simulation.variables.prognostic.clock.time

            # check bit reproducibility of scaling
            SpeedyWeather.scale_prognostic!(simulation.variables, planet.radius)
            @test vor_restarted != simulation.variables.prognostic.vorticity
            SpeedyWeather.unscale!(simulation.variables)
            @test vor_restarted == simulation.variables.prognostic.vorticity

            # with restart half way
            simulation = initialize!(model)
            run!(simulation, steps = s*8*steps)
            run!(simulation, steps = s*8steps)

            # this test is only approximate as bit reproducibility is close but not perfect
            # not sure exactly why, needs further investigation if deemed important
            @test all(vor_restarted .≈  simulation.variables.prognostic.vorticity)
            @test time_restarted == simulation.variables.prognostic.clock.time
        end
    end
end
