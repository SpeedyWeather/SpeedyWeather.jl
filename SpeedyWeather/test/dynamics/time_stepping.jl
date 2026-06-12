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
        w3 = [SpeedyWeather.weight_coefficient(NF, SpeedyWeather.NCycleLorenzABBA(), i, 3) for i in 0:11]
        @test w3 == NF[1.0, 1.5, 3.0, 1.0, 3.0, 1.5, 1.0, 3.0, 1.5, 1.0, 1.5, 3.0]

        w4 = [SpeedyWeather.weight_coefficient(NF, SpeedyWeather.NCycleLorenzABBA(), i, 4) for i in 0:15]
        @test w4 == NF[1.0, 1 + 1/3, 2.0, 4.0,
                        1.0, 4.0, 2.0, 1 + 1/3,
                        1.0, 4.0, 2.0, 1 + 1/3,
                        1.0, 1 + 1/3, 2.0, 4.0]
    end
end

@testset "Set timestep manually" begin
    @testset for trunc in (31, 63, 127)
        spectral_grid = SpectralGrid(; trunc)
        time_stepping = Leapfrog(spectral_grid)
        set!(time_stepping, Minute(10))
        @test time_stepping.Δt_sec == 10 * 60

        set!(time_stepping, Δt = Minute(20))
        @test time_stepping.Δt_sec == 20 * 60
    end
end

@testset "Bit reproducibility" begin
    spectral_grid = SpectralGrid(nlayers = 1)
    leapfrog = Leapfrog(spectral_grid; start_with_euler = true, continue_with_leapfrog = true)
    planet = Earth(spectral_grid, radius = 2^22)  # use radius that is power of 2 to avoid rounding errors in scaling

    ic = RandomVelocity(spectral_grid, seed = 1234)
    model = BarotropicModel(spectral_grid, time_stepping = leapfrog, initial_conditions = ic)

    simulation = initialize!(model)
    @test leapfrog.first_step_euler == true
    run!(simulation, steps = 10)
    @test leapfrog.first_step_euler == false

    vor_restarted = deepcopy(simulation.variables.prognostic.vorticity)
    time_restarted = simulation.variables.prognostic.clock.time

    # do a new simulation from same model
    simulation = initialize!(model)
    @test leapfrog.first_step_euler == true
    run!(simulation, steps = 10)
    @test leapfrog.first_step_euler == false
    @test vor_restarted == simulation.variables.prognostic.vorticity
    @test time_restarted == simulation.variables.prognostic.clock.time

    # check bit reproducibility of scaling
    SpeedyWeather.scale_prognostic!(simulation.variables, planet.radius)
    @test vor_restarted != simulation.variables.prognostic.vorticity
    SpeedyWeather.unscale!(simulation.variables)
    @test vor_restarted == simulation.variables.prognostic.vorticity

    # with restart half way
    simulation = initialize!(model)
    @test model.time_stepping.first_step_euler == true
    run!(simulation, steps = 5)
    @test model.time_stepping.first_step_euler == false
    run!(simulation, steps = 5)

    # this test is flagged as "broken" as bit reproducibility is close but not perfect
    # not sure exactly why, needs further investigation if deemed important
    @test_broken vor_restarted == simulation.variables.prognostic.vorticity
    @test time_restarted == simulation.variables.prognostic.clock.time
end
