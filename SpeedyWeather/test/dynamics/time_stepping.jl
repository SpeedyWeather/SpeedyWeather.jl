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
        spectral_grid = SpectralGrid(; NF, truncation = 6, nlayers = 1)
        L = Leapfrog(spectral_grid, adjust_with_output = false, robert_filter = 0.05, williams_filter = 0.51)
        model = BarotropicModel(spectral_grid; time_stepping = L)
        simulation = initialize!(model)
        (; clock) = simulation.variables.prognostic
        L.Δt = Δt
        clock.step_counter = 2
        clock.time_step_counter = 2

        # INITIAL CONDITIONS
        for space in (:grid, :spectrum)
            if space == :grid
                X = ones(Complex{NF}, spectral_grid.grid, 2)
                dX = zeros(Complex{NF}, spectral_grid.grid, 1)
            elseif space == :spectrum
                X = ones(Complex{NF}, spectral_grid.spectrum, 2)
                dX = zeros(Complex{NF}, spectral_grid.spectrum, 1)
            end

            X[:, 2] .*= exp(im * ω * Δt)    # exact 2nd leapfrog step

            # leapfrog forward
            for i in 1:(n_time_steps - 1)
                dX.data .= F.(X[:, 2], NF(ω))
                SpeedyWeather.update_prognostic!(X, dX, clock, L, nothing, model)
            end

            # absolute error to exact result 1+0i
            error = abs.(X[:, 2] .- 1)
            # @info Leapfrog, error[1], abs(X[1])
            @test all(error .< 1.0e-2)
            @test all(abs.(X) .<= 1)         # stable integration?

            # long term stability
            for i in 1:(10 * n_time_steps)
                dX.data .= F.(X[:, 2], NF(ω))
                SpeedyWeather.update_prognostic!(X, dX, clock, L, nothing, model)
            end

            @test all(abs.(X) .<= 1)         # still stable?
        end
    end
end

@testset "Leapfrog spinup" begin

    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)

    # disable RAW filters
    time_stepping = Leapfrog(spectral_grid, adjust_with_output = false, robert_filter = 0, williams_filter = 1)
    model = BarotropicModel(spectral_grid; time_stepping)
    simulation = initialize!(model)
    (; clock) = simulation.variables.prognostic
    (; NF) = spectral_grid

    Δt = 1
    time_stepping.Δt = Δt

    # initial conditions in step 1, 0 in step 2
    X = simulation.variables.prognostic.vorticity
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
    @test all(X2 .== X0 .+ (Δt // 2) * dX1)
    @test all(X1 .== X0)         # previous time step still initial conditions

    # X2old = deepcopy(X2)

    # new time step, new random tendencies
    dX .= rand(Complex{NF}, spectral_grid.spectrum, 1, 1)
    dX .= im * imag.(dX)          # make them imaginary only to better tell them apart
    SpeedyWeather.update_prognostic!(X, dX, clock, time_stepping, model.implicit, model)
    SpeedyWeather.time_step!(clock, time_stepping)

    # this Leapfrog step does count!
    @test clock.step_counter == 2
    @test clock.time_step_counter == 1

    # with Δt step size
    @test all(X2 .== X1 .+ Δt * dX1)
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
    @test all(X2 .== X0 .+ 2Δt * dX1)

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
    @test all(X2 .== X1old .+ 2Δt * dX1)
end

@testset "NCycleLorenz oscillation" begin

    ω = 1.0             # frequency
    Δt = 2π / 192        # time step, choose 120 as both 3 and 4 are divisors
    n_rotations = 1     # times around the circle
    n_time_steps = round(Int, 2π * n_rotations / (ω * Δt))
    scale = 1

    # loop over different precisions
    @testset for NF in (Float32, Float64)
        @testset for Variant in (
                SpeedyWeather.NCycleLorenzA,
                SpeedyWeather.NCycleLorenzB,
                SpeedyWeather.NCycleLorenzAB,
                SpeedyWeather.NCycleLorenzABBA,
            )
            @testset for steps in (3, 4)
                spectral_grid = SpectralGrid(; NF, truncation = 6, nlayers = 1)
                L = NCycleLorenz(spectral_grid; steps = steps, variant = Variant(), adjust_with_output = false)
                model = BarotropicModel(spectral_grid; time_stepping = L)
                simulation = initialize!(model)
                (; clock) = simulation.variables.prognostic
                L.Δt = Δt

                # INITIAL CONDITIONS
                X = ones(LowerTriangularArray{Complex{NF}}, spectral_grid.spectrum, 1)
                dX = zeros(LowerTriangularArray{Complex{NF}}, spectral_grid.spectrum, 2)

                for i in 1:n_time_steps
                    dX[:, 1] .= F.(X, NF(ω))        # tendency in the first step (second is multistep averaged tendency)
                    SpeedyWeather.update_prognostic!(X, dX, clock, L, nothing, model)
                    SpeedyWeather.time_step!(clock, L)
                end

                # absolute error to exact result 1+0i
                error = abs.(X .- 1)
                # @info (steps, Variant, error[1], abs(X[1]))
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
        w3 = [SpeedyWeather.weight_coefficient(NF, SpeedyWeather.NCycleLorenzABBA(), i - 1, 3) for i in 1:12]
        @test w3 == NF[1.0, 1.5, 3.0, 1.0, 3.0, 1.5, 1.0, 3.0, 1.5, 1.0, 1.5, 3.0]

        w4 = [SpeedyWeather.weight_coefficient(NF, SpeedyWeather.NCycleLorenzABBA(), i - 1, 4) for i in 1:16]
        @test w4 == NF[
            1.0, 1 + 1 / 3, 2.0, 4.0,
            1.0, 4.0, 2.0, 1 + 1 / 3,
            1.0, 4.0, 2.0, 1 + 1 / 3,
            1.0, 1 + 1 / 3, 2.0, 4.0,
        ]
    end
end

@testset "Default time step inside kernels" begin
    # time steppers are adapted to a core struct inside kernels (GPU, Reactant), which has to
    # return the same default time step, e.g. 2Δt for leapfrog used by the surface flux damping
    @testset for TS in (Leapfrog, NCycleLorenz)
        time_stepping = TS(SpectralGrid())
        adapted = SpeedyWeather.Adapt.adapt(Array, time_stepping)
        @test typeof(adapted) != typeof(time_stepping)
        @test SpeedyWeather.default_time_step(adapted) == SpeedyWeather.default_time_step(time_stepping)
    end
end

@testset "Set timestep manually" begin
    @testset for TS in (Leapfrog, NCycleLorenz)
        @testset for truncation in (32, 64, 128)
            @testset for Δt in (Minute(10), Minute(20))
                spectral_grid = SpectralGrid(; truncation)
                time_stepping = TS(spectral_grid)
                set!(time_stepping, Δt = Δt)
                @test time_stepping.Δt == Second(Δt).value
                @test time_stepping.Δt_millisec == Millisecond(Second(time_stepping.Δt))

                # Δt = Δt_at_T32 * (32/truncation) (see get_Δt_millisec), so Δt_at_T32 is the
                # requested Δt divided by that factor, i.e. LONGER than Δt for truncation > 32
                factor = SpeedyWeather.resolution_factor(truncation, SpeedyWeather.DEFAULT_RADIUS)
                @test time_stepping.Δt_at_T32.value ≈ Second(Δt).value / factor rtol = 0.05
                @test time_stepping.Δt_at_T32 >= Second(Δt)   # never shorter than Δt itself
            end
        end
    end
end

@testset "Bit reproducibility with NCycleLorenz" begin
    @testset for steps in (3, 4)
        @testset for Variant in (
                SpeedyWeather.NCycleLorenzA,
                SpeedyWeather.NCycleLorenzB,
                SpeedyWeather.NCycleLorenzAB,
                SpeedyWeather.NCycleLorenzABBA,
            )
            s = 4      # run longer? As testing for approximate below, s can't be too large

            spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
            time_stepping = NCycleLorenz(spectral_grid; steps, variant = Variant())
            planet = Earth(spectral_grid, radius = 2^22)  # use radius that is power of 2 to avoid rounding errors in scaling

            ic = RandomVelocity(spectral_grid, seed = 1234)
            model = BarotropicModel(spectral_grid; time_stepping, initial_conditions = ic)
            model.feedback.verbose = false

            simulation = initialize!(model)
            run!(simulation, steps = 2 * s * 8 * steps)

            vor_restarted = deepcopy(simulation.variables.prognostic.vorticity)
            time_restarted = simulation.variables.prognostic.clock.time

            # do a new simulation from same model
            simulation = initialize!(model)
            run!(simulation, steps = 2 * s * 8 * steps)
            @test vor_restarted == simulation.variables.prognostic.vorticity
            @test time_restarted == simulation.variables.prognostic.clock.time

            # check bit reproducibility of scaling
            SpeedyWeather.scale_prognostic!(simulation.variables, planet.radius)
            @test vor_restarted != simulation.variables.prognostic.vorticity
            SpeedyWeather.unscale!(simulation.variables)
            @test vor_restarted == simulation.variables.prognostic.vorticity

            # with restart half way
            simulation = initialize!(model)
            run!(simulation, steps = s * 8 * steps)
            run!(simulation, steps = s * 8steps)

            # this test is only approximate as bit reproducibility is close but not perfect
            # not sure exactly why, needs further investigation if deemed important
            @test all(isapprox.(vor_restarted, simulation.variables.prognostic.vorticity; atol = 1.0e-5, rtol = 1.0e-5))
            @test time_restarted == simulation.variables.prognostic.clock.time
        end
    end
end

@testset "scale_tendencies! / unscale_tendencies! with fused tendency parents" begin
    # Regression test: the dycore's fused tendency members (e.g. grid u/v/temperature/…)
    # are views into a shared parent buffer. `scale_tendencies!` must scale each parent
    # buffer exactly once — not once per member view (which corrupts the data under Reactant,
    # and could double-scale) — and `scale ∘ unscale` must round-trip to the identity.
    @testset for Model in (BarotropicModel, ShallowWaterModel, PrimitiveWetModel)
        nlayers = Model == PrimitiveWetModel ? 4 : 1
        spectral_grid = SpectralGrid(; truncation = 9, nlayers)
        model = Model(spectral_grid)
        vars = Variables(model)
        r = 3.0f0

        # known nonzero tendency state across the whole fused grid_tendencies parent buffer
        grid_parent = parent(vars.fused.grid_tendencies)
        fill!(grid_parent.data, 2)
        original = copy(Array(grid_parent.data))

        vars.prognostic.scale[] = r                     # scale factor read by (un)scale_tendencies!
        SpeedyWeather.scale_tendencies!(vars, model)
        # scaled exactly once: 2*r everywhere (not 2*r² from double-scaling, not corrupted)
        @test Array(grid_parent.data) ≈ original .* r

        SpeedyWeather.unscale_tendencies!(vars)
        @test Array(grid_parent.data) ≈ original         # round-trips to the identity
    end
end

@testset "unscale! leaves the diagnosed vertical velocity scaled" begin
    # `unscale!` undoes the radius scaling of vor, div but not of w = radius*σ̇, which the
    # dycore uses scaled. ParticleAdvection3D accounts for it with a 1/radius in its vertical time step.
    spectral_grid = SpectralGrid(truncation = 9, nlayers = 4)
    model = PrimitiveWetModel(spectral_grid)
    vars = Variables(model)

    fill!(vars.dynamics.w.data, 2)
    vars.prognostic.scale[] = 3
    SpeedyWeather.unscale!(vars)
    @test all(Array(vars.dynamics.w.data) .== 2)
end

@testset "Single time steps match run! from a reset clock" begin
    # `initialize!(simulation; steps)` resets the clock, so the next steps are the leapfrog
    # start (Δt/2, Δt, then 2Δt) while the implicit solver still has 2Δt from the spin-up.
    # A loop over `time_step!(vars, time_stepping, model)` (as used for differentiation)
    # has to reinitialize the solver like `run!` does, or it integrates a different model.
    @testset for Model in (PrimitiveDryModel, PrimitiveWetModel)
        # two identical simulations, both spun up so that the implicit solver is at 2Δt
        function spun_up()
            spectral_grid = SpectralGrid(truncation = 22, nlayers = 4, NF = Float64)
            simulation = initialize!(Model(spectral_grid))
            run!(simulation, steps = 5)
            return simulation
        end

        steps = 4
        sim_run = spun_up()
        run!(sim_run; steps)

        sim_step = spun_up()
        initialize!(sim_step; steps)
        (; variables, model) = sim_step
        for _ in 1:steps
            SpeedyWeather.time_step!(variables, model.time_stepping, model)
            SpeedyWeather.time_step!(variables.prognostic.clock, model.time_stepping)
        end
        SpeedyWeather.unscale!(variables)

        @test sim_step.variables.prognostic.clock.time == sim_run.variables.prognostic.clock.time
        for var in (:vorticity, :divergence, :temperature, :pressure)
            @test sim_step.variables.prognostic[var] ≈ sim_run.variables.prognostic[var]
        end
    end
end

@testset "Implicit reference temperature profile: $reference_temperature" for reference_temperature in (
        GlobalMeanReferenceTemperature(offset = 0.0f0),
        GlobalMeanReferenceTemperature(offset = 20.0f0),
        IsothermalReferenceTemperature(temperature = 280.0f0),
        MaximumReferenceTemperature(margin = 30.0f0),
    )
    spectral_grid = SpectralGrid(truncation = 22, nlayers = 4)
    implicit = ImplicitPrimitiveEquation(spectral_grid; reference_temperature)
    model = PrimitiveWetModel(spectral_grid; implicit)
    simulation = initialize!(model)
    (; implicit) = model
    (; clock) = simulation.variables.prognostic
    σ_levels = Array(model.geometry.σ_levels_full)
    n = implicit.schedule.every_n_steps

    for (i, ΔT) in enumerate((0, 5))
        # new temperature: global mean 200 + 100σ + ΔT per layer, plus a wave with zero global mean
        set!(simulation, temperature = (λ, φ, σ) -> 200 + 100σ + ΔT + 10 * cosd(λ) * cosd(φ))
        SpeedyWeather.transform!(simulation.variables, model, initialize = true)   # computes the global mean

        clock.time_step_counter = i * n     # a time step on which the reference profile is recomputed
        SpeedyWeather.reinitialize!(implicit, model, simulation.variables)
        if reference_temperature isa GlobalMeanReferenceTemperature
            (; offset) = reference_temperature
            @test Array(implicit.temp_profile) ≈ 200 .+ 100σ_levels .+ ΔT .+ offset atol = 0.01
        elseif reference_temperature isa MaximumReferenceTemperature
            # isothermal at max(T) + margin
            temp_grid = Array(SpeedyWeather.get_prognostic_step(simulation.variables.grid.temperature, model.time_stepping, implicit))
            expected = maximum(temp_grid) + reference_temperature.margin
            @test all(==(expected), Array(implicit.temp_profile))
        else
            @test all(==(reference_temperature.temperature), Array(implicit.temp_profile))
        end
    end

    # a constant isothermal reference profile is only reinitialized on a new time step, not on the schedule
    @test SpeedyWeather.rescheduled(reference_temperature) == !(reference_temperature isa IsothermalReferenceTemperature)
end

@testset "Leapfrog: previous grid step matches previous spectral step" begin
    # during the Leapfrog start (Euler step, first leapfrog step) the previous spectral step
    # is not advanced, so the previous grid step must not be either
    spectral_grid = SpectralGrid(truncation = 22, nlayers = 4)
    model = PrimitiveDryModel(spectral_grid; time_stepping = Leapfrog(spectral_grid, robert_filter = 0))
    simulation = initialize!(model)
    initialize!(simulation; steps = 4)
    (; variables) = simulation
    for _ in 1:4
        SpeedyWeather.time_step!(variables, model.time_stepping, model)
        SpeedyWeather.time_step!(variables.prognostic.clock, model.time_stepping)
        temperature_grid = transform(get_step(variables.prognostic.temperature, 1), model.spectral_transform)
        @test get_step(variables.grid.temperature, 1) ≈ temperature_grid rtol = 1.0e-5
    end
end

@testset "NCycleLorenz implicit centering per substep" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
    implicit = ImplicitShallowWater(spectral_grid, centering = 0.5)
    clock = Clock()

    # default centerings for the 3-cycle variants A and B, implicit.centering otherwise
    for (variant, centering) in ((NCycleLorenzA(), (0.116, 1.039, 0.462)), (NCycleLorenzB(), (1.386, 0.468, -0.464)))
        L = NCycleLorenz(spectral_grid; steps = 3, variant)
        @test L.implicit_centering == Float32.(centering)
        for (i, θ) in enumerate(L.implicit_centering)
            clock.step_counter = i - 1
            @test SpeedyWeather.implicit_centering(L, implicit, clock) == θ
        end
    end

    for L in (NCycleLorenz(spectral_grid, steps = 4), NCycleLorenz(spectral_grid, variant = NCycleLorenzAB()))
        @test isnothing(L.implicit_centering)
        @test SpeedyWeather.implicit_centering(L, implicit, clock) == implicit.centering
    end

    L = NCycleLorenz(spectral_grid, steps = 3, implicit_centering = (0.5, 0.5, 0.5))
    @test SpeedyWeather.implicit_centering(L, implicit, clock) == 0.5f0
end

@testset "ShallowWater gravity waves with NCycleLorenz" begin
    # variance of η (excluding the mean) relative to the initial conditions
    function η_variance(η)
        (; l_indices, m_indices) = η.spectrum
        return sum((m == 1 ? 1 : 2) * abs2(η[lm]) for (lm, (l, m)) in enumerate(zip(l_indices, m_indices)) if l > 1)
    end

    variance_ratio = Dict()
    for (name, TS) in (("Leapfrog", Leapfrog), ("NCycleLorenzA", sg -> NCycleLorenz(sg)), ("NCycleLorenzB", sg -> NCycleLorenz(sg, variant = NCycleLorenzB())))
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
        initial_conditions = RandomWaves(spectral_grid, lmin = 5, lmax = 15)
        time_stepping = TS(spectral_grid)
        model = ShallowWaterModel(spectral_grid; initial_conditions, time_stepping)
        simulation = initialize!(model)
        variance0 = η_variance(get_step(simulation.variables.prognostic.η, 1))
        run!(simulation, period = Day(1))
        @test !simulation.model.feedback.nans_detected
        variance_ratio[name] = η_variance(get_step(simulation.variables.prognostic.η, 1)) / variance0
    end

    # NCycleLorenz used to damp gravity waves much stronger than Leapfrog
    @test variance_ratio["NCycleLorenzA"] > variance_ratio["Leapfrog"] / 2
    @test variance_ratio["NCycleLorenzB"] > variance_ratio["Leapfrog"] / 2
end
