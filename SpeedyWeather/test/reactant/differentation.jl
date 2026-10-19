# Differentiability with Reactant: Enzyme reverse mode through one time step, and through several
# time steps in one traced loop (`time_steps!`) with its checkpointing options. Enzyme-MLIR can only differentiate KernelAbstractions kernels that are
# raised to StableHLO, hence `raise = true`. The Reactant gradients are compared against Enzyme on
# the CPU for the same state.
# The models use the SpectralTransform (on Reactant the default, Fourier transform as truncated DFT) unless a
# `Transform` is given. CI runs this for the BarotropicModel (one time step with the SpectralTransform and with the
# MatrixSpectralTransform, several time steps) and, without the comparison against Enzyme on the CPU (too slow to
# compile), the gradient of several time steps of the PrimitiveWetModel; `differentiation_primitive_wet.jl` uses the
# functions defined here for the full comparison of the PrimitiveWetModel.
#
# Needs Reactant >= v0.2.294: raising kernels that store a complex constant (e.g. `div[I] = 0`) fails before
# 0.2.293 ("cannot raise op to stablehlo llvm.mlir.constant", Reactant.jl#3450), the gradient of the traced time
# loop aborts Julia ("LLVM ERROR: Failed to infer result type(s)") before 0.2.294 (Reactant.jl#3478, #3479).
# (The `UV_from_vor(div)!` kernels are branch-free as their branches were raised to a reduction
# without a reverse rule in Enzyme-MLIR: "Unsupported operation in reduction rev autodiff")

# The compile options of every compilation in SpeedyWeather plus `raise_first = true` for the differentiation.
# Further options (e.g. `excluded_passes` to work around a Reactant bug) can be added here, they are only
# accepted through `compile_options`, not as keywords of `@jit`.
const REACTANT_AD_COMPILE_OPTIONS = SpeedyWeather.reactant_compile_options(; raise_first = true)

"""Scalar loss to differentiate: squared norm of the prognostic vorticity after one time step."""
function loss_one_step(vars, model)
    SpeedyWeather.time_step!(vars, model.time_stepping, model)
    return sum(abs2, vars.prognostic.vorticity.data)
end

"""Reverse-mode gradient of `loss_one_step` with respect to the prognostic vorticity on the CPU."""
function gradient_vorticity_cpu(vars, model)
    dvars = Enzyme.make_zero(vars)
    Enzyme.autodiff(set_runtime_activity(Reverse), loss_one_step, Active, Duplicated(vars, dvars), Const(model))
    return Array(dvars.prognostic.vorticity.data)
end

"""Reverse-mode gradient of `loss_one_step` with respect to the prognostic vorticity with Reactant,
compiled with all kernels raised to StableHLO."""
function gradient_vorticity_reactant(vars, model)
    dvars = Enzyme.make_zero(vars)
    @jit compile_options = REACTANT_AD_COMPILE_OPTIONS Enzyme.autodiff(Reverse, loss_one_step, Active, Duplicated(vars, dvars), Const(model))
    return on_architecture(SpeedyWeather.CPU(), dvars.prognostic.vorticity.data)
end

"""Scalar loss to differentiate: squared norm of the prognostic vorticity after `n` time steps,
with Reactant in one traced loop whose reverse-mode differentiation uses the keyword `options` of
`time_steps!` (a `NamedTuple` of `checkpointing` and `mincut`, a positional argument as `Enzyme.autodiff`
takes no keywords; the defaults for the forward pass, where they don't matter)."""
function loss_n_steps(vars, model, n, options = (;))
    SpeedyWeather.time_steps!(Simulation(vars, model), n; options...)
    return sum(abs2, vars.prognostic.vorticity.data)
end

"""Reverse-mode gradient of `loss_n_steps` with respect to the prognostic vorticity on the CPU."""
function gradient_vorticity_cpu(vars, model, n)
    dvars = Enzyme.make_zero(vars)
    Enzyme.autodiff(
        set_runtime_activity(Reverse), loss_n_steps, Active,
        Duplicated(vars, dvars), Const(model), Const(n),
    )
    return Array(dvars.prognostic.vorticity.data)
end

"""Reverse-mode gradient of `loss_n_steps` with respect to the prognostic vorticity with Reactant."""
function gradient_vorticity_reactant(vars, model, n, options)
    dvars = Enzyme.make_zero(vars)
    @jit compile_options = REACTANT_AD_COMPILE_OPTIONS Enzyme.autodiff(
        Reverse, loss_n_steps, Active,
        Duplicated(vars, dvars), Const(model), Const(n), Const(options),
    )
    return on_architecture(SpeedyWeather.CPU(), dvars.prognostic.vorticity.data)
end

"""CPU and Reactant simulations of `ModelType` in the same state, ready for one time step."""
function differentiability_simulations(ModelType::Type; Transform = SpectralTransform, spinup_steps = 5)
    simulation_cpu = initialize!(create_cpu_model(ModelType; Transform))
    simulation_reactant = initialize!(create_reactant_model(ModelType; Transform))

    # spin up both (a simulation's first run differs), then the same state on both, scaled for the time step
    run!(simulation_cpu, steps = spinup_steps)
    run!(simulation_reactant, steps = spinup_steps)
    sync_variables!(simulation_cpu, simulation_reactant)
    initialize!(simulation_cpu, steps = 1)
    initialize!(simulation_reactant, steps = 1)
    sync_variables!(simulation_cpu, simulation_reactant)
    return simulation_cpu, simulation_reactant
end

"""Test loss and gradient of one time step of `ModelType` with Reactant against Enzyme on the CPU.
`compare_with_enzyme = false` skips the gradient with vanilla Enzyme on the CPU (slow to compile for the
primitive equation models) and only checks that the Reactant gradient compiles and is finite and non-zero;
the loss is always compared with the CPU. Both models use a spectral transform of type `Transform`."""
function test_differentiability(
        ModelType::Type; Transform = SpectralTransform, compare_with_enzyme = true,
        rtol_loss = 1.0e-4, rtol_gradient = 1.0e-3
    )
    return @testset "Differentiability with Reactant: $ModelType with $(nameof(Transform))" begin
        simulation_cpu, simulation_reactant = differentiability_simulations(ModelType; Transform)
        vars_cpu, model_cpu = SpeedyWeather.unpack(simulation_cpu)
        vars_reactant, model_reactant = SpeedyWeather.unpack(simulation_reactant)

        # forward: one time step with raised kernels
        loss_cpu = loss_one_step(deepcopy(vars_cpu), model_cpu)
        loss_reactant = @jit compile_options = REACTANT_AD_COMPILE_OPTIONS loss_one_step(deepcopy(vars_reactant), model_reactant)
        @test Float64(loss_reactant) ≈ loss_cpu rtol = rtol_loss

        # reverse: gradient with respect to the vorticity before the time step
        gradient_reactant = gradient_vorticity_reactant(deepcopy(vars_reactant), model_reactant)
        @test all(isfinite, gradient_reactant)
        @test any(!iszero, gradient_reactant)

        if compare_with_enzyme
            gradient_cpu = gradient_vorticity_cpu(deepcopy(vars_cpu), model_cpu)
            @test all(isfinite, gradient_cpu)
            @test any(!iszero, gradient_cpu)
            @test gradient_reactant ≈ gradient_cpu rtol = rtol_gradient
        end
    end
end

"""Test loss and gradient of `n` time steps of `ModelType` in one traced loop with Reactant, for
several checkpointing options, against Enzyme on the CPU (`compare_with_enzyme = false`: only
that the gradients compile and are finite and non-zero, see `test_differentiability`). Reactant unrolls
traced loops of at most 4 iterations, `n > 4` keeps the `stablehlo.while` that the checkpointing options
act on. `options` are the keyword arguments of `time_steps!` to test (as `NamedTuple`s). Both models use a spectral
transform of type `Transform`."""
function test_differentiability_time_steps(
        ModelType::Type;
        Transform = SpectralTransform,
        n = 6,
        compare_with_enzyme = true,
        options = ((;), (; checkpointing = Reactant.Periodic(2)), (; checkpointing = Reactant.Binomial(2))),
        rtol_loss = 1.0e-4,
        rtol_gradient = 1.0e-3,
    )
    return @testset "Differentiability of $n time steps with Reactant: $ModelType with $(nameof(Transform))" begin
        simulation_cpu, simulation_reactant = differentiability_simulations(ModelType; Transform)
        vars_cpu, model_cpu = SpeedyWeather.unpack(simulation_cpu)
        vars_reactant, model_reactant = SpeedyWeather.unpack(simulation_reactant)

        # forward: n time steps, one program
        loss_cpu = loss_n_steps(deepcopy(vars_cpu), model_cpu, n)
        loss_reactant = @jit compile_options = REACTANT_AD_COMPILE_OPTIONS loss_n_steps(deepcopy(vars_reactant), model_reactant, n)
        @test Float64(loss_reactant) ≈ loss_cpu rtol = rtol_loss

        # reverse: gradient with respect to the vorticity before the n time steps
        if compare_with_enzyme
            gradient_cpu = gradient_vorticity_cpu(deepcopy(vars_cpu), model_cpu, n)
            @test all(isfinite, gradient_cpu)
            @test any(!iszero, gradient_cpu)
        end

        # `(;)` stores every step; mincut = false isn't tested: its compilation of the loop gradient doesn't
        # finish (> 30 min)
        @testset "$option" for option in options
            gradient_reactant = gradient_vorticity_reactant(deepcopy(vars_reactant), model_reactant, n, option)
            @test all(isfinite, gradient_reactant)
            @test any(!iszero, gradient_reactant)
            compare_with_enzyme && @test gradient_reactant ≈ gradient_cpu rtol = rtol_gradient
        end
    end
end

# The tests run in CI, skipped when another script includes this file for its functions only
if !(@isdefined(DIFFERENTIATION_FUNCTIONS_ONLY) && DIFFERENTIATION_FUNCTIONS_ONLY)
    @testset "make_zero(::Variables) with Reactant" begin
        # variables are views into fused parents, make_zero relies on deepcopy preserving that aliasing
        simulation = initialize!(create_reactant_model(BarotropicModel))
        run!(simulation, steps = 1)     # to have non-zero grid variables too
        vars, _ = SpeedyWeather.unpack(simulation)
        dvars = Enzyme.make_zero(vars)
        for (primal, shadow) in ((vars.prognostic.vorticity, dvars.prognostic.vorticity), (vars.grid.u, dvars.grid.u))
            @test any(!iszero, Array(primal.data))
            @test all(iszero, Array(shadow.data))
        end
    end

    test_differentiability(BarotropicModel)                                       # SpectralTransform (truncated DFT)
    test_differentiability(BarotropicModel; Transform = MatrixSpectralTransform)
    test_differentiability_time_steps(BarotropicModel)

    # PrimitiveWetModel: only that the gradient of the traced loop compiles (and is finite), the comparison
    # with Enzyme on the CPU takes too long to compile for CI (`differentiation_primitive_wet.jl`). Without
    # `Reactant.Periodic(2)`: with Reactant 0.2.294 its compilation doesn't finish within an hour (on 0.2.293
    # with the then needed pattern exclusions it did), the other options compile in ~3.5 min each
    test_differentiability_time_steps(
        PrimitiveWetModel; compare_with_enzyme = false, options = ((;), (; checkpointing = Reactant.Binomial(2))),
    )
end
