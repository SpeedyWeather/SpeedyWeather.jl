# Differentiability with Reactant: Enzyme reverse mode through one time step. Enzyme-MLIR can only
# differentiate KernelAbstractions kernels that are raised to StableHLO, hence `raise = true`.
# The Reactant gradient is compared against Enzyme on the CPU for the same state.
# CI runs this for the BarotropicModel, `differentiation_primitive_wet.jl` uses the functions
# defined here for the PrimitiveWetModel (too slow for CI).
#
# Needs Reactant_jll >= v0.0.419: raising kernels that store a complex constant (e.g. `div[I] = 0`)
# fails before ("cannot raise op to stablehlo llvm.mlir.constant", Reactant.jl#3450).
# (The `UV_from_vor(div)!` kernels are branch-free as their branches were raised to a reduction
# without a reverse rule in Enzyme-MLIR: "Unsupported operation in reduction rev autodiff")

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
    @jit raise = true raise_first = true Enzyme.autodiff(Reverse, loss_one_step, Active, Duplicated(vars, dvars), Const(model))
    return on_architecture(SpeedyWeather.CPU(), dvars.prognostic.vorticity.data)
end

"""CPU and Reactant simulations of `ModelType` in the same state, ready for one time step."""
function differentiability_simulations(ModelType::Type; spinup_steps = 5)
    simulation_cpu = initialize!(create_cpu_model(ModelType))
    simulation_reactant = initialize!(create_reactant_model(ModelType))

    # spin up both (a simulation's first run differs), then the same state on both, scaled for the time step
    run!(simulation_cpu, steps = spinup_steps)
    run!(simulation_reactant, steps = spinup_steps)
    sync_variables!(simulation_cpu, simulation_reactant)
    initialize!(simulation_cpu, steps = 1)
    initialize!(simulation_reactant, steps = 1)
    sync_variables!(simulation_cpu, simulation_reactant)
    return simulation_cpu, simulation_reactant
end

"""Test loss and gradient of one time step of `ModelType` with Reactant against Enzyme on the CPU."""
function test_differentiability(ModelType::Type; rtol_loss = 1.0e-4, rtol_gradient = 1.0e-3)
    return @testset "Differentiability with Reactant: $ModelType" begin
        simulation_cpu, simulation_reactant = differentiability_simulations(ModelType)
        vars_cpu, model_cpu = SpeedyWeather.unpack(simulation_cpu)
        vars_reactant, model_reactant = SpeedyWeather.unpack(simulation_reactant)

        # forward: one time step with raised kernels
        loss_cpu = loss_one_step(deepcopy(vars_cpu), model_cpu)
        loss_reactant = @jit raise = true raise_first = true loss_one_step(deepcopy(vars_reactant), model_reactant)
        @test Float64(loss_reactant) ≈ loss_cpu rtol = rtol_loss

        # reverse: gradient with respect to the vorticity before the time step
        gradient_cpu = gradient_vorticity_cpu(deepcopy(vars_cpu), model_cpu)
        @test all(isfinite, gradient_cpu)
        @test any(!iszero, gradient_cpu)

        gradient_reactant = gradient_vorticity_reactant(deepcopy(vars_reactant), model_reactant)
        @test gradient_reactant ≈ gradient_cpu rtol = rtol_gradient
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

    test_differentiability(BarotropicModel)
end
