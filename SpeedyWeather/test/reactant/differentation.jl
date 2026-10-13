# Differentiability with Reactant: Enzyme reverse mode through one time step. Enzyme-MLIR can only
# differentiate KernelAbstractions kernels that are raised to StableHLO, hence `raise = true`.
# The Reactant gradient is compared against Enzyme on the CPU for the same state.
#
# Currently broken upstream (Reactant v0.2.292, Reactant_jll v0.0.418): raising fails for kernels
# that store a complex constant, e.g. `div[I] = 0` in `_divergence_kernel!`: "cannot raise op to
# stablehlo llvm.mlir.constant". Fixed on Reactant's main, to be released. With that kernel storing
# a computed zero instead the gradient agrees with the CPU to a relative error of 1e-6.
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

"""Evaluate `f()` but return `nothing` and log the first line of the error if it throws,
to mark the test that depends on it as broken instead of erroring the test set."""
function try_or_nothing(f, description)
    try
        return f()
    catch err
        @info "$description failed: $(first(split(sprint(showerror, err), '\n')))"
        return nothing
    end
end

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

@testset "Differentiability with Reactant: BarotropicModel" begin
    simulation_cpu = initialize!(create_cpu_model(BarotropicModel))
    simulation_reactant = initialize!(create_reactant_model(BarotropicModel))

    # spin up both (a simulation's first run differs), then the same state on both, scaled for the time step
    run!(simulation_cpu, steps = 5)
    run!(simulation_reactant, steps = 5)
    sync_variables!(simulation_cpu, simulation_reactant)
    initialize!(simulation_cpu, steps = 1)
    initialize!(simulation_reactant, steps = 1)
    sync_variables!(simulation_cpu, simulation_reactant)

    vars_cpu, model_cpu = SpeedyWeather.unpack(simulation_cpu)
    vars_reactant, model_reactant = SpeedyWeather.unpack(simulation_reactant)

    # forward: one time step with raised kernels
    loss_cpu = loss_one_step(deepcopy(vars_cpu), model_cpu)
    loss_reactant = try_or_nothing("Raising the time step") do
        @jit raise = true raise_first = true loss_one_step(deepcopy(vars_reactant), model_reactant)
    end
    @test_broken !isnothing(loss_reactant) && isapprox(Float64(loss_reactant), loss_cpu; rtol = 1.0e-4)

    # reverse: gradient with respect to the vorticity before the time step
    gradient_cpu = gradient_vorticity_cpu(deepcopy(vars_cpu), model_cpu)
    @test all(isfinite, gradient_cpu)
    @test any(!iszero, gradient_cpu)

    gradient_reactant = try_or_nothing("Reverse-mode differentiation of the time step") do
        gradient_vorticity_reactant(deepcopy(vars_reactant), model_reactant)
    end
    @test_broken !isnothing(gradient_reactant) && isapprox(gradient_reactant, gradient_cpu; rtol = 1.0e-3)
end
