"""$(TYPEDSIGNATURES)
Compute tendencies for u, v, temp, humid from physical parametrizations."""
function parameterization_tendencies!(
        vars::Variables,
        model::PrimitiveEquation,
    )
    reset_variables!(vars)

    # parameterizations with their own kernel
    global_parameterizations!(vars, model)

    # parameterizations fused into a single kernel over ij
    column_parameterizations!(vars, model)

    return nothing
end

function global_parameterizations!(vars::Variables, model::PrimitiveEquation)
    _global_parameterizations!(vars, get_parameterizations(model), model)
    return nothing
end

# Use @generated to unroll NamedTuple iteration at compile time
@generated function _global_parameterizations!(vars, parameterizations::NamedTuple{names}, model) where {names}
    calls = [:(parameterization!(vars, parameterizations.$name, model)) for name in names]
    return quote
        $(Expr(:block, calls...))
    end
end

# COLUMN-BASED PARAMETERIZATIONS
function column_parameterizations!(vars, model)
    (; architecture, npoints) = model.spectral_grid
    if architecture isa Architectures.AbstractCPU
        # bypass kernel launch on CPU
        column_parameterizations_cpu!(vars, model)
    else
        # GPU: all other parameterizations are fused into a single kernel over horizontal grid point index ij
        launch!(
            architecture, LinearWorkOrder, (npoints,), column_parameterizations_kernel!,
            vars, get_parameterizations(model), model
        )
    end
    return nothing
end

# GPU kernel, unrolling NamedTuple iteration at compile time, fuses all parameterizations
@kernel inbounds = true function column_parameterizations_kernel!(vars, parameterizations, model)

    ij = @index(Global, Linear)     # every horizontal grid point ij

    # manually unroll loop over all parameterizations (NamedTuple iteration not GPU-compatible)
    column_parameterizations!(ij, vars, parameterizations, model)
end

# Use @generated to unroll NamedTuple iteration at compile time for GPU compatibility
@generated function column_parameterizations!(ij, vars, parameterizations::NamedTuple{names}, model) where {names}
    calls = [:(parameterization!(ij, vars, parameterizations.$name, model)) for name in names]
    return quote
        Base.@_propagate_inbounds_meta
        $(Expr(:block, calls...))
    end
end

# CPU without kernel, just a loop, change loop order compared to GPU though:
# outer loop over parameterizations, inner loop over horizontal grid points
# this yields a more contiguous memory access pattern on CPU
function column_parameterizations_cpu!(vars, model)
    if SPLIT_COLUMN_LOOPS       # workaround for llvm/llvm-project#219885, see below
        _column_parameterizations_cpu_by_name!(vars, model, Val(keys(get_parameterizations(model))))
    else
        _column_parameterizations_cpu!(vars, get_parameterizations(model), model)
    end
    return nothing
end

# >>> WORKAROUND llvm/llvm-project#219885 (LLVM GVN bug with `!invariant.group` loads) >>>>>>>>>>>>>>>>>>>>>>
# Inlined, all column loops end up in one reverse function of Enzyme large enough to trigger an LLVM bug
# in Enzyme's post-AD GVN (stale MemDep cache for the `!invariant.group` loads that Enzyme emits): a broken
# module ("Broken module found") or a segfault while compiling the reverse pass, see #1290 and
# docs/dev/2026-11/2026-11-17-extended-diff-tests-failure.md. With SPLIT_COLUMN_LOOPS every parameterization
# runs in its own non-inlined function, also in Enzyme's reverse pass. No measurable cost (the call is per
# parameterization and time step, not per grid point). Fixed upstream by llvm/llvm-project#219929, which
# is not in the LLVM of Julia 1.10, 1.11 or 1.12 (2026-11).
# TO REVERT once the Julia versions used for Enzyme CI ship an LLVM with that fix: a quick check is
# `opt -passes=gvn` on docs/dev/2026-11/2026-11-17-gvn-invariant-group{,-upstream}.ll with the `opt` of the
# LLVM_full_jll matching Julia's LLVM (see the headers of these files, valid output = fixed). Decisive: set
# SPLIT_COLUMN_LOOPS = false and run test/differentiability/primitivewet.jl (the "parameter AD" testset aborted).
# Then delete this block, the branch above and test/parameterizations/column_loops.jl.
# xref Enzyme.jl issue: EnzymeAD/Enzyme.jl#4074
const SPLIT_COLUMN_LOOPS = true

# Unrolled at compile time over the names of the model's parameterizations, one call per parameterization
@generated function _column_parameterizations_cpu_by_name!(vars, model, ::Val{names}) where {names}
    calls = [:(_column_parameterization_cpu_by_name!(vars, Val($(QuoteNode(name))), model)) for name in names]
    return Expr(:block, calls..., :(return nothing))
end

# The parameterization is selected by name inside: passing the struct to a non-inlined function makes
# Enzyme copy it untyped (EnzymeNoTypeError with a Duplicated model)
@noinline function _column_parameterization_cpu_by_name!(vars, ::Val{name}, model) where {name}
    column_parameterization_cpu!(vars, getfield(model, name), model)
    return nothing
end
# <<< END WORKAROUND llvm/llvm-project#219885 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Use @generated to unroll NamedTuple iteration at compile time also on CPU for performance
@generated function _column_parameterizations_cpu!(vars, parameterizations::NamedTuple{names}, model) where {names}
    # parameterizations outer loop, grid points inner loop inside column_parameterization_cpu!
    calls = [:(column_parameterization_cpu!(vars, parameterizations.$name, model)) for name in names]
    return quote
        Base.@_propagate_inbounds_meta
        $(Expr(:block, calls...))
    end
end

"""$(TYPEDSIGNATURES)
Loop over all horizontal grid points `ij` for a single column `parameterization` on CPU.
Parameterizations that bundle several schemes (e.g. `Radiation`) extend this method to loop
over the grid points for every scheme separately, keeping the memory access contiguous per
scheme and the loop bodies small."""
@propagate_inbounds function column_parameterization_cpu!(vars, parameterization, model)
    for ij in 1:model.geometry.npoints      # horizontal grid points inner loop
        parameterization!(ij, vars, parameterization, model)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Flux `flux` into surface layer with surface pressure `pₛ` [Pa] and gravity `g` [m/s^2]
converted to tendency [?/s]."""
@propagate_inbounds surface_flux_to_tendency(flux::Real, pₛ::Real, model) =
    flux_to_tendency(flux, pₛ, model.geometry.nlayers, model)

"""$(TYPEDSIGNATURES)
Damping factor for a surface flux into the lowermost layer with mass `conductance` [kg/m²/s],
i.e. `-∂flux/∂x` with respect to the lowermost-layer variable x for the flux in units of
x⋅kg/m²/s (so without cₚ for heat fluxes), surface pressure `pₛ` [Pa] and
time integration `scheme`. Uses the time step of the time stepper, e.g. 2Δt for Leapfrog."""
@propagate_inbounds function surface_flux_damping(scheme::AbstractRelaxationIntegration, conductance::Real, pₛ::Real, model)
    Δp = pressure_thickness(model.geometry.nlayers, pₛ, model.geometry.vertical_coordinates)
    Δt = default_time_step(model.time_stepping)
    z = Δt * model.planet.gravity * conductance / Δp        # = Δt/τ, relaxation time scale τ
    return relaxation(scheme, z)
end

"""$(TYPEDSIGNATURES)
Flux `flux` (upward positive) into the lowermost layer with mass `conductance` [kg/m²/s]
converted to tendency [?/s], time-integrated following `scheme`, see `surface_flux_damping`."""
@propagate_inbounds surface_flux_to_tendency(scheme::AbstractRelaxationIntegration, flux::Real, conductance::Real, pₛ::Real, model) =
    surface_flux_damping(scheme, conductance, pₛ, model) * surface_flux_to_tendency(flux, pₛ, model)

"""$(TYPEDSIGNATURES)
Flux `flux` into layer `k` of thickness `Δσ`  converted to tendency [?/s].
Using surface pressure `pₛ` [Pa] and gravity `g` [m/s^2]."""
@propagate_inbounds flux_to_tendency(flux::Real, pₛ::Real, g::Real, Δσ_k::Real) = g / (pₛ * Δσ_k) * flux
@propagate_inbounds function flux_to_tendency(flux::Real, pₛ::Real, k::Integer, model)
    Δp = pressure_thickness(k, pₛ, model.geometry.vertical_coordinates)
    g = model.planet.gravity
    return g * flux / Δp
end

# hacky, temporary placement, and also modularize this?
function reset_variables!(vars::Variables)
    nlayers = size(vars.prognostic.temperature, 2)
    resets = (
        cloud_top = nlayers + 1,                # reset to below top layer
        rain_rate = 0,
        snow_rate = 0,
        rain_rate_convection = 0,
        snow_rate_convection = 0,
        surface_humidity_flux = 0,
        sensible_heat_flux = 0,
    )
    return reset_variables!(vars.parameterizations, resets)
end

"""$(TYPEDSIGNATURES)
Reset every variable in `vars` that is named in `resets` to its value there, `vars.name .= resets.name`,
skipping names that `vars` doesn't have. Generated from the names of both NamedTuples, so the variables
are picked at compile time, one concretely inferred broadcast each, however many variables `vars` holds."""
@generated function reset_variables!(vars::NamedTuple{names}, resets::NamedTuple{reset_names}) where {names, reset_names}
    exprs = [
        :(getfield(vars, $(QuoteNode(name))) .= getfield(resets, $(QuoteNode(name))))
            for name in reset_names if name in names
    ]
    return Expr(:block, exprs..., :(return nothing))
end
