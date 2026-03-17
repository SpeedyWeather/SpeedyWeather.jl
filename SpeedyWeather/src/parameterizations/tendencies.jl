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

    # tendencies have to be scaled by the radius for the dynamical core
    scale_tendencies!(ij, vars.tendencies.grid, model.planet.radius)
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
    _column_parameterizations_cpu!(vars, get_parameterizations(model), model)

    # tendencies have to be scaled by the radius for the dynamical core
    @inbounds scale_tendencies!(vars.tendencies.grid, model.planet.radius)
    return nothing
end

# Use @generated to unroll NamedTuple iteration at compile time also on CPU for performance
@generated function _column_parameterizations_cpu!(vars, parameterizations::NamedTuple{names}, model) where {names}
    # runic: off
    calls = [
        quote
            for ij in 1:model.geometry.npoints      # horizontal grid points inner loop
                parameterization!(ij, vars, parameterizations.$name, model)
            end
        end for name in names                       # parameterizations outer loop
    ]
    # runic: on
    return quote
        Base.@_propagate_inbounds_meta
        $(Expr(:block, calls...))
    end
end

"""$(TYPEDSIGNATURES)
Flux `flux` into surface layer with surface pressure `pₛ` [Pa] and gravity `g` [m/s^2]
converted to tendency [?/s]."""
@propagate_inbounds surface_flux_to_tendency(flux::Real, pₛ::Real, model) =
    flux_to_tendency(flux, pₛ, model.planet.gravity, model.geometry.σ_levels_thick[end])

"""$(TYPEDSIGNATURES)
Flux `flux` into layer `k` of thickness `Δσ`  converted to tendency [?/s].
Using surface pressure `pₛ` [Pa] and gravity `g` [m/s^2]."""
@propagate_inbounds flux_to_tendency(flux::Real, pₛ::Real, g::Real, Δσ_k::Real) = g / (pₛ * Δσ_k) * flux
@propagate_inbounds flux_to_tendency(flux::Real, pₛ::Real, k::Int, model) =
    flux_to_tendency(flux, pₛ, model.planet.gravity, model.geometry.σ_levels_thick[k])

# hacky, temporary placement, and also modularize this?
function reset_variables!(vars::Variables)
    nlayers = size(vars.prognostic.temp, 2)
    reset_variable!(vars.parameterizations, :cloud_top, nlayers + 1)   # reset to below top layer
    reset_variable!(vars.parameterizations, :rain_rate, 0)
    reset_variable!(vars.parameterizations, :snow_rate, 0)
    reset_variable!(vars.parameterizations, :surface_humidity_flux, 0)
    reset_variable!(vars.parameterizations, :sensible_heat_flux, 0)
    return nothing
end

function reset_variable!(vars, var::Symbol, reset_value)
    if haskey(vars, var)
        field = getfield(vars, var)
        field .= reset_value
    end
    return nothing
end
