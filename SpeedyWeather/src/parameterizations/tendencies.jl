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
    _column_parameterizations_cpu!(vars, get_parameterizations(model), model)
    return nothing
end

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

"""Abstract supertype for the time integration of surface fluxes that relax the
lowermost layer towards a surface value, `flux = -conductance * (x - x_surface)`.
Subtypes define `flux_damping(scheme, z)` with `z = Δt/τ` the ratio of time step to
relaxation time scale, a factor that multiplies the explicit flux."""
abstract type AbstractFluxIntegration end
Base.broadcastable(scheme::AbstractFluxIntegration) = Ref(scheme)

export ExplicitFluxes, ImplicitFluxes, ExponentialFluxes, RationalFluxes

"""Explicit (forward) time integration of surface fluxes, unstable for `Δt/τ > 2`."""
struct ExplicitFluxes <: AbstractFluxIntegration end

"""Implicit (backward Euler) time integration of surface fluxes, damping the explicit
flux by `1/(1 + Δt/τ)`. Unconditionally stable."""
struct ImplicitFluxes <: AbstractFluxIntegration end

"""Exponential time differencing of surface fluxes, damping the explicit flux by
`(1 - exp(-Δt/τ))/(Δt/τ)`. Exact for a fixed surface value, unconditionally stable."""
struct ExponentialFluxes <: AbstractFluxIntegration end

"""Rational approximation of `ExponentialFluxes`, damping the explicit flux by
`(6 + z)/(6 + 4z + z²)` with `z = Δt/τ`. Matches `(1 - exp(-z))/z` to second order
in `z` and as `1/z` for large `z`, at the cost of `ImplicitFluxes`. Unconditionally
stable but overshoots the surface value slightly (≤ 10%) for `z > 3`."""
struct RationalFluxes <: AbstractFluxIntegration end

const DEFAULT_FLUX_INTEGRATION = RationalFluxes

"""$(TYPEDSIGNATURES)
Factor φ(z) multiplying an explicit relaxation flux for `z = Δt/τ` given time integration `scheme`."""
@inline flux_damping(::ExplicitFluxes, z) = one(z)
@inline flux_damping(::ImplicitFluxes, z) = 1 / (1 + z)
@inline flux_damping(::ExponentialFluxes, z) = ifelse(z < sqrt(eps(z)), one(z), -expm1(-z) / z)
@inline flux_damping(::RationalFluxes, z) = (6 + z) / (6 + 4z + z^2)

"""$(TYPEDSIGNATURES)
Damping factor for a surface flux into the lowermost layer with `conductance = -∂flux/∂x`
[kg/m²/s] with respect to the lowermost-layer variable x, surface pressure `pₛ` [Pa] and
time integration `scheme`. Uses the time step of the time stepper, e.g. 2Δt for Leapfrog."""
@propagate_inbounds function surface_flux_damping(scheme::AbstractFluxIntegration, conductance::Real, pₛ::Real, model)
    Δp = pressure_thickness(model.geometry.nlayers, pₛ, model.geometry.vertical_coordinates)
    Δt = default_time_step(model.time_stepping)
    z = Δt * model.planet.gravity * conductance / Δp        # = Δt/τ, relaxation time scale τ
    return flux_damping(scheme, z)
end

"""$(TYPEDSIGNATURES)
Flux `flux` (upward positive) into the lowermost layer with `conductance = -∂flux/∂x`
converted to tendency [?/s], time-integrated following `scheme`, see `surface_flux_damping`."""
@propagate_inbounds surface_flux_to_tendency(scheme::AbstractFluxIntegration, flux::Real, conductance::Real, pₛ::Real, model) =
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
    reset_variable!(vars.parameterizations, :cloud_top, nlayers + 1)   # reset to below top layer
    reset_variable!(vars.parameterizations, :rain_rate, 0)
    reset_variable!(vars.parameterizations, :snow_rate, 0)
    reset_variable!(vars.parameterizations, :rain_rate_convection, 0)
    reset_variable!(vars.parameterizations, :snow_rate_convection, 0)
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
