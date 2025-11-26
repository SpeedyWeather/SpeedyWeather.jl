"""$(TYPEDSIGNATURES)
Compute tendencies for u, v, temp, humid from physical parametrizations."""
function parameterization_tendencies!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    model::PrimitiveEquation,
)
    # parameterizations with their own kernel
    (; time) = progn.clock
    cos_zenith!(diagn, time, model)

    model_parameters = get_model_parameters(model)      # subset of GPU-compatible model components
    parameterizations = get_parameterizations(model)    # subset of model: parameterizations only

    # all other parameterizations are fused into a single kernel over horizontal grid point index ij
    (; architecture, npoints) = model.spectral_grid
    if architecture isa Architectures.AbstractCPU
        parameterization_tendencies_cpu!(diagn, progn, parameterizations, model_parameters)
    else
        launch!(architecture, LinearWorkOrder, (npoints,), parameterization_tendencies_kernel!,
            diagn, progn, parameterizations, model_parameters)
    end
    return nothing
end

# GPU kernel
@kernel function parameterization_tendencies_kernel!(diagn, progn, @Const(parameterizations), @Const(model))
    
    ij = @index(Global, Linear)     # every horizontal grid point ij

    # manually unroll loop over all parameterizations (NamedTuple iteration not GPU-compatible)
    _call_parameterizations!(ij, diagn, progn, parameterizations, model)

    # tendencies have to be scaled by the radius for the dynamical core
    scale!(ij, diagn.tendencies, model.planet.radius)
end

# CPU without kernel, just a loop
function parameterization_tendencies_cpu!(diagn, progn, parameterizations, model)
    _call_parameterizations_cpu!(diagn, progn, parameterizations, model)

    radius = model.planet.radius
    for ij in 1:model.geometry.npoints
        # tendencies have to be scaled by the radius for the dynamical core
        scale!(ij, diagn.tendencies, radius)
    end
end

# Use @generated to unroll NamedTuple iteration at compile time for GPU compatibility
@generated function _call_parameterizations!(ij, diagn, progn, parameterizations::NamedTuple{names}, model) where {names}
    calls = [:(parameterization!(ij, diagn, progn, parameterizations.$name, model)) for name in names]
    return Expr(:block, calls...)
end

# Use @generated to unroll NamedTuple iteration at compile time also on CPU for performance
@generated function _call_parameterizations_cpu!(diagn, progn, parameterizations::NamedTuple{names}, model) where {names}
    calls = [:(
        for ij in 1:model.geometry.npoints
            parameterization!(ij, diagn, progn, parameterizations.$name, model)
        end
        ) for name in names]
    return Expr(:block, calls...)
end

"""$(TYPEDSIGNATURES)
Flux `flux` into surface layer with surface pressure `pₛ` [Pa] and gravity `g` [m/s^2]
converted to tendency [?/s]."""
@inline surface_flux_to_tendency(flux::Real, pₛ::Real, model) =
    flux_to_tendency(flux, pₛ, model.planet.gravity, model.geometry.σ_levels_thick[end])

"""$(TYPEDSIGNATURES)
Flux `flux` into layer `k` of thickness `Δσ`  converted to tendency [?/s].
Using surface pressure `pₛ` [Pa] and gravity `g` [m/s^2]."""
@inline flux_to_tendency(flux::Real, pₛ::Real, g::Real, Δσ_k::Real) = g/(pₛ*Δσ_k) * flux