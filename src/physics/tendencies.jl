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
    launch!(architecture, LinearWorkOrder, (npoints,), parameterization_tendencies_kernel!,
        diagn, progn, parameterizations, model_parameters)
    return nothing
end

@kernel function parameterization_tendencies_kernel!(diagn, progn, @Const(parameterizations), @Const(model_parameters))
    
    ij = @index(Global, Linear)     # every horizontal grid point ij

    # loop over all parameterizations in order
    for parameterization in parameterizations
        parameterization!(ij, diagn, progn, parameterization, model_parameters)
    end

    # tendencies have to be scaled by the radius for the dynamical core
    scale!(ij, diagn.tendencies, model_parameters.planet.radius)
end

"""$(TYPEDSIGNATURES)
Flux `flux` into surface layer with surface pressure `pₛ` [Pa] and gravity `g` [m/s^2]
converted to tendency [?/s]."""
surface_flux_to_tendency(flux::Real, pₛ::Real, model) =
    flux_to_tendency(flux, pₛ, model.planet.gravity, model.geometry.σ_levels_thick[end])

"""$(TYPEDSIGNATURES)
Flux `flux` into layer `k` of thickness `Δσ`  converted to tendency [?/s].
Using surface pressure `pₛ` [Pa] and gravity `g` [m/s^2]."""
flux_to_tendency(flux::Real, pₛ::Real, g::Real, Δσ_k::Real) = g/(pₛ*Δσ_k) * flux