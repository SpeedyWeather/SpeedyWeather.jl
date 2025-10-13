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

    model_parameters = get_model_parameters(model)      # subsection of GPU-compatible model components
    parameterizations = get_parameterizations(model)    # subsection of model: parameterizations only

    # all other parameterizations are fused into a single kernel over horizontal grid point index ij
    (; architecture, npoints) = model.spectral_grid
    launch!(architecture, LinearWorkOrder, (npoints,), parameterization_tendencies_kernel!,
        diagn, progn, parameterizations, model_parameters)
    return nothing
end

@kernel function parameterization_tendencies_kernel!(diagn, progn, parameterizations, model_parameters)
    
    ij = @index(Global, Linear)     # every horizontal grid point ij

    # loop over all parameterizations in order
    for parameterization in parameterizations
        parameterization!(ij, diagn, progn, parameterization, model_parameters)
    end

    # tendencies have to be scaled by the radius for the dynamical core
    scale!(ij, diagn, model_parameters.planet.radius)
end

# function barrier
flux_to_tendency(k, flux, pₛ, model) =
    flux_to_tendency(k, flux, pₛ, model.planet.gravity, model.geometry.σ_levels_thick)

"""$(TYPEDSIGNATURES)
Flux `flux` into layer `k` converted to tendency [?/s]"""
flux_to_tendency(k, flux, pₛ, g, Δσ) = g/(pₛ*Δσ[k]) * flux