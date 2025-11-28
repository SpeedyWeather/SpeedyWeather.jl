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

    #model_parameters = get_model_parameters(model)      # subset of GPU-compatible model components
    #parameterizations = get_parameterizations(model)    # subset of model: parameterizations only

    (; architecture, npoints) = model.spectral_grid
    if architecture isa Architectures.AbstractCPU
        # bypass kernel launch on CPU
        parameterization_tendencies_cpu!(diagn, progn, model)
    else
        # GPU: all other parameterizations are fused into a single kernel over horizontal grid point index ij
        launch!(architecture, LinearWorkOrder, (npoints,), parameterization_tendencies_kernel!,
            diagn, progn, model)
    end
    return nothing
end

# GPU kernel, unrolling NamedTuple iteration at compile time, fuses all parameterizations
@kernel function parameterization_tendencies_kernel!(diagn, progn, model)
    
    ij = @index(Global, Linear)     # every horizontal grid point ij

    # manually unroll loop over all parameterizations (NamedTuple iteration not GPU-compatible)
    _call_parameterizations!(ij, diagn, progn, model)

    # tendencies have to be scaled by the radius for the dynamical core
    scale!(ij, diagn.tendencies, model.planet.radius)
end

# CPU without kernel, just a loop, change loop order compared to GPU though:
# outer loop over parameterizations, inner loop over horizontal grid points
# this yields a more contiguous memory access pattern on CPU
function parameterization_tendencies_cpu!(diagn, progn, model)
    _call_parameterizations_cpu!(diagn, progn, model)

    radius = model.planet.radius
    for ij in 1:model.geometry.npoints
        # tendencies have to be scaled by the radius for the dynamical core
        scale!(ij, diagn.tendencies, radius)
    end
end

# GPU kernel version - unroll at compile time using @generated
@generated function _call_parameterizations!(ij, diagn, progn, model::ModelType) where ModelType <: PrimitiveEquation
    # Get the type of the parameterizations field - it's an NTuple with length in the type
    parameterizations_type = fieldtype(ModelType, :parameterizations)
    N = length(parameterizations_type.parameters)
    
    # Generate a loop with compile-time known bounds
    return quote
        Base.@nexprs $N i -> parameterization!(ij, diagn, progn, model.parameterizations[i], model)
        return nothing
    end
end

# CPU version - unroll outer loop over parameterizations at compile time
@generated function _call_parameterizations_cpu!(diagn, progn, model::ModelType) where ModelType <: PrimitiveEquation
    # Extract the Val type parameter from the params field
    params_type = fieldtype(ModelType, :params)
    params_tuple = params_type.parameters[1]
    
    # Generate nested loops: outer over parameterizations, inner over grid points
    calls = Expr[]
    for param_name in params_tuple
        push!(calls, quote
            for ij in 1:model.geometry.npoints
                parameterization!(ij, diagn, progn, model.$param_name, model)
            end
        end)
    end
    
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