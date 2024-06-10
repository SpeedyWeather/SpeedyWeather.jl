"""$(TYPEDSIGNATURES)
Calculates the virtual temperature Tᵥ as

    Tᵥ = T(1+μq)

With absolute temperature T, specific humidity q and

    μ = (1-ξ)/ξ, ξ = R_dry/R_vapour.
    
in grid-point space."""
function virtual_temperature!(
    diagn::DiagnosticVariables,
    model::PrimitiveWet,
    )
    
    (; temp_grid, humid_grid, temp_virt_grid) = diagn.grid
    μ = model.atmosphere.μ_virt_temp
    # @. temp_virt_grid = temp_grid * (1 + μ*humid_grid)    # same as loop but a bit slower (?)

    @inbounds for ijk in eachindex(temp_virt_grid, temp_grid, humid_grid)
        temp_virt_grid[ijk] = temp_grid[ijk]*(1 + μ*humid_grid[ijk])
    end
end

"""
$(TYPEDSIGNATURES)
Virtual temperature in grid-point space: For the PrimitiveDry temperature
and virtual temperature are the same (humidity=0). Just copy over the arrays."""
function virtual_temperature!(
    diagn::DiagnosticVariables,
    model::PrimitiveDry,
)   
    (; temp_grid, temp_virt_grid) = diagn.grid
    # Tᵥ = T(1 + μ*q) with humid=q=0
    temp_virt_grid .= temp_grid
    return nothing
end

"""
$(TYPEDSIGNATURES)
Linear virtual temperature for `model::PrimitiveDry`: Just copy over
arrays from `temp` to `temp_virt` at timestep `lf` in spectral space
as humidity is zero in this `model`."""
function linear_virtual_temperature!(   
    diagn::DiagnosticVariablesLayer,
    progn::PrognosticLayerTimesteps,
    model::PrimitiveDry,
    lf::Integer,
)
    (; temp_virt) = diagn.dynamics_variables
    (; temp) = progn.timesteps[lf]
    copyto!(temp_virt, temp)
end

# function barrier
function linear_virtual_temperature!(  
    diagn::DiagnosticVariablesLayer,
    progn::PrognosticLayerTimesteps,
    model::PrimitiveWet,
    lf::Integer,
)
    linear_virtual_temperature!(diagn, progn, model.atmosphere, lf)
end

"""
$(TYPEDSIGNATURES)
Calculates a linearised virtual temperature Tᵥ as

    Tᵥ = T + Tₖμq

With absolute temperature T, layer-average temperarture Tₖ (computed in temperature_average!),
specific humidity q and

    μ = (1-ξ)/ξ, ξ = R_dry/R_vapour.
    
in spectral space."""
function linear_virtual_temperature!(   diagn::DiagnosticVariablesLayer,
                                        progn::PrognosticLayerTimesteps,
                                        atmosphere::AbstractAtmosphere,
                                        lf::Int)
    
    (; temp_virt) = diagn.dynamics_variables
    μ = atmosphere.μ_virt_temp
    Tₖ = diagn.temp_average[]   
    (; temp, humid) = progn.timesteps[lf]

    # TODO check that doing a non-linear virtual temperature in grid-point space
    # but a linear virtual temperature in spectral space to avoid another transform
    # does not cause any problems. Alternative do the transform or have a linear
    # virtual temperature in both grid and spectral space
    # spectral!(temp_virt, temp_virt_grid, S)
    @. temp_virt = temp + (Tₖ*μ)*humid
end