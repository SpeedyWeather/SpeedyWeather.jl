# function barrier
function virtual_temperature!(  diagn::DiagnosticVariablesLayer,
                                temp::LowerTriangularMatrix,    # only needed for dispatch compat with DryCore
                                model::PrimitiveWet)
    virtual_temperature!(diagn, temp, model.atmosphere)
end

"""
$(TYPEDSIGNATURES)
Calculates the virtual temperature Tᵥ as

    Tᵥ = T(1+μq)

With absolute temperature T, specific humidity q and

    μ = (1-ξ)/ξ, ξ = R_dry/R_vapour.
    
in grid-point space."""
function virtual_temperature!(
    diagn::DiagnosticVariablesLayer,
    temp::LowerTriangularMatrix,    # only needed for dispatch compat with DryCore
    atmosphere::AbstractAtmosphere,
    )
    
    (;temp_grid, humid_grid, temp_virt_grid) = diagn.grid_variables
    (;temp_virt) = diagn.dynamics_variables
    μ = atmosphere.μ_virt_temp

    @inbounds for ij in eachgridpoint(temp_virt_grid, temp_grid, humid_grid)
        temp_virt_grid[ij] = temp_grid[ij]*(1 + μ*humid_grid[ij])
    end
    # TODO check that doing a non-linear virtual temperature in grid-point space
    # but a linear virtual temperature in spectral space to avoid another transform
    # does not cause any problems. Alternative do the transform or have a linear
    # virtual temperature in both grid and spectral space
    # spectral!(temp_virt,temp_virt_grid,S)
end

"""
$(TYPEDSIGNATURES)
Virtual temperature in grid-point space: For the PrimitiveDry temperature
and virtual temperature are the same (humidity=0). Just copy over the arrays."""
function virtual_temperature!(  diagn::DiagnosticVariablesLayer,
                                temp::LowerTriangularMatrix,
                                model::PrimitiveDry)
    
    (;temp_grid, temp_virt_grid) = diagn.grid_variables
    (;temp_virt) = diagn.dynamics_variables

    copyto!(temp_virt_grid,temp_grid)
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
    (;temp_virt) = diagn.dynamics_variables
    (;temp) = progn.timesteps[lf]
    copyto!(temp_virt,temp)
end

# function barrier
function linear_virtual_temperature!(  
    diagn::DiagnosticVariablesLayer,
    progn::PrognosticLayerTimesteps,
    model::PrimitiveWet,
    lf::Integer,
)
    linear_virtual_temperature!(diagn,progn,model.atmosphere,lf)
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
    
    (;temp_virt) = diagn.dynamics_variables
    μ = atmosphere.μ_virt_temp
    Tₖ = diagn.temp_average[]   
    (;temp,humid) = progn.timesteps[lf]

    @. temp_virt = temp + (Tₖ*μ)*humid
end