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

function virtual_temperature!(
    column::ColumnVariables,
    model::PrimitiveEquation,
)
    (; temp, temp_virt, humid) = column
    μ = model.atmosphere.μ_virt_temp

    @. temp_virt = temp*(1 + μ*humid)
    return nothing
end

function virtual_temperature!(
    column::ColumnVariables,
    model::PrimitiveDry,
)
    (; temp, temp_virt) = column
    @. temp_virt = temp             # temp = temp_virt for PrimitiveDry
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

"""
$(TYPEDSIGNATURES)
Calculates a linearised virtual temperature Tᵥ as

    Tᵥ = T + Tₖμq

With absolute temperature T, layer-average temperarture Tₖ (computed in temperature_average!),
specific humidity q and

    μ = (1-ξ)/ξ, ξ = R_dry/R_vapour.
    
in spectral space."""
function linear_virtual_temperature!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    model::PrimitiveEquation,
    lf::Integer,
)
    (; temp_virt) = diagn.dynamics
    μ = model.atmosphere.μ_virt_temp
    (; temp_average) = diagn
    temp = progn.temp[lf]
    humid = progn.humid[lf]

    # TODO check that doing a non-linear virtual temperature in grid-point space
    # but a linear virtual temperature in spectral space to avoid another transform
    # does not cause any problems. Alternative do the transform or have a linear
    # virtual temperature in both grid and spectral space
    # spectral!(temp_virt, temp_virt_grid, S)

    for k in eachmatrix(temp_virt, temp, humid)
        Tₖ = temp_average[k]
        for lm in eachharmonic(temp_virt, temp, humid)
            temp_virt[lm, k] = temp[lm, k] + (Tₖ*μ)*humid[lm, k]
        end
    end
end