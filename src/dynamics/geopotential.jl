"""
$(TYPEDSIGNATURES)
Precomputes constants for the vertical integration of the geopotential, defined as

`Φ_{k+1/2} = Φ_{k+1} + R*T_{k+1}*(ln(p_{k+1}) - ln(p_{k+1/2}))` (half levels)
`Φ_k = Φ_{k+1/2} + R*T_k*(ln(p_{k+1/2}) - ln(p_k))` (full levels)

Same formula but `k → k-1/2`."""
function initialize_geopotential(   σ_levels_full::Vector,
                                    σ_levels_half::Vector,
                                    R_dry::Real)

    nlev = length(σ_levels_full)    # number of vertical levels
    @assert nlev+1 == length(σ_levels_half) "σ half levels must have length nlev+1"

    Δp_geopot_half = zeros(nlev)    # allocate arrays
    Δp_geopot_full = zeros(nlev)

    # 1. integration onto half levels
    for k in 1:nlev-1               # k is full level index, 1=top, nlev=bottom
        # used for: Φ_{k+1/2} = Φ_{k+1} + R*T_{k+1}*(ln(p_{k+1}) - ln(p_{k+1/2}))
        Δp_geopot_half[k+1] = R_dry*log(σ_levels_full[k+1]/σ_levels_half[k+1])
    end

    # 2. integration onto full levels (same formula but k -> k-1/2)
    for k in 1:nlev
        # used for: Φ_k = Φ_{k+1/2} + R*T_k*(ln(p_{k+1/2}) - ln(p_k))
        Δp_geopot_full[k] = R_dry*log(σ_levels_half[k+1]/σ_levels_full[k])
    end

    return Δp_geopot_half, Δp_geopot_full
end

"""
$(TYPEDSIGNATURES)
Compute spectral geopotential `geopot` from spectral temperature `temp`
and spectral surface geopotential `geopot_surf` (orography*gravity).
"""
function geopotential!( 
    diagn::DiagnosticVariables,
    O::AbstractOrography,       # contains surface geopotential
    C::DynamicsConstants,       # contains precomputed layer-thickness arrays
)

    (;geopot_surf) = O                      # = orography*gravity
    (;Δp_geopot_half, Δp_geopot_full) = C   # = R*Δlnp either on half or full levels
    (;nlev) = diagn                         # number of vertical levels

    @boundscheck nlev == length(Δp_geopot_full) || throw(BoundsError)

    # for PrimitiveDry virtual temperature = absolute temperature here
    # note these are not anomalies here as they are only in grid-point fields
    
    # BOTTOM FULL LAYER
    temp = diagn.layers[end].dynamics_variables.temp_virt
    geopot = diagn.layers[end].dynamics_variables.geopot
    
    @inbounds for lm in eachharmonic(geopot,geopot_surf,temp)
        geopot[lm] = geopot_surf[lm] + temp[lm]*Δp_geopot_full[end]
    end

    # OTHER FULL LAYERS, integrate two half-layers from bottom to top
    @inbounds for k in nlev-1:-1:1
        temp_k    = diagn.layers[k].dynamics_variables.temp_virt
        temp_k1   = diagn.layers[k+1].dynamics_variables.temp_virt
        geopot_k  = diagn.layers[k].dynamics_variables.geopot
        geopot_k1 = diagn.layers[k+1].dynamics_variables.geopot

        for lm in eachharmonic(temp_k,temp_k1,geopot_k,geopot_k1)
            geopot_k½ = geopot_k1[lm] + temp_k1[lm]*Δp_geopot_half[k+1] # 1st half layer integration
            geopot_k[lm] = geopot_k½  + temp_k[lm]*Δp_geopot_full[k]    # 2nd onto full layer
        end      
    end
end

"""
$(TYPEDSIGNATURES)
Calculate the geopotential based on `temp` in a single column.
This exclues the surface geopotential that would need to be added to the returned vector.
Function not used in the dynamical core but for post-processing and analysis."""
function geopotential!( geopot::AbstractVector,
                        temp::AbstractVector,
                        C::DynamicsConstants,
                        geopot_surf::Real = 0)

    nlev = length(geopot)
    (;Δp_geopot_half, Δp_geopot_full) = C  # = R*Δlnp either on half or full levels

    @boundscheck length(temp) >= nlev || throw(BoundsError)
    @boundscheck length(Δp_geopot_full) >= nlev || throw(BoundsError)
    @boundscheck length(Δp_geopot_half) >= nlev || throw(BoundsError)

    # bottom layer
    geopot[nlev] = geopot_surf + temp[nlev]*Δp_geopot_full[end]

    # OTHER FULL LAYERS, integrate two half-layers from bottom to top
    @inbounds for k in nlev-1:-1:1
        geopot[k] = geopot[k+1] + temp[k+1]*Δp_geopot_half[k+1] + temp[k]*Δp_geopot_full[k]
    end

    return geopot
end

function geopotential(  temp::Vector,
                        C::DynamicsConstants) 
    geopot = zero(temp)
    geopotential!(geopot,temp,C)
    return geopot
end

"""
$(TYPEDSIGNATURES)
calculates the geopotential in the ShallowWaterModel as g*η,
i.e. gravity times the interface displacement (field `pres`)"""
function geopotential!( diagn::DiagnosticVariablesLayer,
                        pres::LowerTriangularMatrix,
                        C::DynamicsConstants)
    (;gravity) = C
    (;geopot) = diagn.dynamics_variables
    geopot .= pres*gravity
end 

# function barrier
function virtual_temperature!(  diagn::DiagnosticVariablesLayer,
                                temp::LowerTriangularMatrix,    # only needed for dispatch compat with DryCore
                                model::PrimitiveWet)
    virtual_temperature!(diagn,temp,model.constants)
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
    constants::DynamicsConstants,
    )
    
    (;temp_grid, humid_grid, temp_virt_grid) = diagn.grid_variables
    (;temp_virt) = diagn.dynamics_variables
    μ = constants.μ_virt_temp

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
    linear_virtual_temperature!(diagn,progn,model.constants,lf)
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
                                        constants::DynamicsConstants,
                                        lf::Int)
    
    (;temp_virt) = diagn.dynamics_variables
    μ = constants.μ_virt_temp
    Tₖ = diagn.temp_average[]   
    (;temp,humid) = progn.timesteps[lf]

    @. temp_virt = temp + (Tₖ*μ)*humid
end