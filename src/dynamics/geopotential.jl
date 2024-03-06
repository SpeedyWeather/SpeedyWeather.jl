abstract type AbstractGeopotential <: AbstractModelComponent end

export Geopotential
Base.@kwdef struct Geopotential{NF} <: AbstractGeopotential
    nlev::Int
    Δp_geopot_half::Vector{NF} = zeros(NF, nlev)
    Δp_geopot_full::Vector{NF} = zeros(NF, nlev)
end

Geopotential(SG::SpectralGrid) = Geopotential{SG.NF}(; nlev=SG.nlev)

"""
$(TYPEDSIGNATURES)
Precomputes constants for the vertical integration of the geopotential, defined as

`Φ_{k+1/2} = Φ_{k+1} + R*T_{k+1}*(ln(p_{k+1}) - ln(p_{k+1/2}))` (half levels)
`Φ_k = Φ_{k+1/2} + R*T_k*(ln(p_{k+1/2}) - ln(p_k))` (full levels)

Same formula but `k → k-1/2`."""
function initialize!(
    geopotential::Geopotential,
    model::PrimitiveEquation
)
    (; Δp_geopot_half, Δp_geopot_full, nlev) = geopotential
    (; R_dry) = model.atmosphere
    (; σ_levels_full, σ_levels_half) = model.geometry

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
end

"""
$(TYPEDSIGNATURES)
Compute spectral geopotential `geopot` from spectral temperature `temp`
and spectral surface geopotential `geopot_surf` (orography*gravity)."""
function geopotential!( 
    diagn::DiagnosticVariables,
    geopotential::Geopotential,
    orography::AbstractOrography,
)

    (; geopot_surf) = orography                          # = orography*gravity
    (; Δp_geopot_half, Δp_geopot_full) = geopotential    # = R*Δlnp either on half or full levels
    (; nlev) = diagn                                     # number of vertical levels

    @boundscheck nlev == length(Δp_geopot_full) || throw(BoundsError)

    # for PrimitiveDry virtual temperature = absolute temperature here
    # note these are not anomalies here as they are only in grid-point fields
    
    # BOTTOM FULL LAYER
    temp = diagn.layers[end].dynamics_variables.temp_virt
    geopot = diagn.layers[end].dynamics_variables.geopot
    
    @inbounds for lm in eachharmonic(geopot, geopot_surf, temp)
        geopot[lm] = geopot_surf[lm] + temp[lm]*Δp_geopot_full[end]
    end

    # OTHER FULL LAYERS, integrate two half-layers from bottom to top
    @inbounds for k in nlev-1:-1:1
        temp_k    = diagn.layers[k].dynamics_variables.temp_virt
        temp_k1   = diagn.layers[k+1].dynamics_variables.temp_virt
        geopot_k  = diagn.layers[k].dynamics_variables.geopot
        geopot_k1 = diagn.layers[k+1].dynamics_variables.geopot

        for lm in eachharmonic(temp_k, temp_k1, geopot_k, geopot_k1)
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
function geopotential!( 
    geopot::AbstractVector,
    temp::AbstractVector,
    G::Geopotential,
    geopot_surf::Real = 0
)

    nlev = length(geopot)
    (; Δp_geopot_half, Δp_geopot_full) = G  # = R*Δlnp either on half or full levels

    @boundscheck length(temp) >= nlev || throw(BoundsError)
    @boundscheck length(Δp_geopot_full) >= nlev || throw(BoundsError)
    @boundscheck length(Δp_geopot_half) >= nlev || throw(BoundsError)

    # bottom layer
    geopot[nlev] = geopot_surf + temp[nlev]*Δp_geopot_full[end]

    # OTHER FULL LAYERS, integrate two half-layers from bottom to top
    @inbounds for k in nlev-1:-1:1
        geopot[k] = geopot[k+1] + temp[k+1]*Δp_geopot_half[k+1] + temp[k]*Δp_geopot_full[k]
    end
end

function geopotential(  temp::Vector,
                        G::Geopotential) 
    geopot = zero(temp)
    geopotential!(geopot, temp, G)
    return geopot
end

"""
$(TYPEDSIGNATURES)
calculates the geopotential in the ShallowWaterModel as g*η,
i.e. gravity times the interface displacement (field `pres`)"""
function geopotential!( diagn::DiagnosticVariablesLayer,
                        pres::LowerTriangularMatrix,
                        planet::AbstractPlanet)
    (; geopot) = diagn.dynamics_variables
    geopot .= pres * planet.gravity
end 