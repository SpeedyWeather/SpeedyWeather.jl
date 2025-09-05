abstract type AbstractGeopotential <: AbstractModelComponent end

export Geopotential
@kwdef struct Geopotential{NF, VectorType} <: AbstractGeopotential

    "Number of vertical layers"
    nlayers::Int

    "Used to compute the geopotential on half layers"
    Δp_geopot_half::VectorType = zeros(NF, nlayers)

    "Used to compute the geopotential on full layers"
    Δp_geopot_full::VectorType = zeros(NF, nlayers)
end

Geopotential(SG::SpectralGrid) = Geopotential{SG.NF, SG.VectorType}(; nlayers=SG.nlayers)

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
    (; Δp_geopot_half, Δp_geopot_full, nlayers) = geopotential
    (; R_dry) = model.atmosphere
    (; σ_levels_full, σ_levels_half) = model.geometry

    # 1. integration onto half levels
    for k in 1:nlayers-1               # k is full level index, 1=top, nlayers=bottom
        # used for: Φ_{k+1/2} = Φ_{k+1} + R*T_{k+1}*(ln(p_{k+1}) - ln(p_{k+1/2}))
        Δp_geopot_half[k+1] = R_dry*log(σ_levels_full[k+1]/σ_levels_half[k+1])
    end

    # 2. integration onto full levels (same formula but k -> k-1/2)
    for k in 1:nlayers
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
    (; geopot, temp_virt) = diagn.dynamics
    (; geopot_surf) = orography                          # = orography*gravity
    (; Δp_geopot_half, Δp_geopot_full) = geopotential    # = R*Δlnp either on half or full levels
    (; nlayers) = diagn                                 # number of vertical levels

    @boundscheck nlayers == length(Δp_geopot_full) || throw(BoundsError)

    # for PrimitiveDry virtual temperature = absolute temperature here
    # note these are not anomalies here as they are only in grid-point fields
    
    # BOTTOM FULL LAYER
    local k::Int = nlayers
    @inbounds for lm in eachharmonic(geopot, geopot_surf, temp_virt)
        geopot[lm, k] = geopot_surf[lm] + temp_virt[lm, k]*Δp_geopot_full[k]
    end

    # OTHER FULL LAYERS, integrate two half-layers from bottom to top
    @inbounds for k in nlayers-1:-1:1
        for lm in eachharmonic(geopot, temp_virt)
            geopot_k½ = geopot[lm, k+1] + temp_virt[lm, k+1]*Δp_geopot_half[k+1] # 1st half layer integration
            geopot[lm, k] = geopot_k½  + temp_virt[lm, k]*Δp_geopot_full[k]      # 2nd onto full layer
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

    nlayers = length(geopot)
    (; Δp_geopot_half, Δp_geopot_full) = G  # = R*Δlnp either on half or full levels

    @boundscheck length(temp) >= nlayers || throw(BoundsError)
    @boundscheck length(Δp_geopot_full) >= nlayers || throw(BoundsError)
    @boundscheck length(Δp_geopot_half) >= nlayers || throw(BoundsError)

    # bottom layer
    geopot[nlayers] = geopot_surf + temp[nlayers]*Δp_geopot_full[end]

    # OTHER FULL LAYERS, integrate two half-layers from bottom to top
    @inbounds for k in nlayers-1:-1:1
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
function geopotential!( diagn::DiagnosticVariables,
                        pres::LowerTriangularArray,
                        planet::AbstractPlanet)
    (; geopot) = diagn.dynamics

    # don't use broadcasting as geopot will have size Nxnlayers but pres N
    # [lm] indexing bypasses the incompatible sizes (necessary for primitive models)
    for lm in eachindex(geopot, pres)
        geopot[lm] = pres[lm] * planet.gravity
    end
end