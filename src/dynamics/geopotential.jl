abstract type AbstractGeopotential <: AbstractModelComponent end

export Geopotential

"""Geopotential calculation component for primitive equation models. Precomputes
integration constants for vertical hydrostatic integration from temperature to
geopotential on both full and half pressure levels using the hypsometric equation.
Fields are $(TYPEDFIELDS)"""
@kwdef struct Geopotential{VectorType} <: AbstractGeopotential
    # number of vertical levels
    nlayers::Int 

    "Used to compute the geopotential on half layers"
    Δp_geopot_half::VectorType = zero(VectorType(undef, nlayers))

    "Used to compute the geopotential on full layers"
    Δp_geopot_full::VectorType = zero(VectorType(undef, nlayers))
end

Adapt.@adapt_structure Geopotential

# generator
Geopotential(SG::SpectralGrid) = Geopotential(
    SG.nlayers,
    on_architecture(SG.architecture, zeros(SG.NF, SG.nlayers)),
    on_architecture(SG.architecture, zeros(SG.NF, SG.nlayers))
)

# function barrier to unpack only model components needed
function initialize!(
    geopotential::Geopotential,
    model::PrimitiveEquation
)
    model_parameters = (atmosphere=model.atmosphere, geometry=model.geometry)
    initialize!(geopotential, model_parameters)
end

"""
$(TYPEDSIGNATURES)
Precomputes constants for the vertical integration of the geopotential, defined as

`Φ_{k+1/2} = Φ_{k+1} + R*T_{k+1}*(ln(p_{k+1}) - ln(p_{k+1/2}))` (half levels)
`Φ_k = Φ_{k+1/2} + R*T_k*(ln(p_{k+1/2}) - ln(p_k))` (full levels)

Same formula but `k → k-1/2`."""
function initialize!(G::Geopotential, model)
    (; Δp_geopot_half, Δp_geopot_full) = G
    (; R_dry) = model.atmosphere
    (; σ_levels_full, σ_levels_half) = model.geometry
    nlayers = length(σ_levels_full)

    # 1. integration onto half levels
    # used for: Φ_{k+1/2} = Φ_{k+1} + R*T_{k+1}*(ln(p_{k+1}) - ln(p_{k+1/2}))
    @. Δp_geopot_half[2:nlayers] = R_dry * log(σ_levels_full[2:nlayers] / σ_levels_half[2:nlayers])

    # 2. integration onto full levels (same formula but k -> k-1/2)
    # used for: Φ_k = Φ_{k+1/2} + R*T_k*(ln(p_{k+1/2}) - ln(p_k))
    @. Δp_geopot_full = R_dry * log(σ_levels_half[2:nlayers+1] / σ_levels_full)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Compute spectral geopotential `geopotential` from spectral temperature `temp`
and spectral `surface_geopotential` (orography*gravity)."""
function geopotential!( 
    diagn::DiagnosticVariables,
    G::Geopotential,
    orography::AbstractOrography,
)
    (; geopotential, temp_virt) = diagn.dynamics
    (; surface_geopotential) = orography                          # = orography*gravity
    (; Δp_geopot_half, Δp_geopot_full) = G    # = R*Δlnp either on half or full levels
    (; nlayers) = diagn                                 # number of vertical levels

    @boundscheck nlayers == length(Δp_geopot_full) || throw(BoundsError)

    # for PrimitiveDry virtual temperature = absolute temperature here
    # note these are not anomalies here as they are only in grid-point fields
    
    # BOTTOM FULL LAYER
    local k::Int = nlayers
    @inbounds for lm in eachharmonic(geopotential, surface_geopotential, temp_virt)
        geopotential[lm, k] = surface_geopotential[lm] + temp_virt[lm, k]*Δp_geopot_full[k]
    end

    # OTHER FULL LAYERS, integrate two half-layers from bottom to top
    @inbounds for k in nlayers-1:-1:1
        for lm in eachharmonic(geopotential, temp_virt)
            geopot_k½ = geopotential[lm, k+1] + temp_virt[lm, k+1]*Δp_geopot_half[k+1] # 1st half layer integration
            geopotential[lm, k] = geopot_k½  + temp_virt[lm, k]*Δp_geopot_full[k]      # 2nd onto full layer
        end      
    end
end

"""
$(TYPEDSIGNATURES)
Calculate the geopotential based on `temp` in the column `ij`.
Used for parameterizations that need the geopotential in a single column,
convection and vertical diffusion."""
function geopotential!(
    ij,                         # horizontal grid point index
    geopotential,                     # ::AbstractField3D, geopotential to be filled
    temp,                       # ::AbstractField3D, temperature field (virtual temperature for PrimitiveWet)
    orography,                  # ::AbstractField2D, orography field
    gravity,                    # ::Real, gravity [m/s^2]
    G::Geopotential,            # precomputed integration constants
)
    (; nlayers, Δp_geopot_half, Δp_geopot_full) = G  # = R*Δlnp either on half or full levels

    # bottom layer
    geopotential[ij, end] = gravity*orography[ij] + temp[ij, nlayers]*Δp_geopot_full[end]

    # OTHER FULL LAYERS, integrate two half-layers from bottom to top
    @inbounds for k in nlayers-1:-1:1
        geopotential[ij, k] = geopotential[ij, k+1] + temp[ij, k+1]*Δp_geopot_half[k+1] + temp[ij, k]*Δp_geopot_full[k]
    end
end

# TODO provide a convenience function to calculate geopotential?

"""
$(TYPEDSIGNATURES)
calculates the geopotential in the ShallowWaterModel as g*η,
i.e. gravity times the interface displacement (field `pres`)"""
function geopotential!( diagn::DiagnosticVariables,
                        pres::LowerTriangularArray,
                        planet::AbstractPlanet)

    # note this only works in 2D models with nlayers=1 otherwise gepotential is actually 3D and not just 2D x 1
    # one would need to use lta_view(gepotential, :, nlayers) to make it 2D again
    (; geopotential) = diagn.dynamics
    geopotential .= planet.gravity .* pres
    return nothing
end