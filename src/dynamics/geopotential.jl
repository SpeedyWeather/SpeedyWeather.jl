abstract type AbstractGeopotential <: AbstractModelComponent end

export Geopotential

"""Geopotential calculation component for primitive equation models. Precomputes
integration constants for vertical hydrostatic integration from temperature to
geopotential on both full and half pressure levels using the hypsometric equation.
Fields are $(TYPEDFIELDS)"""
@kwdef struct Geopotential{VectorType} <: AbstractGeopotential
    "[DERIVED] Used to compute the geopotential on half layers"
    Δp_geopot_half::VectorType

    "[DERIVED] Used to compute the geopotential on full layers"
    Δp_geopot_full::VectorType
end

Adapt.@adapt_structure Geopotential

# generator
Geopotential(SG::SpectralGrid) = Geopotential(
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

# calculate geopotential in grid-point space for parameterizations
function geopotential!(
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)
    temp = diagn.grid.temp_grid
    humid = diagn.grid.humid_grid
    g = model.planet.gravity
    G = model.geopotential
    (; geopotential) = diagn.grid
    (; orography) = model.orography

    # this triggers a fallback to use virtual temperature = absolute temperature for PrimitiveDry
    A = model isa PrimitiveWet ? model.atmosphere : nothing

    arch = architecture(temp)

    if typeof(arch) <: GPU
        launch!(arch, LinearWorkOrder, (size(temp, 1),), geopotential_kernel!, geopotential, temp, humid, orography, g, G, A)
    else 
        geopotential_cpu!(geopotential, temp, humid, orography, g, G, A)
    end 
end

function geopotential_cpu!(geopotential, temp, humid, orography, gravity, Geopotential, atmosphere)
    nlayers = size(temp, 2)
    @inbounds for ij in eachgridpoint(geopotential)
        geopotential_compute!(ij, geopotential, temp, humid, orography, gravity, Geopotential.Δp_geopot_half, Geopotential.Δp_geopot_full, nlayers, atmosphere)
    end
end 

@kernel function geopotential_kernel!(geopotential, temp, humid, orography, gravity, Geopotential, atmosphere)
    ij = @index(Global, Linear)
    nlayers = size(temp, 2)

    @inbounds geopotential_compute!(ij, geopotential, temp, humid, orography, gravity, Geopotential.Δp_geopot_half, Geopotential.Δp_geopot_full, nlayers, atmosphere)
end

@propagate_inbounds function geopotential_compute!(ij, geopotential, temp, humid, orography, gravity, Δp_geopot_half, Δp_geopot_full, nlayers, atmosphere)
    # bottom layer
    Tᵥ = virtual_temperature(temp[ij, nlayers], humid[ij, nlayers], atmosphere)
    geopotential[ij, nlayers] = gravity*orography[ij] + Tᵥ*Δp_geopot_full[nlayers]

    # OTHER FULL LAYERS, integrate two half-layers from bottom to top
    for k in nlayers-1:-1:1
        Tᵥ_below = Tᵥ
        Tᵥ = virtual_temperature(temp[ij, k], humid[ij, k], atmosphere)
        geopotential[ij, k] = geopotential[ij, k+1] + Tᵥ_below*Δp_geopot_half[k+1] + Tᵥ*Δp_geopot_full[k]
    end
end



"""
$(TYPEDSIGNATURES)
Compute spectral geopotential `geopot` from spectral temperature `temp`
and spectral surface geopotential `geopot_surf` (orography*gravity)."""
function geopotential!( 
    diagn::DiagnosticVariables,
    G::Geopotential,
    orography::AbstractOrography,
)
    (; temp_virt) = diagn.dynamics
    geopot = diagn.dynamics.geopotential                # spectral geopotential to fill
    geopot_surf = orography.surface_geopotential        # = orography*gravity
    (; Δp_geopot_half, Δp_geopot_full) = G              # = R*Δlnp either on half or full levels
    (; nlayers) = diagn                                 # number of vertical levels

    @boundscheck nlayers == length(Δp_geopot_full) || throw(BoundsError)

    # for PrimitiveDry virtual temperature = absolute temperature here
    # note these are not anomalies here as they are only in grid-point fields
    
    # BOTTOM FULL LAYER
    geopot[:, nlayers] .= geopot_surf .+ temp_virt[:, nlayers] .* Δp_geopot_full[nlayers]

    # OTHER FULL LAYERS, integrate two half-layers from bottom to top
    arch = architecture(geopot)
    launch!(arch, SpectralWorkOrder, (size(geopot, 1),), geopotential_spectral_kernel!,
            geopot, temp_virt, Δp_geopot_half, Δp_geopot_full, nlayers)
end

@kernel inbounds = true function geopotential_spectral_kernel!(
    geopot,                     # Output: spectral geopotential
    temp_virt,                  # Input: spectral virtual temperature
    @Const(Δp_geopot_half),     # Input: integration constant for half levels
    @Const(Δp_geopot_full),     # Input: integration constant for full levels
    @Const(nlayers),            # Input: number of vertical layers
)
    lm = @index(Global, Linear)  # global index: harmonic lm

    # Integrate from bottom to top over all layers except bottom (already computed)
    for k in nlayers-1:-1:1
        geopot_k½ = geopot[lm, k+1] + temp_virt[lm, k+1] * Δp_geopot_half[k+1]  # 1st half layer integration
        geopot[lm, k] = geopot_k½ + temp_virt[lm, k] * Δp_geopot_full[k]        # 2nd onto full layer
    end
end

"""
$(TYPEDSIGNATURES)
calculates the geopotential in the ShallowWaterModel as g*η,
i.e. gravity times the interface displacement (field `pres`)"""
function geopotential!( diagn::DiagnosticVariables,
                        planet::AbstractPlanet)

    # note this only works in 2D models with nlayers=1 otherwise gepotential is actually 3D and not just 2D x 1
    # one would need to use lta/field_view(gepotential, :, nlayers) to make it 2D again
    (; geopotential, pres_grid) = diagn.grid
    geopotential .= planet.gravity .* pres_grid
    return nothing
end