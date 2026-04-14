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

function variables(::Geopotential)
    return (
        GridVariable(:geopotential, Grid3D(), desc = "Geopotential", units = "m^2/s^2"),
        # TODO only the grid should be necessary, remove this when adapting the geopotential calculation
        # in the dynamical core to only need grid geopotential
        DynamicsVariable(:geopotential, Spectral3D(), desc = "Geopotential", units = "m^2/s^2"),
    )
end

"""
$(TYPEDSIGNATURES)
Precomputes constants for the vertical integration of the geopotential, defined as

`Φ_{k+1/2} = Φ_{k+1} + R*T_{k+1}*(ln(p_{k+1}) - ln(p_{k+1/2}))` (half levels)
`Φ_k = Φ_{k+1/2} + R*T_k*(ln(p_{k+1/2}) - ln(p_k))` (full levels)

Same formula but `k → k-1/2`."""
function initialize!(G::Geopotential, model::PrimitiveEquation)
    (; Δp_geopot_half, Δp_geopot_full) = G
    (; R_dry) = model.atmosphere
    (; σ_levels_full, σ_levels_half) = model.geometry
    nlayers = length(σ_levels_full)

    # 1. integration onto half levels
    # used for: Φ_{k+1/2} = Φ_{k+1} + R*T_{k+1}*(ln(p_{k+1}) - ln(p_{k+1/2}))
    @. Δp_geopot_half[2:nlayers] = R_dry * log(σ_levels_full[2:nlayers] / σ_levels_half[2:nlayers])

    # 2. integration onto full levels (same formula but k -> k-1/2)
    # used for: Φ_k = Φ_{k+1/2} + R*T_k*(ln(p_{k+1/2}) - ln(p_k))
    @. Δp_geopot_full = R_dry * log(σ_levels_half[2:(nlayers + 1)] / σ_levels_full)
    return nothing
end

# calculate geopotential in grid-point space for parameterizations
function geopotential!(
        vars::Variables,
        model::PrimitiveEquation,
    )
    T = vars.grid.temperature
    Φ = vars.grid.geopotential

    # use zero scratch for humidity in dry models to not distinguish in kernels below
    vars.scratch.grid.a .= 0
    q = haskey(vars.grid, :humidity) ? vars.grid.humidity : vars.scratch.grid.a

    (; orography) = model.orography
    g = model.planet.gravity
    G = model.geopotential
    (; orography) = model.orography
    (; atmosphere) = model

    arch = architecture(T)
    if typeof(arch) <: GPU
        launch!(arch, LinearWorkOrder, (size(T, 1),), geopotential_kernel!, Φ, T, q, orography, g, G, atmosphere)
    else
        geopotential_cpu!(Φ, T, q, orography, g, G, atmosphere)
    end
    return nothing
end

function geopotential_cpu!(Φ, T, q, orography, g, G, atmosphere)
    nlayers = size(T, 2)
    @inbounds for ij in eachgridpoint(Φ)
        geopotential_compute!(
            ij, Φ, T, q, orography, g, G.Δp_geopot_half, G.Δp_geopot_full, nlayers, atmosphere
        )
    end
    return nothing
end

@kernel function geopotential_kernel!(Φ, T, q, orography, g, G, atmosphere)
    ij = @index(Global, Linear)
    nlayers = size(T, 2)
    @inbounds geopotential_compute!(
        ij, Φ, T, q, orography, g, G.Δp_geopot_half, G.Δp_geopot_full, nlayers, atmosphere
    )
end

@propagate_inbounds function geopotential_compute!(
        ij, geopotential, temp, humid, orography,
        gravity, Δp_geopot_half, Δp_geopot_full, nlayers, atmosphere
    )
    # bottom layer
    Tᵥ = virtual_temperature(temp[ij, nlayers], humid[ij, nlayers], atmosphere)
    geopotential[ij, nlayers] = gravity * orography[ij] + Tᵥ * Δp_geopot_full[nlayers]

    # OTHER FULL LAYERS, integrate two half-layers from bottom to top
    for k in (nlayers - 1):-1:1
        Tᵥ_below = Tᵥ
        Tᵥ = virtual_temperature(temp[ij, k], humid[ij, k], atmosphere)
        geopotential[ij, k] = geopotential[ij, k + 1] + Tᵥ_below * Δp_geopot_half[k + 1] + Tᵥ * Δp_geopot_full[k]
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Compute spectral geopotential `geopot` from spectral temperature `temp`
and spectral surface geopotential `geopot_surf` (orography*gravity)."""
function geopotential!(
        vars::Variables,
        G::Geopotential,
        orography::AbstractOrography,
    )
    Tᵥ = vars.dynamics.virtual_temperature
    Φ = vars.dynamics.geopotential              # spectral geopotential to fill
    Φₛ = orography.surface_geopotential         # = orography*gravity
    (; Δp_geopot_half, Δp_geopot_full) = G      # = R*Δlnp either on half or full levels
    nlayers = size(Φ, 2)                        # number of vertical levels

    @boundscheck nlayers == length(Δp_geopot_full) || throw(BoundsError)

    # for PrimitiveDry virtual temperature = absolute temperature here
    # note these are not anomalies here as they are only in grid-point fields

    # BOTTOM FULL LAYER
    # TODO: broadcasting with LTA issue here
    Φ.data[:, nlayers] .= Φₛ.data .+ Tᵥ.data[:, nlayers] .* Δp_geopot_full[nlayers:nlayers]

    # OTHER FULL LAYERS, integrate two half-layers from bottom to top
    arch = architecture(Φ)
    launch!(
        arch, SpectralWorkOrder, (size(Φ, 1),), geopotential_spectral_kernel!,
        Φ, Tᵥ, Δp_geopot_half, Δp_geopot_full, nlayers
    )
    return nothing
end

@kernel inbounds = true function geopotential_spectral_kernel!(
        geopotential,               # Output: spectral geopotential
        virtual_temperature,        # Input: spectral virtual temperature
        Δp_geopot_half,             # Input: integration constant for half levels
        Δp_geopot_full,             # Input: integration constant for full levels
        nlayers,                    # Input: number of vertical layers
    )
    lm = @index(Global, Linear)     # global index: harmonic lm

    Φ = geopotential                # for brevity
    Tᵥ = virtual_temperature

    # Integrate from bottom to top over all layers except bottom (already computed)
    for k in (nlayers - 1):-1:1
        Φ_k½ = Φ[lm, k + 1] + Tᵥ[lm, k + 1] * Δp_geopot_half[k + 1]     # 1st half layer integration
        Φ[lm, k] = Φ_k½ + Tᵥ[lm, k] * Δp_geopot_full[k]                 # 2nd onto full layer
    end
end

"""
$(TYPEDSIGNATURES)
calculates the geopotential in the ShallowWaterModel as g*η,
i.e. gravity times the interface displacement (field `pres`)"""
function geopotential!(vars::Variables, planet::AbstractPlanet)
    # note this only works in 2D models with nlayers=1 otherwise gepotential is actually 3D and not just 2D x 1
    # one would need to use lta/field_view(gepotential, :, nlayers) to make it 2D again
    (; geopotential, η) = vars.grid
    geopotential .= planet.gravity .* η
    return geopotential
end
