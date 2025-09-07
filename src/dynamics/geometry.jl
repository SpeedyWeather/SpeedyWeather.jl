abstract type AbstractGeometry end
export Geometry

"""
$(TYPEDSIGNATURES)
Construct Geometry struct containing parameters and arrays describing an iso-latitude grid `<:AbstractGrid`
and the vertical levels. Pass on `SpectralGrid` to calculate the following fields
$(TYPEDFIELDS)
"""
@kwdef struct Geometry{
    NF,
    Grid,
    VectorType,
    VectorFloat64Type,
} <: AbstractGeometry

    "SpectralGrid that defines spectral and grid resolution"
    spectral_grid::SpectralGrid

    "resolution parameter nlat_half of Grid, # of latitudes on one hemisphere (incl Equator)"
    nlat_half::Int = spectral_grid.nlat_half

    "maximum number of longitudes (at/around Equator)"
    nlon_max::Int = get_nlon_max(Grid, nlat_half)

    "number of latitude rings"
    nlat::Int = spectral_grid.nlat

    "number of vertical levels"
    nlayers::Int = spectral_grid.nlayers

    "total number of horizontal grid points"
    npoints::Int = spectral_grid.npoints

    # ARRAYS OF LANGITUDES/LONGITUDES
    "array of longitudes in degrees (0...360˚), empty for non-full grids"
    lond::VectorFloat64Type = get_lond(Grid, nlat_half)

    "array of latitudes in degrees (90˚...-90˚)"
    latd::VectorFloat64Type = get_latd(Grid, nlat_half)

    "array of latitudes in radians (π...-π)"
    lat::VectorType = get_lat(Grid, nlat_half)

    "array of colatitudes in radians (0...π)"
    colat::VectorType = get_colat(Grid, nlat_half)

    "longitude (0˚...360˚) for each grid point in ring order"
    londs::VectorType = get_londlatds(Grid, nlat_half)[1]

    "latitude (-90˚...˚90) for each grid point in ring order"
    latds::VectorType = get_londlatds(Grid, nlat_half)[2]

    "longitude (0...2π) for each grid point in ring order"
    lons::VectorType = RingGrids.get_lonlats(Grid, nlat_half)[1]

    "latitude (-π/2...π/2) for each grid point in ring order"
    lats::VectorType = RingGrids.get_lonlats(Grid, nlat_half)[2]

    "sin of latitudes"
    sinlat::VectorType = sind.(latd)

    "cos of latitudes"
    coslat::VectorType = cosd.(latd)

    "= 1/cos(lat)"
    coslat⁻¹::VectorType = 1 ./ coslat

    "= cos²(lat)"
    coslat²::VectorType = coslat.^2

    "# = 1/cos²(lat)"
    coslat⁻²::VectorType = 1 ./ coslat²

    # VERTICAL SIGMA COORDINATE σ = p/p0 (fraction of surface pressure)
    "σ at half levels, σ_k+1/2"
    σ_levels_half::VectorType = default_sigma_coordinates(nlayers)

    "σ at full levels, σₖ"
    σ_levels_full::VectorType = 0.5*(σ_levels_half[2:end] + σ_levels_half[1:end-1])

    "σ level thicknesses, σₖ₊₁ - σₖ"
    σ_levels_thick::VectorType = σ_levels_half[2:end] - σ_levels_half[1:end-1]

    "log of σ at full levels, include surface (σ=1) as last element"
    ln_σ_levels_full::VectorType = log.(vcat(σ_levels_full, 1))

    "Full to half levels interpolation"
    full_to_half_interpolation::VectorType = σ_interpolation_weights(σ_levels_full, σ_levels_half)
end

"""
$(TYPEDSIGNATURES)
Generator function for `Geometry` struct based on `spectral_grid`."""
function Geometry(SG::SpectralGrid; vertical_coordinates=SigmaCoordinates(SG.nlayers))

    (; nlayers) = SG
    error_message = "nlayers=$(SG.nlayers) does not match length nlayers="*
        "$(vertical_coordinates.nlayers) in spectral_grid.vertical_coordinates."
    @assert nlayers == vertical_coordinates.nlayers error_message

    (; NF, grid, VectorType) = SG
    VectorTypeFloat64 = SG.ArrayType{Float64, 1}

    (; σ_half) = vertical_coordinates
    return Geometry{NF, typeof(grid), VectorType, VectorTypeFloat64}(; spectral_grid=SG, σ_levels_half=σ_half)
end

function Base.show(io::IO, G::Geometry)
    print(io, "Geometry for $(G.spectral_grid)")
end

"""
$(TYPEDSIGNATURES)
Interpolation weights for full to half level interpolation
on sigma coordinates. Following Fortran SPEEDY documentation eq. (1)."""
function σ_interpolation_weights(
    σ_levels_full::AbstractVector,
    σ_levels_half::AbstractVector)

    weights = zero(σ_levels_full)
    nlayers = length(weights)
    nlayers == 1 && return weights     # escape early for 1 layer to avoid out-of-bounds access

    for k in 1:nlayers-1
        weights[k] = (log(σ_levels_half[k+1]) - log(σ_levels_full[k])) /
                        (log(σ_levels_full[k+1]) - log(σ_levels_full[k]))
    end
    # was log(0.99) in Fortran SPEEDY code but doesn't make sense to me
    weights[end] =  (log(σ_levels_half[nlayers+1]) - log(σ_levels_full[nlayers])) /
                        (log(σ_levels_full[nlayers]) - log(σ_levels_full[nlayers-1]))

    return weights
end


"""
$(TYPEDSIGNATURES)
Given a vector in column defined at full levels, do a linear interpolation in
log(σ) to calculate its values at half-levels, skipping top (k=1/2), extrapolating to bottom (k=nlayers+1/2).
"""
function vertical_interpolate!(
    A_half::Vector,             # quantity A on half levels (excl top)
    A_full::Vector,             # quantity A on full levels
    G::Geometry,
)
    nlayers = length(A_half)
    weights = G.full_to_half_interpolation

    # full levels contain one more for surface
    # TODO this is currently confusing because the surface fluxes use full[end]
    # as surface value which is technically on half levels though!
    @boundscheck nlayers <= length(A_full) || throw(BoundsError)
    @boundscheck nlayers <= length(weights) || throw(BoundsError)

    # For A at each full level k, compute A at the half-level below, i.e. at the boundary
    # between the full levels k and k+1. Fortran SPEEDY documentation eq. (1)
    for k = 1:nlayers-1
        A_half[k] = A_full[k] + weights[k]*(A_full[k+1] - A_full[k])
    end

    # Compute the values at the surface separately
    A_half[nlayers] = A_full[nlayers] + weights[nlayers]*(A_full[nlayers] - A_full[nlayers-1])

    return nothing
end