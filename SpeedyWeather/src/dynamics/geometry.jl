abstract type AbstractGeometry <: AbstractModelComponent end
export Geometry

"""
$(TYPEDSIGNATURES)
Construct Geometry struct containing parameters and arrays describing an iso-latitude grid `<:AbstractGrid`
and the vertical levels. Pass on `SpectralGrid` to calculate the following fields
$(TYPEDFIELDS)
"""
@kwdef struct Geometry{
        SpectralGridType,   # <:Union{SpectralGrid, Nothing}, the latter only matter inside GPU kernels
        RefValueNF,         # <:Union{Base.RefValue{NF}, CUDA.RefValue{NF}}
        VectorIntType,
        VectorType,
        IntType,            # <: Integer
    } <: AbstractGeometry

    "SpectralGrid that defines spectral and grid resolution"
    spectral_grid::SpectralGridType

    "Resolution parameter nlat_half of Grid, # of latitudes on one hemisphere (incl Equator)"
    nlat_half::IntType = spectral_grid.nlat_half

    "Maximum number of longitudes (at/around Equator)"
    nlon_max::IntType = get_nlon_max(spectral_grid.Grid, nlat_half)

    "Number of latitude rings"
    nlat::IntType = spectral_grid.nlat

    "Number of vertical levels"
    nlayers::IntType = spectral_grid.nlayers

    "Total number of horizontal grid points"
    npoints::IntType = spectral_grid.npoints

    "Planet's radius [m], set from model.planet during initialize!"
    radius::RefValueNF = RefValueNF(DEFAULT_RADIUS)

    # INDEXING
    "Latitude ring index j for every grid point ij"
    whichring::VectorIntType = RingGrids.whichring(spectral_grid.grid)

    # ARRAYS OF LANGITUDES/LONGITUDES
    "Array of longitudes in degrees (0...360˚), empty for non-full grids"
    lond::VectorType = get_lond(spectral_grid.Grid, nlat_half)

    "Array of latitudes in degrees (90˚...-90˚)"
    latd::VectorType = get_latd(spectral_grid.Grid, nlat_half)

    "Array of latitudes in radians (π...-π)"
    lat::VectorType = get_lat(spectral_grid.Grid, nlat_half)

    "Array of colatitudes in radians (0...π)"
    colat::VectorType = get_colat(spectral_grid.Grid, nlat_half)

    "Longitude (0˚...360˚) for each grid point in ring order"
    londs::VectorType = get_londlatds(spectral_grid.Grid, nlat_half)[1]

    "Latitude (-90˚...˚90) for each grid point in ring order"
    latds::VectorType = get_londlatds(spectral_grid.Grid, nlat_half)[2]

    "Longitude (0...2π) for each grid point in ring order"
    lons::VectorType = RingGrids.get_lonlats(spectral_grid.Grid, nlat_half)[1]

    "Latitude (-π/2...π/2) for each grid point in ring order"
    lats::VectorType = RingGrids.get_lonlats(spectral_grid.Grid, nlat_half)[2]

    "sin of latitudes"
    sinlat::VectorType = sind.(latd)

    "cos of latitudes"
    coslat::VectorType = cosd.(latd)

    "= 1/cos(lat)"
    coslat⁻¹::VectorType = 1 ./ coslat

    "= cos²(lat)"
    coslat²::VectorType = coslat .^ 2

    "= 1/cos²(lat)"
    coslat⁻²::VectorType = 1 ./ coslat²

    # VERTICAL SIGMA COORDINATE σ = p/p0 (fraction of surface pressure)
    "σ at half levels, σ_k+1/2"
    σ_levels_half::VectorType = default_sigma_coordinates(nlayers)

    "σ at full levels, σₖ"
    σ_levels_full::VectorType = 0.5 * (σ_levels_half[2:end] + σ_levels_half[1:(end - 1)])

    "σ level thicknesses, σₖ₊₁ - σₖ"
    σ_levels_thick::VectorType = σ_levels_half[2:end] - σ_levels_half[1:(end - 1)]
end

Adapt.@adapt_structure Geometry

"""
$(TYPEDSIGNATURES)
Generator function for `Geometry` struct based on `spectral_grid`."""
function Geometry(SG::SpectralGrid; vertical_coordinates = SigmaCoordinates(SG.nlayers))

    (; nlayers) = SG
    error_message = "nlayers=$(SG.nlayers) does not match length nlayers=" *
        "$(vertical_coordinates.nlayers) in spectral_grid.vertical_coordinates."
    @assert nlayers == vertical_coordinates.nlayers error_message

    (; NF, VectorIntType, VectorType) = SG
    (; σ_half) = vertical_coordinates
    return Geometry{typeof(SG), Base.RefValue{NF}, VectorIntType, VectorType, typeof(nlayers)}(; spectral_grid = SG, σ_levels_half = σ_half)
end

function Base.show(io::IO, G::Geometry)
    return print(io, "Geometry for $(G.spectral_grid)")
end

# take over radius from model.planet
function initialize!(geometry::Geometry, model::AbstractModel)
    geometry.radius[] = model.planet.radius
    return geometry
end