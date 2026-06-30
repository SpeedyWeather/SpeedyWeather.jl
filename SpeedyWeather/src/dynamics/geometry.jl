abstract type AbstractGeometry <: AbstractModelComponent end
export Geometry

"""
$(TYPEDSIGNATURES)
Construct Geometry struct containing parameters and arrays describing an iso-latitude grid `<:AbstractGrid`
and the vertical levels. Pass on `SpectralGrid` to calculate the following fields
$(TYPEDFIELDS)
"""
@kwdef struct Geometry{
        SpectralGridType,   # <: Union{SpectralGrid, Nothing}, the latter only matter inside GPU kernels
        IntType,            # <: Integer
        RefValueNF,         # <: Union{Base.RefValue{NF}, CUDA.RefValue{NF}}
        VectorIntType,
        VC,
        VT1,                # These seemingly redundant extra VectorTypes are needed for Reactant compat 
        VTLat,              # which needs different types for different sizes of the array
        VTGrid,
        VTVert1,
        VTVert2,
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
    lond::VT1 = get_lond(spectral_grid.Grid, nlat_half)

    "Array of latitudes in degrees (90˚...-90˚)"
    latd::VTLat = get_latd(spectral_grid.Grid, nlat_half)

    "Array of latitudes in radians (π...-π)"
    lat::VTLat = get_lat(spectral_grid.Grid, nlat_half)

    "Array of colatitudes in radians (0...π)"
    colat::VTLat = get_colat(spectral_grid.Grid, nlat_half)

    "Longitude (0˚...360˚) for each grid point in ring order"
    londs::VTGrid = get_londlatds(spectral_grid.Grid, nlat_half)[1]

    "Latitude (-90˚...˚90) for each grid point in ring order"
    latds::VTGrid = get_londlatds(spectral_grid.Grid, nlat_half)[2]

    "Longitude (0...2π) for each grid point in ring order"
    lons::VTGrid = RingGrids.get_lonlats(spectral_grid.Grid, nlat_half)[1]

    "Latitude (-π/2...π/2) for each grid point in ring order"
    lats::VTGrid = RingGrids.get_lonlats(spectral_grid.Grid, nlat_half)[2]

    "sin of latitudes"
    sinlat::VTLat = sind.(latd)

    "cos of latitudes"
    coslat::VTLat = cosd.(latd)

    "= 1/cos(lat)"
    coslat⁻¹::VTLat = 1 ./ coslat

    "= cos²(lat)"
    coslat²::VTLat = coslat .^ 2

    "= 1/cos²(lat)"
    coslat⁻²::VectorType = 1 ./ coslat²

    "Vertical coordinates used"
    vertical_coordinates::VC

    # VERTICAL SIGMA COORDINATE σ = p/p0 (fraction of surface pressure)
    "σ at half levels, σ_k+1/2"
    σ_levels_half::VTVert1

    "σ at full levels, σₖ"
    σ_levels_full::VTVert2 = (σ_levels_half[2:end] + σ_levels_half[1:(end - 1)]) / 2

    "σ level thicknesses, σₖ₊₁ - σₖ"
    σ_levels_thick::VTVert2 = σ_levels_half[2:end] - σ_levels_half[1:(end - 1)]
end

Adapt.@adapt_structure Geometry

"""
$(TYPEDSIGNATURES)
Generator function for `Geometry` struct based on `spectral_grid`."""
function Geometry(SG::SpectralGrid; vertical_coordinates = SigmaCoordinates(SG))

    (; nlayers) = SG
    error_message = "nlayers=$(SG.nlayers) does not match length nlayers=" *
        "$(get_nlayers(vertical_coordinates)) in spectral_grid.vertical_coordinates."
    @assert nlayers == get_nlayers(vertical_coordinates) error_message

    (; NF, VectorIntType, VectorType) = SG
    σ_half = get_σ_half(vertical_coordinates)
    return Geometry{typeof(SG), typeof(nlayers), Base.RefValue{NF}, VectorIntType, typeof(vertical_coordinates), VectorType, VectorType, VectorType, VectorType, VectorType}(;
        spectral_grid = SG,
        vertical_coordinates,
        σ_levels_half = σ_half,
    )
end

function Base.show(io::IO, G::Geometry)
    (; grid, nlat, npoints, nlayers) = G.spectral_grid
    Grid = nonparametric_type(grid)

    params = "{Spectrum{...}, $Grid{...}}"
    println(io, styled"{warning:Geometry} for SpectralGrid{note:$params}")
    println(io, styled"├ {info:Grid}: $nlat-ring $Grid, $npoints grid points")
    print(io, styled"└ {info:Vertical}: $(G.vertical_coordinates)")
    return nothing
end

# take over radius from model.planet
function initialize!(geometry::Geometry, model::AbstractModel)
    geometry.radius[] = model.planet.radius

    if hasproperty(geometry.vertical_coordinates, :reference_pressure)
        p_ref_coord = geometry.vertical_coordinates.reference_pressure
        p_ref_atmos = model.atmosphere.reference_pressure

        p_ref_coord != p_ref_atmos &&
            @warn "Reference pressure of vertical coordinates and atmosphere differ. "*
                "$p_ref_coord Pa vs $p_ref_atmos Pa"
    end

    return geometry
end

@inline pressure(k::Integer, surface_pressure::Number, geometry::Geometry) = 
    pressure(k, surface_pressure, geometry.vertical_coordinate)

@inline pressure_thickness(k::Integer, surface_pressure::Number, geometry::Geometry) = 
    pressure_thickness(k, surface_pressure, geometry.vertical_coordinate)