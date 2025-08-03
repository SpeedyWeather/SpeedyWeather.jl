isreduced(::Type{<:AbstractGrid}) = true
isreduced(::Type{<:AbstractFullGrid}) = false
isreduced(grid) = isreduced(typeof(grid))
isfull(x) = ~isreduced(x)

function Base.show(io::IO, grid::AbstractGrid)
    Grid_ = nonparametric_type(grid)
    nlat_half = get_nlat_half(grid)
    nlat = get_nlat(grid)
    npoints = get_npoints(grid)
    full_or_reduced = isfull(grid) ? "full" : "reduced"
    res = sqrt(4π/npoints)*360/2π
    digits = round(Int, log10(1/res)) + 2
    average_resolution = Printf.@sprintf("%.*f˚", digits, res)

    println(io, "$nlat-ring $Grid_")
    println(io, "├ nlat_half=$nlat_half ($npoints points, ~$average_resolution, $full_or_reduced)")
    print(io,   "└ architecture: $(typeof(grid.architecture))")
end

## TYPES
grid_type(::Type{Grid}) where {Grid<:AbstractGrid} = nonparametric_type(Grid) 
full_grid_type(grid::AbstractGrid) = full_grid_type(typeof(grid))

"""$(TYPEDSIGNATURES) For any instance of `AbstractGrid` type its n-dimensional type
(*Grid{T, N, ...} returns *Array) but without any parameters `{T, N, ArrayType}`"""
Architectures.nonparametric_type(grid::AbstractGrid) = nonparametric_type(typeof(grid))

"""$(TYPEDSIGNATURES) Resolution paraemeters `nlat_half` of a `grid`.
Number of latitude rings on one hemisphere, Equator included."""
get_nlat_half(grid::AbstractGrid) = grid.nlat_half

"""$(TYPEDSIGNATURES) True for a `grid` that has an odd number of
latitude rings `nlat` (both hemispheres)."""
nlat_odd(grid::AbstractGrid) = nlat_odd(typeof(grid))

# get total number of latitude rings, *(nlat_half > 0) to return 0 for nlat_half = 0
get_nlat(Grid::Type{<:AbstractGrid}, nlat_half::Integer) = 2nlat_half - nlat_odd(Grid)*(nlat_half > 0)

"""$(TYPEDSIGNATURES) Get number of latitude rings, pole to pole."""
get_nlat(grid::AbstractGrid) = get_nlat(typeof(grid), grid.nlat_half)

"""$(TYPEDSIGNATURES) Total number of grid points in all dimensions of `grid`.
Equivalent to length of the underlying data array."""
get_npoints(grid::Grid, args...) where {Grid<:AbstractGrid} = get_npoints(Grid, grid.nlat_half, args...)
get_npoints(G::Type{<:AbstractGrid}, nlat_half::Integer, k::Integer...) = prod(k) * get_npoints(G, nlat_half)
get_quadrature_weights(grid::AbstractGrid) = get_quadrature_weights(typeof(grid), grid.nlat_half)

"""$(TYPEDSIGNATURES) Size of the matrix of the horizontal grid if representable as such (not all grids)."""
matrix_size(grid::Grid) where {Grid<:AbstractGrid} = matrix_size(Grid, get_nlat_half(grid))

# CONSTRUCTORS
"""$(TYPEDSIGNATURES) Create a new `grid` of type `Grid` with resolution parameter `nlat_half`.
`architecture` is the device type (CPU/GPU). Precomputes the ring indices `rings`."""
function (::Type{Grid})(nlat_half::Integer, architecture=DEFAULT_ARCHITECTURE()) where {Grid<:AbstractGrid}
    Grid_ = nonparametric_type(Grid)                # strip away parameters of type, obtain from arguments
    rings = eachring(Grid, nlat_half)               # precompute indices to access the variable-length rings
    w = whichring(Grid, nlat_half, rings)           # precompute ring indices for each grid point
    return on_architecture(architecture, Grid_(nlat_half, architecture, rings, w))
end

# also allow to construct a field with Grid(data)
function (::Type{Grid})(data::AbstractArray; input_as=Vector, kwargs...) where {Grid<:AbstractGrid}
    return Grid(data, input_as, kwargs...)          # make input_as a positional argument
end

# change the architecture of a grid, keep all other fields 
function (::Type{Grid})(grid::Grid, architecture::AbstractArchitecture) where {Grid<:AbstractGrid}
    Grid_ = nonparametric_type(Grid)                # strip away parameters of type, obtain from arguments
    return Grid_(grid.nlat_half, architecture, adapt(array_type(architecture), grid.rings), adapt(array_type(architecture), grid.whichring))
end

function (::Type{Grid})(
    data::AbstractArray,
    input_as::Type{Vector};
    architecture=DEFAULT_ARCHITECTURE(),
) where {Grid<:AbstractGrid}
    # create a grid based on the size of data
    npoints = size(data, 1)
    nlat_half = get_nlat_half(Grid, npoints)
    grid = Grid(nlat_half, architecture)
    return Field(data, grid)
end

function (::Type{Grid})(
    data::AbstractArray,
    input_as::Type{Matrix};
    architecture=DEFAULT_ARCHITECTURE(),
) where {Grid<:AbstractGrid}
    npoints = size(data, 1)*size(data, 2)
    nlat_half = get_nlat_half(Grid, npoints)
    grid = Grid(nlat_half, architecture)
    data_flat = reshape(data, :, size(data)[3:end]...)
    return Field(data_flat, grid)
end

## COORDINATES

"""$(TYPEDSIGNATURES) Longitudes (degrees, 0-360˚E), latitudes (degrees, 90˚N to -90˚N) for
every (horizontal) grid point in `grid` in ring order (0-360˚E then north to south)."""
get_londlatds(grid::Grid) where {Grid<:AbstractGrid} = get_londlatds(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Longitudes (radians, 0-2π), latitudes (degrees, π/2 to -π/2) for
every (horizontal) grid point in `grid` in ring order (0-360˚E then north to south)."""
get_lonlats(grid::Grid) where {Grid<:AbstractGrid} = get_lonlats(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Longitudes (radians, 0-2π), colatitudes (degrees, 0 to π) for
every (horizontal) grid point in `grid` in ring order (0-360˚E then north to south)."""
get_loncolats(grid::Grid) where {Grid<:AbstractGrid} = get_loncolats(Grid, grid.nlat_half)

# radians
function get_lonlats(Grid::Type{<:AbstractGrid}, nlat_half::Integer)
    londs, latds = get_londlatds(Grid, nlat_half)  # longitudes, latitudes in degrees
    return londs * π/180, latds * π/180             # to radians
end

# radians and colatitudes
function get_loncolats(Grid::Type{<:AbstractGrid}, nlat_half::Integer)
    lons, lats = get_lonlats(Grid, nlat_half)
    colats = π/2 .- lats
    return lons, colats
end

"""$(TYPEDSIGNATURES) Latitude (radians) for each ring in `grid`, north to south."""
get_lat(grid::Grid) where {Grid<:AbstractGrid} = get_lat(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Latitude (degrees) for each ring in `grid`, north to south."""
get_latd(grid::Grid) where {Grid<:AbstractGrid} = get_latd(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Colatitudes (radians) for each ring in `grid`, north to south."""
get_colat(grid::Grid) where {Grid<:AbstractGrid} = get_colat(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Longitude (degrees). Full grids only."""
get_lond(grid::Grid) where {Grid<:AbstractGrid} = get_lond(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Longitude (radians). Full grids only."""
get_lon(grid::Grid) where {Grid<:AbstractGrid} = get_lon(Grid, grid.nlat_half)

get_nlon_max(grid::Grid) where {Grid<:AbstractGrid} = get_nlon_max(Grid, grid.nlat_half)
get_nlons(grid::AbstractGrid) = get_nlons(typeof(grid), grid.nlat_half)

get_lat(Grid::Type{<:AbstractGrid}, nlat_half::Integer) = get_latd(Grid, nlat_half) * (π/180)
get_colat(Grid::Type{<:AbstractGrid}, nlat_half::Integer) = π/2 .- get_lat(Grid, nlat_half)

"""$(TYPEDSIGNATURES)
Returns a vector `nlons` for the number of longitude points per latitude ring, north to south.
Provide grid `Grid` and its resolution parameter `nlat_half`. For keyword argument
`both_hemispheres=false` only the northern hemisphere (incl Equator) is returned."""
function get_nlons(Grid::Type{<:AbstractGrid}, nlat_half::Integer; both_hemispheres::Bool=true)
    n = both_hemispheres ? get_nlat(Grid, nlat_half) : nlat_half
    return [get_nlon_per_ring(Grid, nlat_half, j) for j in 1:n]
end

"""$(TYPEDSIGNATURES)
Number of longitude points per latitude ring `j`."""
function get_nlon_per_ring(grid::AbstractGrid, j::Integer)
    return get_nlon_per_ring(typeof(grid), grid.nlat_half, j)
end

get_solid_angles(grid::AbstractGrid) = get_solid_angles(typeof(grid), grid.nlat_half)

## ITERATOR

"""$(TYPEDSIGNATURES)
Vector{UnitRange} `rings` to loop over every ring of `grid`
and then each grid point per ring. To be used like

    rings = eachring(grid)
    for ring in rings
        for ij in ring
            field[ij]

Accesses precomputed `grid.rings`."""
@inline eachring(grid::AbstractGrid) = grid.rings
@inline eachring(grid::AbstractGrid{<:GPU}) = Vector(grid.rings) # on GPU transfer indices back to CPU first 

"""$(TYPEDSIGNATURES)
Computes the ring indices `i0:i1` for start and end of every longitudinal point
on a given ring `j` of `Grid` at resolution `nlat_half`. Used to loop
over rings of a grid. These indices are also precomputed in every `grid.rings`."""
function eachring(Grid::Type{<:AbstractGrid}, nlat_half::Integer)
    nlat = get_nlat(Grid, nlat_half)
    rings = Vector{UnitRange{Int}}(undef, nlat)    # allocate
    each_index_in_ring!(rings, Grid, nlat_half)    # calculate iteratively
    return rings                                                    
end

"""$(TYPEDSIGNATURES) Same as `eachring(grid)` but performs a bounds check to assess
that all `grids` according to `grids_match` (non-parametric grid type, nlat_half and length)."""
function eachring(grid1::AbstractGrid, grids::AbstractGrid...)
    grids_match(grid1, grids...) || throw(DimensionMismatch(grid1, grids...))
    return eachring(grid1)
end

function Base.DimensionMismatch(grid1::AbstractGrid, grids::AbstractGrid...)
    nlat = get_nlat_half(grid1)
    s = "grids do not match; $nlat-ring $(nonparametric_type(grid1))"
    for grid in grids
        nlat = get_nlat_half(grid)
        s *= ", $nlat-ring $(nonparametric_type(grid))"
    end
    return DimensionMismatch(s)
end

# equality and comparison, somehow needed as not covered by broadcasting
Base.:(==)(G1::AbstractGrid, G2::AbstractGrid) = grids_match(G1, G2)

"""$(TYPEDSIGNATURES) True if both `A` and `B` are of the same nonparametric grid type
(e.g. OctahedralGaussianArray, regardless type parameter `T` or underyling array type `ArrayType`)
and of same resolution (`nlat_half`) and total grid points (`length`). Sizes of `(4,)` and `(4,1)`
would match for example, but `(8,1)` and `(4,2)` would not (`nlat_half` not identical)."""
function grids_match(A::AbstractGrid, B::AbstractGrid)
    horizontal_match = get_nlat_half(A) == get_nlat_half(B)
    type_match = grids_match(typeof(A), typeof(B))
    return horizontal_match && type_match
end

"""$(TYPEDSIGNATURES) True if both `A` and `B` are of the same nonparametric grid type."""
function grids_match(A::Type{<:AbstractGrid}, B::Type{<:AbstractGrid})
    return nonparametric_type(A) == nonparametric_type(B)
end

"""$(TYPEDSIGNATURES) True if all grids `A, B, C, ...` provided as arguments
match according to `grids_match` wrt to `A` (and therefore all)."""
function grids_match(A::AbstractGrid, Bs::AbstractGrid...)
    match = true    # single grid A always matches itself
    for B in Bs     # check for all matching respectively with A
        match &= grids_match(A, B)
    end
    return match
end

"""$(TYPEDSIGNATURES) UnitRange to access data on grid `grid` on ring `j`."""
each_index_in_ring(grid::AbstractGrid, j::Integer) = grid.rings[j]

""" $(TYPEDSIGNATURES) UnitRange to access each horizontal grid point on grid `grid`.
For a `NxM` (`N` horizontal grid points, `M` vertical layers) `OneTo(N)` is returned."""
eachgridpoint(grid::AbstractGrid) = Base.OneTo(get_npoints(grid))

""" $(TYPEDSIGNATURES) Like `eachgridpoint(::AbstractGrid)` but checks `grids` match."""
function eachgridpoint(grid1::AbstractGrid, grids::AbstractGrid...)
    grids_match(grid1, grids...) || throw(DimensionMismatch(grid1, grids...))
    return eachgridpoint(grid1)
end

"""$(TYPEDSIGNATURES) Obtain ring index `j` from gridpoint `ij` and `rings`
describing rind indices as obtained from `eachring(::Grid)`"""
function whichring(ij::Integer, rings)
    @boundscheck 0 < ij <= rings[end][end] || throw(BoundsError)
    j = 1
    @inbounds while ij > rings[j][end]
        j += 1
    end
    return j
end

whichring(ij::Integer, grid::AbstractGrid) = whichring(ij, grid.rings)

"""$(TYPEDSIGNATURES) Vector of ring indices for every grid point in `grid`."""
function whichring(Grid::Type{<:AbstractGrid}, nlat_half, rings)
    w = zeros(Int, get_npoints(Grid, nlat_half))
    @inbounds for (j, ring) in enumerate(rings)
        w[ring] .= j
    end
    return w
end

whichring(grid::AbstractGrid) = whichring(typeof(grid), grid.nlat_half, grid.rings)
whichring(Grid::Type{<:AbstractGrid}, nlat_half::Integer) = whichring(Grid, nlat_half, eachring(Grid, nlat_half))

# for architectures / adapt 
Architectures.ismatching(grid::AbstractGrid, array_type::Type{<:AbstractArray}) = ismatching(grid.architecture, array_type)
Architectures.ismatching(grid::AbstractGrid, array::AbstractArray) = ismatching(grid.architecture, typeof(array))

Architectures.architecture(grid::AbstractGrid) = grid.architecture

function Architectures.on_architecture(arch::AbstractArchitecture, grid::Grid) where Grid<:AbstractGrid 
    Grid_ = nonparametric_type(Grid)
    return Grid_(grid, arch)
end 