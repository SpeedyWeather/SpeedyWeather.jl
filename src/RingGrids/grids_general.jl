abstract type AbstractGrid{T} <: AbstractVector{T} end

# all AbstractGrids have their grid points stored in a vector field `data`
# propagate length, size, getindex, setindex! for that
Base.length(G::AbstractGrid) = length(G.data)
Base.size(G::AbstractGrid) = size(G.data)
@inline function Base.getindex(G::AbstractGrid,k::Integer)
    @boundscheck 0 < k <= length(G.data) || throw(BoundsError(G,k))
    @inbounds r = G.data[k]
    return r
end

@inline Base.setindex!(G::AbstractGrid,x,k::Integer) = setindex!(G.data,x,k)

# with ranges
@inline Base.getindex(G::AbstractGrid,r::AbstractRange) = G.data[r]
@inline Base.setindex!(G::AbstractGrid,x::AbstractVector,r::AbstractRange) = setindex!(G.data,x,r)

function grids_match(A::AbstractGrid,B::AbstractGrid)
    length(A) == length(B) && return _grids_match(typeof(A),typeof(B))
    return false
end

function _grids_match(A::Type{<:AbstractGrid},B::Type{<:AbstractGrid})
    # throws an error for non-parametric types...
    return A.name.wrapper == B.name.wrapper
end

"""
    abstract type AbstractGrid{T} <: AbstractVector{T} end

The abstract supertype for all spatial grids on the sphere supported by SpeedyWeather.jl.
Every new grid has to be of the form

    abstract type AbstractGridClass{T} <: AbstractGrid{T} end
    struct MyNewGrid{T} <: AbstractGridClass{T}
        data::Vector{T}     # all grid points unravelled into a vector
        nlat_half::Int      # resolution: latitude rings on one hemisphere (Equator incl)
    end
    
`MyNewGrid` should belong to a grid class like `AbstractFullGrid`, `AbstractOctahedralGrid` or
`AbstractHEALPixGrid` (that already exist but you may introduce a new class of grids) that share
certain features such as the number of longitude points per latitude ring and indexing, but may
have different latitudes or offset rotations. Each new grid `Grid` (or grid class) then has to
implement the following methods (as an example, see octahedral.jl)
    
Fundamental grid properties
    get_npoints         # total number of grid points
    nlat_odd            # does the grid have an odd number of latitude rings?
    get_nlat            # total number of latitude rings
    get_nlat_half       # number of latitude rings on one hemisphere incl Equator
    
Indexing
    get_nlon_max        # maximum number of longitudes points (at the Equator)
    get_nlon_per_ring   # number of longitudes on ring j
    each_index_in_ring  # a unit range that indexes all longitude points on a ring
    
Coordinates
    get_colat           # vector of colatitudes (radians)
    get_colatlon        # vectors of colatitudes, longitudes (both radians)

Spectral truncation
    truncation_order    # linear, quadratic, cubic = 1,2,3 for grid 
    get_truncation      # spectral truncation given a grid resolution
    get_resolution      # grid resolution given a spectral truncation

Quadrature weights and solid angles
    get_quadrature_weights  # = sinθ Δθ for grid points on ring j for meridional integration
    get_solid_angle         # = sinθ Δθ Δϕ, solid angle of grid points on ring j
"""
AbstractGrid

# Define methods that are universally applicable to any G<:AbstractGrid here
# generator functions for grid
Base.zeros(::Type{Grid},nlat_half::Integer) where {Grid<:AbstractGrid{T}} where T =
    Grid(zeros(T,get_npoints(Grid,nlat_half)),nlat_half)
# use Float64 if not provided
Base.zeros(::Type{Grid},nlat_half::Integer) where {Grid<:AbstractGrid} = zeros(Grid{Float64},nlat_half)
# zero element of an AbstractGrid instance grid by packing a zero(::Vector) into grid
Base.zero(grid::Grid) where {Grid<:AbstractGrid} = Grid(zero(grid.data))

# initialise with ones
Base.ones(::Type{Grid},nlat_half::Integer) where {Grid<:AbstractGrid{T}} where T =
    Grid(ones(T,get_npoints(Grid,nlat_half)),nlat_half)
Base.ones(::Type{Grid},nlat_half::Integer) where {Grid<:AbstractGrid} = ones(Grid{Float64},nlat_half)

# in case type parameter T in Grid{T} is provided
function (::Type{Grid})(::UndefInitializer,nlat_half::Integer) where {Grid<:AbstractGrid{T}} where T
    return Grid(Vector{T}(undef,get_npoints(Grid,nlat_half)),nlat_half)
end

# use Float64 if not
function (::Type{Grid})(::UndefInitializer,nlat_half::Integer) where {Grid<:AbstractGrid}
    return Grid(Vector{Float64}(undef,get_npoints(Grid,nlat_half)),nlat_half)
end

# randn initializer, use Float64 if T in Grid{T} not provided
Base.randn(::Type{Grid},nlat_half::Integer) where {Grid<:AbstractGrid} = randn(Grid{Float64},nlat_half)
Base.randn(::Type{Grid},nlat_half::Integer) where {Grid<:AbstractGrid{T}} where T = 
    Grid(randn(T,get_npoints(Grid,nlat_half)),nlat_half)

# rand initializer, use Float64 if T in Grid{T} not provided
Base.rand(::Type{Grid},nlat_half::Integer) where {Grid<:AbstractGrid} = rand(Grid{Float64},nlat_half)
Base.rand(::Type{Grid},nlat_half::Integer) where {Grid<:AbstractGrid{T}} where T = 
    Grid(rand(T,get_npoints(Grid,nlat_half)),nlat_half)

# truncation is the spectral truncation corresponding to size of grid and lin/quad/cubic truncation
get_resolution(grid::AbstractGrid) = get_nlat_half(grid)
get_nlat_half(grid::AbstractGrid) = grid.nlat_half

# does the grid have an odd number of latitudes?
nlat_odd(grid::AbstractGrid) = nlat_odd(typeof(grid))

# get total number of latitude rings, *(nlat_half > 0) to return 0 for nlat_half = 0
get_nlat(Grid::Type{<:AbstractGrid},nlat_half::Integer) = 2nlat_half - nlat_odd(Grid)*(nlat_half > 0)
get_nlat(grid::Grid) where {Grid<:AbstractGrid} = get_nlat(Grid,grid.nlat_half)

# get total number of grid pointst
get_npoints(grid::Grid) where {Grid<:AbstractGrid} = get_npoints(Grid,grid.nlat_half)

# coordinates
get_latdlonds(grid::Grid) where {Grid<:AbstractGrid} = get_latdlonds(Grid,grid.nlat_half)

function get_latdlonds(Grid::Type{<:AbstractGrid},nlat_half::Integer)
    colats,lons = get_colatlons(Grid,nlat_half)     # colatitudes, longitudes in radians
    latds = colats                                  # flat copy rename before conversion
    londs = lons
    latds .= π/2 .- colats                          # colatiudes to latitudes in radians
    latds .= latds .* (360/2π)                      # now in degrees 90˚...-90˚
    londs .*= (360/2π)

    return latds, londs
end

function get_latlons(Grid::Type{<:AbstractGrid},nlat_half::Integer)
    colats,lons = get_colatlons(Grid,nlat_half)     # colatitudes, longitudes in radians
    lats = colats                                   # flat copy rename before conversion
    lats .= π/2 .- colats                           # colatiudes to latitudes in radians

    return lats, lons
end

get_lat(grid::Grid) where {Grid<:AbstractGrid} = get_lat(Grid,grid.nlat_half)
get_latd(grid::Grid) where {Grid<:AbstractGrid} = get_latd(Grid,grid.nlat_half)
get_lond(grid::Grid) where {Grid<:AbstractGrid} = get_lond(Grid,grid.nlat_half)
get_lon(grid::Grid) where {Grid<:AbstractGrid} = get_lon(Grid,grid.nlat_half)
get_colat(grid::Grid) where {Grid<:AbstractGrid} = get_colat(Grid,grid.nlat_half)
get_colatlons(grid::Grid) where {Grid<:AbstractGrid} = get_colatlons(Grid,grid.nlat_half)

function get_lat(Grid::Type{<:AbstractGrid},nlat_half::Integer)
    return π/2 .- get_colat(Grid,nlat_half)
end

function get_latd(Grid::Type{<:AbstractGrid},nlat_half::Integer)
    return get_lat(Grid,nlat_half) * (360/2π)
end

# only defined for full grids, empty vector as fallback
get_lon(::Type{<:AbstractGrid},nlat_half::Integer) = Float64[]
get_lond(::Type{<:AbstractGrid},nlat_half::Integer) = Float64[]

"""
    i = each_index_in_ring(grid,j)

UnitRange `i` to access data on grid `grid` on ring `j`."""
function each_index_in_ring(grid::Grid,j::Integer) where {Grid<:AbstractGrid}
    return each_index_in_ring(Grid,j,grid.nlat_half)
end

"""
    ijs = eachgridpoint(grid)

UnitRange `ijs` to access each grid point on grid `grid`."""
eachgridpoint(grid::AbstractGrid) = Base.OneTo(get_npoints(grid))
function eachgridpoint(grid1::Grid,grids::Grid...) where {Grid<:AbstractGrid}
    n = length(grid1)
    Base._all_match_first(X->length(X),n,grid1,grids...) || throw(BoundsError)
    return eachgridpoint(grid1)
end

"""
$(TYPEDSIGNATURES)
Vector{UnitRange} `rings` to loop over every ring of grid `grid`
and then each grid point per ring. To be used like

    rings = eachring(grid)
    for ring in rings
        for ij in ring
            grid[ij]"""
eachring(grid::AbstractGrid) = eachring(typeof(grid),grid.nlat_half)

function eachring(Grid::Type{<:AbstractGrid},nlat_half::Integer)
    rings = Vector{UnitRange{Int}}(undef,get_nlat(Grid,nlat_half))
    each_index_in_ring!(rings,Grid,nlat_half)
    return rings
end

"""
$(TYPEDSIGNATURES)
Same as `eachring(grid)` but performs a bounds check to assess that all grids
in `grids` are of same size."""
function eachring(grid1::Grid,grids::Grid...) where {Grid<:AbstractGrid}
    @inline 
    n = length(grid1)
    Base._all_match_first(X->length(X),n,grid1,grids...) || throw(BoundsError)
    return eachring(grid1)
end

"""
$(TYPEDSIGNATURES)
Obtain ring index j from gridpoint ij and Vector{UnitRange} describing rind indices
as obtained from eachring(::Grid)"""
function whichring(ij::Integer,rings::Vector{UnitRange{Int}})
    @boundscheck 0 < ij <= rings[end][end] || throw(BoundsError)
    j = 1
    @inbounds while ij > rings[j][end]
        j += 1
    end
    return j
end

"""
$(TYPEDSIGNATURES)
Returns a vector `nlons` for the number of longitude points per latitude ring, north to south.
Provide grid `Grid` and its resolution parameter `nlat_half`. For both_hemisphere==false only
the northern hemisphere (incl Equator) is returned."""
function get_nlons(Grid::Type{<:AbstractGrid},nlat_half::Integer;both_hemispheres::Bool=false)
    n = both_hemispheres ? get_nlat(Grid,nlat_half) : nlat_half
    return [get_nlon_per_ring(Grid,nlat_half,j) for j in 1:n]
end

get_nlon_max(grid::Grid) where {Grid<:AbstractGrid} = get_nlon_max(Grid,grid.nlat_half)