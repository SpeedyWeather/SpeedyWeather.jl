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

"""
    abstract type AbstractGrid{T} <: AbstractVector{T} end

The abstract supertype for all spatial grids on the sphere supported by SpeedyWeather.jl.
Every new grid has to be of the form

    abstract type AbstractGridClass{T} <: AbstractGrid{T} end
    struct MyNewGrid{T} <: AbstractGridClass{T}
        data::Vector{T}     # all grid points unravelled into a vector
        n::Int              # resolution parameter, like `nlat_half`, `nside` etc
    end
    
`MyNewGrid` should belong to a grid class like `AbstractFullGrid`, `AbstractOctahedralGrid` or
`AbstractHEALPixGrid` (that already exist but you may introduce a new class of grids) that share
certain features such as the number of longitude points per latitude ring and indexing, but may
have different latitudes or offset rotations. Each new grid `Grid` (or grid class) then has to
implement the following methods (see grids.jl)
    
Fundamental grid properties
    get_nresolution     # returns the resolution parameter nlat_half, nside, etc.
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
Base.zeros(::Type{Grid},n::Integer) where {Grid<:AbstractGrid{T}} where T =
    Grid(zeros(T,get_npoints(Grid,n)),n)
# use Float64 if not provided
Base.zeros(::Type{Grid},n::Integer) where {Grid<:AbstractGrid} = zeros(Grid{Float64},n)
# zero element of an AbstractGrid instance grid by packing a zero(::Vector) into grid
Base.zero(grid::Grid) where {Grid<:AbstractGrid} = Grid(zero(grid.data))

get_truncation(grid::Grid) where {Grid<:AbstractGrid} = get_truncation(Grid,get_nresolution(grid))
get_nlat_half(::Grid,n::Integer) where {Grid<:AbstractGrid} = get_nlat_half(Grid,n) 
nlat_odd(grid::AbstractGrid) = nlat_odd(typeof(grid))
get_nlat_half(::Type{<:AbstractGrid},nlat_half::Integer) = nlat_half
get_nlat(Grid::Type{<:AbstractGrid},n::Integer) = 2get_nlat_half(Grid,n) - nlat_odd(Grid)
get_nlat(grid::Grid) where {Grid<:AbstractGrid} = get_nlat(Grid,get_nresolution(grid))

function each_index_in_ring(grid::G,j::Integer) where {G<:AbstractGrid}
    return each_index_in_ring(G,j,get_nresolution(grid))
end

eachgridpoint(grid::G) where {G<:AbstractGrid} = Base.OneTo(get_npoints(G,get_nresolution(grid)))

function eachring(grid::G) where {G<:AbstractGrid}
    rings = [each_index_in_ring(grid,j) for j in 1:get_nlat(grid)]
    return rings                                            # return Vector{UnitRange}
end

function eachring(grids::Grid...) where {Grid<:AbstractGrid}
    for grid in grids
        @boundscheck length(grid) == length(grids[1]) || throw(BoundsError)
    end
    return eachring(grids[1])
end

"""
    nlons = get_nlons(  Grid::Type{<:AbstractGrid},
                        nresolution::Integer;
                        both_hemispheres::Bool=false)

Returns a vector `nlons` for the number of longitude points per latitude ring, north to south.
Provide grid `Grid` and its resolution parameter `nresolution`. For both_hemisphere==false only
the northern hemisphere (incl Equator) is returned."""
function get_nlons(Grid::Type{<:AbstractGrid},nresolution::Integer;both_hemispheres::Bool=false)
    n = both_hemispheres ? get_nlat(Grid,nresolution) : get_nlat_half(Grid,nresolution)
    return [get_nlon_per_ring(Grid,nresolution,j) for j in 1:n]
end