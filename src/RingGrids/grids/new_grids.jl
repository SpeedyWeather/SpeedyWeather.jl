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
    truncation_order    # linear, quadratic, cubic = 1, 2, 3 for grid 
    get_truncation      # spectral truncation given a grid resolution
    get_resolution      # grid resolution given a spectral truncation

Quadrature weights and solid angles
    get_quadrature_weights  # = sinθ Δθ for grid points on ring j for meridional integration
    get_solid_angle         # = sinθ Δθ Δϕ, solid angle of grid points on ring j
"""