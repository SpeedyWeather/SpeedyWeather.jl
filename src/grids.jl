"""
    abstract type AbstractGrid{T} <: AbstractVector{T} end

The abstract supertype for all spatial grids on the sphere supported by SpeedyWeather.jl."""
abstract type AbstractGrid{T} <: AbstractVector{T} end

# all AbstractGrids have their grid points stored in a vector field `v`
# propagate length, size, getindex, setindex! for that
Base.length(G::AbstractGrid) = length(G.v)
Base.size(G::AbstractGrid) = size(G.v)
@inline function Base.getindex(G::AbstractGrid,k::Integer)
    @boundscheck 0 < k <= length(G.v) || throw(BoundsError(G,k))
    @inbounds r = G.v[k]
    return r
end

@inline Base.setindex!(G::AbstractGrid,x,k::Integer) = setindex!(G.v,x,k)

# with ranges
@inline Base.getindex(G::AbstractGrid,r::AbstractRange) = G.v[r]
@inline Base.setindex!(G::AbstractGrid,x::AbstractVector,r::AbstractRange) = setindex!(G.v,x,r)


"""
    abstract type AbstractFullGrid{T} <: AbstractGrid{T} end

An `AbstractFullGrid` is a horizontal grid with a constant number of longitude
points across latitude rings. Different latitudes can be used, Gaussian latitudes,
equi-angle latitdes, or others."""
abstract type AbstractFullGrid{T} <: AbstractGrid{T} end

"""
    abstract type AbstractOctahedralGrid{T} <: AbstractGrid{T} end

An `AbstractOctahedralGrid` is a horizontal grid with 16+4i longitude
points on the latitude ring i starting with i=1 around the pole.
Different latitudes can be used, Gaussian latitudes, equi-angle latitdes, or others."""
abstract type AbstractOctahedralGrid{T} <: AbstractGrid{T} end

"""
    abstract type AbstractHEALPixGrid{T} <: AbstractGrid{T} end

An `AbstractHEALPixGrid` is a horizontal grid similar to the standard HEALPixGrid,
but different latitudes can be used, the default HEALPix latitudes or others."""
abstract type AbstractHEALPixGrid{T} <: AbstractGrid{T} end

"""
    G = FullClenshawGrid{T}

A FullClenshawGrid is a regular latitude-longitude grid with an odd number of `nlat` equi-spaced
latitudes, the central latitude ring is on the Equator. The same `nlon` longitudes for every latitude ring.
The grid points are closer in zonal direction around the poles. The values of all grid points are stored
in a vector field `v` that unravels the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct FullClenshawGrid{T} <: AbstractFullGrid{T}
    v::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere (incl Equator)

    FullClenshawGrid{T}(v,nlat_half) where T = length(v) == npoints_clenshaw(nlat_half) ?
    new(v,nlat_half) : error("$(length(v))-element Vector{$(eltype(v))} cannot be used to create a "*
        "L$nlat_half ($(4nlat_half)x$(2nlat_half - 1)) FullClenshawGrid{$T}.")
end

# subtract the otherwise double-counted 4nlat_half equator points
npoints_clenshaw(nlat_half::Integer) = 8nlat_half^2 - 4nlat_half
nlat_half_clenshaw(npoints::Integer) = round(Int,1/4 + sqrt(1/16 + npoints/8))  # inverse

"""
    G = FullGaussianGrid{T}

A full Gaussian grid is a regular latitude-longitude grid that uses `nlat` Gaussian latitudes,
and the same `nlon` longitudes for every latitude ring. The grid points are closer in zonal direction
around the poles. The values of all grid points are stored in a vector field `v` that unravels
the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct FullGaussianGrid{T} <: AbstractFullGrid{T}
    v::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere

    FullGaussianGrid{T}(v,nlat_half) where T = length(v) == 8nlat_half^2 ?
    new(v,nlat_half) : error("$(length(v))-element Vector{$(eltype(v))} cannot be used to create a "*
        "F$nlat_half ($(4nlat_half)x$(2nlat_half)) FullGaussianGrid{$T}.")
end

npoints_gaussian(nlat_half::Integer) = 8nlat_half^2
nlat_half_gaussian(npoints::Integer) = round(Int,sqrt(npoints/8))

"""
    G = OctahedralGaussianGrid{T}

An Octahedral Gaussian grid that uses `nlat` Gaussian latitudes, but a decreasing number of longitude
points per latitude ring towards the poles. Starting with 20 equi-spaced longitude points (starting at 0˚E)
on the rings around the poles, each latitude ring towards the equator has consecuitively 4 more points,
one for each face of the octahedron. E.g. 20,24,28,32,...nlon-4,nlon,nlon,nlon-4,...,32,28,24,20.
The maximum number of longitue points is `nlon`. The values of all grid points are stored in a vector
field `v` that unravels the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct OctahedralGaussianGrid{T} <: AbstractOctahedralGrid{T}
    v::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere

    # check that `nlat_half` match the vector `v` length
    OctahedralGaussianGrid{T}(v,nlat_half) where T = length(v) == npoints_octahedral(nlat_half,false) ?
    new(v,nlat_half) : error("$(length(v))-element Vector{$(eltype(v))}"*
    "cannot be used to create a O$(nlat_half) OctahedralGaussianGrid{$T}.")
end

# number of points and longitudes per ring on the octahedral grid
npoints_octahedral(nlat_half::Integer,nlat_oddp::Bool) =
    nlat_oddp ? max(0,4nlat_half^2 + 32nlat_half - 16) : 4nlat_half^2 + 36nlat_half # max(0,...) needed to avoid negative array size when nlat_half==0
nlat_half_octahedral(npoints::Integer,nlat_oddp::Bool) =
    nlat_oddp ? round(Int,-4+sqrt(20 + npoints/4)) : round(Int,-9/2+sqrt((9/2)^2 + npoints/4))  # inverse
nlon_octahedral(j::Integer) = 16+4j


"""
    G = OctahedralClenshawGrid{T}

An Octahedral Clenshaw grid that uses `nlat` equi-spaced latitudes. Like FullClenshawGrid, the central
latitude ring is on the Equator. Like OctahedralGaussianGrid, the number of longitude points per
latitude ring decreases towards the poles. Starting with 20 equi-spaced longitude points (starting at 0˚E)
on the rings around the poles, each latitude ring towards the equator has consecuitively 4 more points,
one for each face of the octahedron. E.g. 20,24,28,32,...nlon-4,nlon,nlon,nlon-4,...,32,28,24,20.
The maximum number of longitue points is `nlon`. The values of all grid points are stored in a vector
field `v` that unravels the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct OctahedralClenshawGrid{T} <: AbstractOctahedralGrid{T}
    v::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere (incl Equator)

    # check that `nlat_half` match the vector `v` length
    OctahedralClenshawGrid{T}(v,nlat_half) where T = length(v) == npoints_octahedral(nlat_half,true) ?
    new(v,nlat_half) : error("$(length(v))-element Vector{$(eltype(v))}"*
    "cannot be used to create a O$(nlat_half) OctahedralClenshawGrid{$T}.")
end

"""
    H = HEALPixGrid{T}

A HEALPix grid with 12 faces, each `nside`x`nside` grid points, each covering the same area.
The values of all grid points are stored in a vector field `v` that unravels the data 0 to 360˚,
then ring by ring, which are sorted north to south."""
struct HEALPixGrid{T} <: AbstractHEALPixGrid{T}
    v::Vector{T}    # data vector, ring by ring, north to south
    nside::Int      # nside^2 is the number of pixel in each of the 12 base pixel

    HEALPixGrid{T}(v,nside) where T = length(v) == npoints_healpix(nside) ?
    new(v,nside) : error("$(length(v))-element Vector{$(eltype(v))}"*
    "cannot be used to create an H$nside HEALPixGrid{$T}.")
end

# number of points and longitudes per ring on the HEALPix grid
function nside_assert(nside::Integer)
    @assert is_power_2_or_0(nside) "HEALPixGrid: nside=$nside is not a power of 2."
    return true
end

npoints_healpix(nside::Integer) = nside_assert(nside) ? 12nside^2 : nothing
nside_healpix(npoints::Integer) = round(Int,sqrt(npoints/12))  # inverse of npoints_healpix
nlat_healpix(nside::Integer) = nside_assert(nside) ? 4nside-1 : nothing
nlon_healpix(nside::Integer,j::Integer) = nside_assert(nside) ? min(4j,4nside) : nothing
nlon_healpix(nside::Integer) = nside_assert(nside) ? 4nside : nothing

# define for all grids that the type T can be infered from the elements in data vector
# whether the resolution parameter n is provided or not (hence the ...)
FullClenshawGrid(v::AbstractVector,n::Integer...) = FullClenshawGrid{eltype(v)}(v,n...)
FullGaussianGrid(v::AbstractVector,n::Integer...) = FullGaussianGrid{eltype(v)}(v,n...)
OctahedralGaussianGrid(v::AbstractVector,n::Integer...) = OctahedralGaussianGrid{eltype(v)}(v,n...)
OctahedralClenshawGrid(v::AbstractVector,n::Integer...) = OctahedralClenshawGrid{eltype(v)}(v,n...)
HEALPixGrid(v::AbstractVector,n::Integer...) = HEALPixGrid{eltype(v)}(v,n...)

# infer resolution parameter nlat_half or nside from length of vector
FullClenshawGrid{T}(v::AbstractVector) where T = FullClenshawGrid(v,nlat_half_clenshaw(length(v)))
FullGaussianGrid{T}(v::AbstractVector) where T = FullGaussianGrid(v,nlat_half_gaussian(length(v)))
OctahedralGaussianGrid{T}(v::AbstractVector) where T = OctahedralGaussianGrid(v,nlat_half_octahedral(length(v),false))
OctahedralClenshawGrid{T}(v::AbstractVector) where T = OctahedralClenshawGrid(v,nlat_half_octahedral(length(v),true))
HEALPixGrid{T}(v::AbstractVector) where T = HEALPixGrid(v,nside_healpix(length(v)))

# convert an AbstractMatrix to the full grids
FullClenshawGrid(M::AbstractMatrix{T}) where T = FullClenshawGrid{T}(vec(M))
FullGaussianGrid(M::AbstractMatrix{T}) where T = FullGaussianGrid{T}(vec(M))

# and vice versa
Base.Matrix(G::AbstractFullGrid{T}) where T = Matrix{T}(reshape(G.v,:,2G.nlat_half - nlat_odd(G)))

# generator functions for grid
Base.zeros(::Type{FullClenshawGrid{T}},nlat_half::Integer) where T = 
                FullClenshawGrid(zeros(T,npoints_clenshaw(nlat_half)),nlat_half)
Base.zeros(::Type{FullGaussianGrid{T}},nlat_half::Integer) where T =
                FullGaussianGrid(zeros(T,npoints_gaussian(nlat_half)),nlat_half)
Base.zeros(::Type{OctahedralGaussianGrid{T}},nlat_half::Integer) where T =
                OctahedralGaussianGrid(zeros(T,npoints_octahedral(nlat_half,false)),nlat_half)
Base.zeros(::Type{OctahedralClenshawGrid{T}},nlat_half::Integer) where T =
                OctahedralClenshawGrid(zeros(T,npoints_octahedral(nlat_half,true )),nlat_half)

Base.zeros(::Type{HEALPixGrid{T}},nside::Integer) where T =
                HEALPixGrid(zeros(T,npoints_healpix(nside)),nside)

# use Float64 if not provided
Base.zeros(::Type{G},n::Integer) where {G<:AbstractGrid} = zeros(G{Float64},n)

# zero element of an AbstractGrid instance G by packing a zero(::Vector) into G
Base.zero(g::G) where {G<:AbstractGrid} = G(zero(g.v))

# MATCHING SPECTRAL TO GRID POINT RESOLUTION
truncation_order(::Type{<:FullClenshawGrid}) = 3            # cubic
truncation_order(::Type{<:FullGaussianGrid}) = 2            # quadratic
truncation_order(::Type{<:OctahedralGaussianGrid}) = 3      # cubic
truncation_order(::Type{<:OctahedralClenshawGrid}) = 3      # cubic
truncation_order(::Type{<:HEALPixGrid}) = 1                 # linear (in longitude)

get_truncation(::Type{<:FullClenshawGrid},nlat_half::Integer) = floor(Int,(4nlat_half-1)/4)
get_truncation(::Type{<:FullGaussianGrid},nlat_half::Integer) = floor(Int,(4nlat_half-1)/3)
get_truncation(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = nlat_half-1
get_truncation(::Type{<:OctahedralClenshawGrid},nlat_half::Integer) = nlat_half-1
get_truncation(::Type{<:HEALPixGrid},nside::Integer) = nside_assert(nside) ? 2nside-1 : nothing
get_truncation(grid::G) where {G<:AbstractGrid} = get_truncation(G,get_nresolution(grid))

get_resolution(::Type{<:FullClenshawGrid},trunc::Integer) = roundup_fft(ceil(Int,(4*trunc+1)/4))
get_resolution(::Type{<:FullGaussianGrid},trunc::Integer) = roundup_fft(ceil(Int,(3*trunc+1)/4))
get_resolution(::Type{<:OctahedralGaussianGrid},trunc::Integer) = roundup_fft(trunc+1)
get_resolution(::Type{<:OctahedralClenshawGrid},trunc::Integer) = roundup_fft(trunc+1)
get_resolution(::Type{<:HEALPixGrid},trunc::Integer) = roundup_fft(ceil(Int,(trunc+1)/2),small_primes=[2])

# common interface for the resolution parameter nlat_half/nside
get_nresolution(G::AbstractFullGrid) = G.nlat_half
get_nresolution(G::AbstractOctahedralGrid) = G.nlat_half
get_nresolution(G::AbstractHEALPixGrid) = G.nside

# define nlat_half for all grids (HEALPixGrid is different as it doesn't use nlat_half as resolution parameter)
get_nlat_half(::Type{<:AbstractGrid},nlat_half::Integer) = nlat_half
get_nlat_half(::Type{<:AbstractHEALPixGrid},nside::Integer) = 2nside
get_nlat_half(::G,n::Integer) where {G<:AbstractGrid} = get_nlat_half(G,n) 

# define whether there's an odd number of latitude rings for a grid
nlat_odd(::Type{<:FullGaussianGrid}) = false
nlat_odd(::Type{<:FullClenshawGrid}) = true
nlat_odd(::Type{<:OctahedralGaussianGrid}) = false
nlat_odd(::Type{<:OctahedralClenshawGrid}) = true
nlat_odd(::Type{<:AbstractHEALPixGrid}) = true
nlat_odd(grid::AbstractGrid) = nlat_odd(typeof(grid))

get_nlat(Grid::Type{<:AbstractGrid},n::Integer) = 2get_nlat_half(Grid,n) - nlat_odd(Grid)
get_nlat(grid::Grid) where {Grid<:AbstractGrid} = get_nlat(Grid,get_nresolution(grid))

# return the maxmimum number of longitude points for a grid and its resolution parameter nlat_half/nside
get_nlon(::Type{<:AbstractFullGrid},nlat_half::Integer) = 4nlat_half
get_nlon(::Type{<:AbstractOctahedralGrid},nlat_half::Integer) = nlon_octahedral(nlat_half)
get_nlon(::Type{<:AbstractHEALPixGrid},nside::Integer) = nside_assert(nside) ? nlon_healpix(nside) : nothing

# get the number of longitude points at given latitude ring j number 1...J from north to south
get_nlon_per_ring(::Type{<:AbstractFullGrid},nlat_half::Integer,j::Integer) = 4nlat_half
function get_nlon_per_ring(G::Type{<:AbstractOctahedralGrid},nlat_half::Integer,j::Integer)
    nlat = 2nlat_half-nlat_odd(G)
    @assert 0 < j <= nlat "Ring $j is outside O$nlat_half grid."
    j = j > nlat_half ? nlat - j + 1 : j      # flip north south due to symmetry
    return nlon_octahedral(j)
end
function get_nlon_per_ring(::Type{<:AbstractHEALPixGrid},nside::Integer,j::Integer)
    nlat = nlat_healpix(nside)
    @assert 0 < j <= nlat "Ring $j is outside H$nside grid."
    nlat_half = (nlat+1)÷2
    j = j > nlat_half ? nlat - j + 1 : j      # flip north south due to symmetry
    return nlon_healpix(nside,j)
end

# total number of grid points per grid type
get_npoints(::Type{<:FullGaussianGrid},nlat_half::Integer) = npoints_gaussian(nlat_half)
get_npoints(::Type{<:FullClenshawGrid},nlat_half::Integer) = npoints_clenshaw(nlat_half)
get_npoints(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = npoints_octahedral(nlat_half,false)
get_npoints(::Type{<:OctahedralClenshawGrid},nlat_half::Integer) = npoints_octahedral(nlat_half,true)
get_npoints(::Type{<:AbstractHEALPixGrid},nside::Integer) = npoints_healpix(nside)

# colatitude [radians] vectors
get_colat(::Type{<:FullClenshawGrid},nlat_half::Integer) = [j/(2nlat_half)*π for j in 1:2nlat_half-1]
get_colat(::Type{<:FullGaussianGrid},nlat_half::Integer) =
            π .- acos.(FastGaussQuadrature.gausslegendre(2nlat_half)[1])
get_colat(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = get_colat(FullGaussianGrid,nlat_half)
get_colat(::Type{<:OctahedralClenshawGrid},nlat_half::Integer) = get_colat(FullClenshawGrid,nlat_half)
get_colat(G::Type{<:HEALPixGrid},nside::Integer) =
            [acos(Healpix.ring2z(Healpix.Resolution(nside),j)) for j in 1:nlat_healpix(nside)]

# lon [radians] vectors for full grids (empty vectors otherwise)
get_lon(::Type{<:AbstractFullGrid},nlat_half::Integer) = 
            collect(range(0,2π,step=2π/get_nlon(FullClenshawGrid,nlat_half))[1:end-1])
get_lon(::Type{<:AbstractOctahedralGrid},nlat_half::Integer) = Float64[]
get_lon(::Type{<:AbstractHEALPixGrid},nside::Integer) = Float64[]

# get coordinates
function get_colatlons(G::Type{<:AbstractFullGrid},nlat_half::Integer)

    colat = get_colat(G,nlat_half)
    lon = get_lon(G,nlat_half)
    nlon = get_nlon(G,nlat_half)

    colats = zeros(get_npoints(G,nlat_half))
    lons = zeros(get_npoints(G,nlat_half))

    for j in 1:2nlat_half-nlat_odd(G)
        for i in 1:nlon
            ij = i + (j-1)*nlon
            colats[ij] = colat[j]
            lons[ij] = lon[i]
        end
    end

    return colats,lons
end

function get_colatlons(G::Type{<:AbstractOctahedralGrid},nlat_half::Integer)
    
    npoints = get_npoints(G,nlat_half)
    colat = get_colat(G,nlat_half)

    colats = zeros(npoints)
    lons = zeros(npoints)

    j = 1
    for i in 1:2nlat_half-nlat_odd(G)
        nlon = get_nlon_per_ring(G,nlat_half,i)
        lon = collect(0:2π/nlon:2π-π/nlon)

        colats[j:j+nlon-1] .= colat[i]
        lons[j:j+nlon-1] .= lon

        j += nlon
    end

    return colats, lons
end

function get_colatlons(G::Type{<:AbstractHEALPixGrid},nside::Integer)
    npoints = get_npoints(G,nside)
    colats_lons = [Healpix.pix2angRing(Healpix.Resolution(nside),i) for i in 1:npoints]
    colats = [colat_lon[1] for colat_lon in colats_lons]
    lons = [colat_lon[2] for colat_lon in colats_lons]
    return colats, lons
end

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

function each_index_in_ring(G::Type{<:AbstractFullGrid},    # function for full grids
                            j::Integer,                     # ring index north to south
                            nlat_half::Integer)             # resolution param

    @boundscheck 0 < j <= (2nlat_half-nlat_odd(G)) || throw(BoundsError)    # valid ring index?
    nlon = 4nlat_half                                       # number of longitudes per ring (const)
    index_1st = (j-1)*nlon + 1                              # first in-ring index i
    index_end = j*nlon                                      # last in-ring index i  
    return index_1st:index_end                              # range of js in ring
end

function each_index_in_ring(G::Type{<:AbstractOctahedralGrid},   # function for octahedral grids
                            j::Integer,                     # ring index north to south
                            nlat_half::Integer)             # resolution param

    @boundscheck 0 < j <= (2nlat_half-nlat_odd(G)) || throw(BoundsError)  # ring index valid?
    if j <= nlat_half                                       # northern hemisphere incl Equator
        index_1st = 2j*(j+7) - 15                           # first in-ring index i
        index_end = 2j*(j+9)                                # last in-ring index i
    else                                                    # southern hemisphere excl Equator
        j = 2nlat_half-nlat_odd(G) - j + 1                  # mirror ring index around Equator
        n = npoints_octahedral(nlat_half,nlat_odd(G))+1     # number of grid points + 1
        index_1st = n - 2j*(j+9)                            # count backwards
        index_end = n - (2j*(j+7) - 15)
    end
    return index_1st:index_end                              # range of i's in ring
end


function each_index_in_ring(::Type{<:AbstractHEALPixGrid},  # function for HEALPix grids
                            j::Integer,                     # ring index north to south
                            nside::Integer)                 # resolution param

    @boundscheck 0 < j < 4nside || throw(BoundsError)       # ring index valid?
    if j < nside                                            # northern polar cap
        index_1st = 2j*(j-1) + 1                            # first in-ring index i
        index_end = 2j*(j+1)                                # last in-ring index i
    elseif j <= 3nside                                      # equatorial zone with const nlon
        n = 2nside^2-2nside                                 # points in polar cap
        nlon = 4nside                                       # points on latitude rings
        j = j - nside + 1                                   # offset ring index into eq zone
        index_1st = n + (j-1)*nlon + 1                      # add constant nlon per ring
        index_end = n + j*nlon
    else                                                    # southern polar cap
        n = 12nside^2                                       # total number of points
        j = 4nside - j                                      # count ring index from south pole
        index_1st = n - 2j*(j+1) + 1                        # count backwards
        index_end = n - 2j*(j-1)
    end
    return index_1st:index_end                              # range of i's in ring
end

function each_index_in_ring(grid::G,j::Integer) where {G<:AbstractGrid}
    return each_index_in_ring(G,j,get_nresolution(grid))
end

eachgridpoint(grid::G) where {G<:AbstractGrid} = Base.OneTo(get_npoints(G,get_nresolution(grid)))

# QUADRATURE WEIGHTS
# gaussian_weights are exact for Gaussian latitudes when nlat > (2T+1)/2
# clenshaw_curtis_weights are exact for equi-angle latitudes when nlat > 2T+1
# riemann_weights not exact but used for HEALPix
gaussian_weights(nlat_half::Integer) = FastGaussQuadrature.gausslegendre(2nlat_half)[2][1:nlat_half]

function clenshaw_curtis_weights(nlat_half::Integer)
    nlat = 2nlat_half - 1
    θs = get_colat(FullClenshawGrid,nlat_half)
    return [4sin(θj)/(nlat+1)*sum([sin(p*θj)/p for p in 1:2:nlat]) for θj in θs[1:nlat_half]]
end

get_quadrature_weights(::Type{<:FullClenshawGrid},nlat_half::Integer) = clenshaw_curtis_weights(nlat_half)
get_quadrature_weights(::Type{<:FullGaussianGrid},nlat_half::Integer) = gaussian_weights(nlat_half)
get_quadrature_weights(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = gaussian_weights(nlat_half)
get_quadrature_weights(::Type{<:OctahedralClenshawGrid},nlat_half::Integer) = clenshaw_curtis_weights(nlat_half)

# HEALPix's solid angle is always 4π/npoints
get_quadrature_weights(G::Type{<:HEALPixGrid},nside::Integer) = 4π/get_npoints(G,nside)*ones(get_nlat_half(G,nside))

