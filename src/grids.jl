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

A FullClenshawGrid is a regular latitude-longitude grid with `nlat` equi-spaced latitudes,
and the same `nlon` longitudes for every latitude ring. The grid points are closer in zonal direction
around the poles. The values of all grid points are stored in a vector field `v` that unravels
the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct FullClenshawGrid{T} <: AbstractFullGrid{T}
    v::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere

    FullClenshawGrid{T}(v,nlat_half) where T = length(v) == 8nlat_half^2 ?
    new(v,nlat_half) : error("$(length(v))-element Vector{$(eltype(v))} cannot be used to create a "*
        "L$nlat_half ($(4nlat_half)x$(2nlat_half)) FullClenshawGrid{$T}.")
end

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
    OctahedralGaussianGrid{T}(v,nlat_half) where T = length(v) == npoints_octahedral(nlat_half) ?
    new(v,nlat_half) : error("$(length(v))-element Vector{$(eltype(v))}"*
    "cannot be used to create a O$(nlat_half) OctahedralGaussianGrid{$T}.")
end

# number of points and longitudes per ring on the octahedral grid
npoints_octahedral(nlat_half::Integer) = 4nlat_half^2 + 36nlat_half
nlon_octahedral(ilat::Integer) = 16+4ilat

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
    @assert is_power_2(nside) "HEALPixGrid: nside=$nside is not a power of 2."
    return true
end

npoints_healpix(nside::Integer) = nside_assert(nside) ? 12nside^2 : nothing
nlat_healpix(nside::Integer) = nside_assert(nside) ? 4nside-1 : nothing
nlon_healpix(nside::Integer,ilat::Integer) = nside_assert(nside) ? min(4ilat,4nside) : nothing
nlon_healpix(nside::Integer) = nside_assert(nside) ? 4nside : nothing

# define for all grids that the type T can be infered from the elements in data vector
# whether the resolution parameter n is provided or not (hence the ...)
FullClenshawGrid(v::AbstractVector,n::Integer...) = FullClenshawGrid{eltype(v)}(v,n...)
FullGaussianGrid(v::AbstractVector,n::Integer...) = FullGaussianGrid{eltype(v)}(v,n...)
OctahedralGaussianGrid(v::AbstractVector,n::Integer...) = OctahedralGaussianGrid{eltype(v)}(v,n...)
HEALPixGrid(v::AbstractVector,n::Integer...) = HEALPixGrid{eltype(v)}(v,n...)

# infer resolution parameter nlat_half or nside from length of vector
FullClenshawGrid{T}(v::AbstractVector) where T = FullClenshawGrid(v,round(Int,sqrt(length(v)/8)))
FullGaussianGrid{T}(v::AbstractVector) where T = FullGaussianGrid(v,round(Int,sqrt(length(v)/8)))
HEALPixGrid{T}(v::AbstractVector) where T = HEALPixGrid(v,round(Int,sqrt(length(v)/12)))

# for octahedral define the inverse of npoints_octahedral first
nlat_half_octahedral(npoints::Integer) = round(Int,-9/2+sqrt((9/2)^2 + npoints/4))
OctahedralGaussianGrid{T}(v::AbstractVector) where T = OctahedralGaussianGrid(v,nlat_half_octahedral(length(v)))

# convert an AbstractMatrix to the full grids
FullClenshawGrid(M::AbstractMatrix{T}) where T = FullClenshawGrid{T}(vec(M))
FullGaussianGrid(M::AbstractMatrix{T}) where T = FullGaussianGrid{T}(vec(M))

# and vice versa
Base.Matrix(G::AbstractFullGrid{T}) where T = Matrix{T}(reshape(G.v,:,2G.nlat_half))

# generator functions for grid
Base.zeros(::Type{FullClenshawGrid{T}},nlat_half::Integer) where T = 
                FullClenshawGrid(zeros(T,8nlat_half^2),nlat_half)
Base.zeros(::Type{FullGaussianGrid{T}},nlat_half::Integer) where T =
                FullGaussianGrid(zeros(T,8nlat_half^2),nlat_half)
Base.zeros(::Type{OctahedralGaussianGrid{T}},nlat_half::Integer) where T =
                OctahedralGaussianGrid(zeros(T,npoints_octahedral(nlat_half)),nlat_half)
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
truncation_order(::Type{<:HEALPixGrid}) = 1                 # linear (in longitude)

get_truncation(::Type{<:FullClenshawGrid},nlat_half::Integer) = floor(Int,(4nlat_half-1)/4)
get_truncation(::Type{<:FullGaussianGrid},nlat_half::Integer) = floor(Int,(4nlat_half-1)/3)
get_truncation(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = nlat_half-1
get_truncation(::Type{<:HEALPixGrid},nside::Integer) = nside_assert(nside) ? 2nside-1 : nothing
get_truncation(grid::G) where {G<:AbstractGrid} = get_truncation(G,get_nresolution(grid))

get_resolution(::Type{<:FullClenshawGrid},trunc::Integer) = roundup_fft(ceil(Int,(4*trunc+1)/4))
get_resolution(::Type{<:FullGaussianGrid},trunc::Integer) = roundup_fft(ceil(Int,(3*trunc+1)/4))
get_resolution(::Type{<:OctahedralGaussianGrid},trunc::Integer) = roundup_fft(trunc+1)
get_resolution(::Type{<:HEALPixGrid},trunc::Integer) = roundup_fft(ceil(Int,(trunc+1)/2),small_primes=[2])

# common interface for the resolution parameter nlat_half/nside
get_nresolution(G::AbstractFullGrid) = G.nlat_half
get_nresolution(G::AbstractOctahedralGrid) = G.nlat_half
get_nresolution(G::AbstractHEALPixGrid) = G.nside

# define nlat_half for all grids (HEALPixGrid is different as it doesn't use nlat_half as resolution parameter)
get_nlat_half(::Type{<:AbstractGrid},nlat_half::Integer) = nlat_half
get_nlat_half(::Type{<:AbstractHEALPixGrid},nside::Integer) = (nlat_healpix(nside)+1)÷2

# define whether there's an odd number of latitude rings for a grid
nlat_odd(::Type{<:AbstractFullGrid}) = false
nlat_odd(::Type{<:AbstractOctahedralGrid}) = false
nlat_odd(::Type{<:AbstractHEALPixGrid}) = true

# return the maxmimum number of longitude points for a grid and its resolution parameter nlat_half/nside
get_nlon(::Type{<:AbstractFullGrid},nlat_half::Integer) = 4nlat_half
get_nlon(::Type{<:AbstractOctahedralGrid},nlat_half::Integer) = nlon_octahedral(nlat_half)
get_nlon(::Type{<:AbstractHEALPixGrid},nside::Integer) = nside_assert(nside) ? nlon_healpix(nside) : nothing

# get the number of longitude points at given latitude ring i number 1...J from north to south
get_nlon_per_ring(::Type{<:AbstractFullGrid},nlat_half::Integer,iring::Integer) = 4nlat_half
function get_nlon_per_ring(::Type{<:AbstractOctahedralGrid},nlat_half::Integer,iring::Integer)
    @assert 0 < iring <= 2nlat_half "Ring $iring is outside O$nlat_half grid."
    iring = iring > nlat_half ? 2nlat_half - iring + 1 : iring      # flip north south due to symmetry
    return nlon_octahedral(iring)
end
function get_nlon_per_ring(::Type{<:AbstractHEALPixGrid},nside::Integer,iring::Integer)
    nlat = nlat_healpix(nside)
    @assert 0 < iring <= nlat "Ring $iring is outside H$nside grid."
    nlat_half = (nlat+1)÷2
    iring = iring > nlat_half ? nlat - iring + 1 : iring      # flip north south due to symmetry
    return nlon_healpix(nside,iring)
end

# total number of grid points per grid type
get_npoints(::Type{<:AbstractFullGrid},nlat_half::Integer) = 8nlat_half^2
get_npoints(::Type{<:AbstractOctahedralGrid},nlat_half::Integer) = npoints_octahedral(nlat_half)
get_npoints(::Type{<:AbstractHEALPixGrid},nside::Integer) = nside_assert(nside) ? npoints_healpix(nside) : nothing

# colatitude [radians] vectors
get_colat(::Type{<:FullClenshawGrid},nlat_half::Integer) = [j/(2nlat_half+1)*π for j in 1:2nlat_half]
get_colat(::Type{<:FullGaussianGrid},nlat_half::Integer) =
            π .- acos.(FastGaussQuadrature.gausslegendre(2nlat_half)[1])
get_colat(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = get_colat(FullGaussianGrid,nlat_half)
get_colat(::Type{<:HEALPixGrid},nside::Integer) =
            [acos(Healpix.ring2z(Healpix.Resolution(nside),i)) for i in 1:nlat_healpix(nside)]

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

    for j in 1:2nlat_half
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
    for i in 1:2nlat_half
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

# INDEXING HELPERS
function get_last_index_per_ring(G::Type{<:AbstractGrid},nresolution::Integer)
    nlat_half = get_nlat_half(G,nresolution)    # contains equator for HEALPix
    nlat = 2nlat_half - nlat_odd(G)             # one less if grids have odd # of latitude rings
    nlons = [get_nlon_per_ring(G,nresolution,i) for i in 1:nlat]
    last_indices = cumsum(nlons)                # last index is always sum of all previous points
    return last_indices
end

function get_first_index_per_ring(G::Type{<:AbstractGrid},nresolution::Integer)
    last_indices = get_last_index_per_ring(G,nresolution)
    first_indices = zero(last_indices)
    first_indices[1] = 1
    for i in 1:length(first_indices)-1
        first_indices[i+1] = last_indices[i]+1
    end
    return first_indices
end

function eachring(grid::G) where {G<:AbstractGrid}
    nlat_half = get_nlat_half(G,get_nresolution(grid))  # contains equator for HEALPix
    nlat = 2nlat_half - nlat_odd(G)                     # -1 for odd # of latitude rings
    return Base.OneTo(nlat)                             # return iterable range
end

function each_index_in_ring(::Type{<:AbstractFullGrid},     # function for full grids
                                    i::Integer,             # ring index north to south
                                    nlat_half::Integer)     # resolution param

    @boundscheck 0 < i <= 2nlat_half || throw(BoundsError)  # valid ring index?
    nlon = 4nlat_half                                       # number of longitudes per ring (const)
    index_1st = (i-1)*nlon + 1                              # first in-ring index j
    index_end = i*nlon                                      # last in-ring index j      
    return index_1st:index_end                              # range of js in ring
end

function each_index_in_ring(::Type{<:AbstractOctahedralGrid},   # function for octahedral
                                    i::Integer,                 # ring index north to south
                                    nlat_half::Integer)         # resolution param

    @boundscheck 0 < i <= 2nlat_half || throw(BoundsError)  # ring index valid?
    if i <= nlat_half                                       # northern hemisphere
        index_1st = 2i*(i+7) - 15                           # last in-ring index j      
        index_end = 2i*(i+9)                                # last in-ring index j
    else                                                    # southern hemisphere
        i = 2nlat_half - i + 1                              # mirror ring index around Equator
        n = npoints_octahedral(nlat_half)+1                 # number of grid points + 1
        index_1st = n - 2i*(i+9)                            # count backwards
        index_end = n - (2i*(i+7) - 15)
    end
    return index_1st:index_end                              # range of js in ring
end

function each_index_in_ring(::Type{<:AbstractHEALPixGrid},  # function for HEALPix grids
                                    i::Integer,             # ring index north to south
                                    nside::Integer)         # resolution param

    @boundscheck 0 < i < 4nside || throw(BoundsError)       # ring index valid?
    if i < nside                                            # northern polar cap
        index_1st = 2i*(i-1) + 1                            # first in-ring index j
        index_end = 2i*(i+1)                                # last in-ring index j
    elseif i <= 3nside                                      # equatorial zone with const nlon
        n = 2nside^2-2nside                                 # points in polar cap
        nlon = 4nside                                       # points on latitude rings
        i = i - nside + 1                                   # offset ring index into eq zone
        index_1st = n + (i-1)*nlon + 1                      # add constant nlon per ring
        index_end = n + i*nlon
    else                                                    # southern polar cap
        n = 12nside^2                                       # total number of points
        i = 4nside - i                                      # count ring index from south pole
        index_1st = n - 2i*(i+1) + 1                        # count backwards
        index_end = n - 2i*(i-1)
    end
    return index_1st:index_end                              # range of js in ring
end

@inline function each_index_in_ring(grid::G,i::Integer) where {G<:AbstractGrid}
    return each_index_in_ring(G,i,get_nresolution(grid))
end

eachgridpoint(grid::G) where {G<:AbstractGrid} = Base.OneTo(get_npoints(G,get_nresolution(grid)))

# QUADRATURE WEIGHTS
# gaussian_weights are exact for Gaussian latitudes when nlat > (2T+1)/2
# clenshaw_curtis_weights are exact for equi-angle latitudes when nlat > 2T+1
# riemann_weights not exact but used for HEALPix
gaussian_weights(nlat_half::Integer) = FastGaussQuadrature.gausslegendre(2nlat_half)[2][1:nlat_half]

function clenshaw_curtis_weights(nlat_half::Integer)
    nlat = 2nlat_half
    θs = get_colat(FullClenshawGrid,nlat_half)
    return [4sin(θj)/(nlat+1)*sum([sin(p*θj)/p for p in 1:2:nlat]) for θj in θs[1:nlat_half]]
end

get_quadrature_weights(::Type{<:FullClenshawGrid},nlat_half::Integer) = clenshaw_curtis_weights(nlat_half)
get_quadrature_weights(::Type{<:FullGaussianGrid},nlat_half::Integer) = gaussian_weights(nlat_half)
get_quadrature_weights(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = gaussian_weights(nlat_half)
get_quadrature_weights(::Type{<:HEALPixGrid},nside::Integer) =
                nside_assert(nside) ? 2/12nside^2*[min(4iring,4nside) for iring in 1:2nside+1] : nothing