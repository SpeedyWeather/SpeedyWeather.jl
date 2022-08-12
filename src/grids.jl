abstract type AbstractGrid{T} <: AbstractVector{T} end

# all Abstract grids have their grid points stored in a vector field `v`
# propagate length, size, getindex, setindex! for that
Base.length(G::AbstractGrid) = length(G.v)
Base.size(G::AbstractGrid) = size(G.v)
@inline function Base.getindex(G::AbstractGrid,k::Integer)
    @boundscheck 0 < k <= length(G.v) || throw(BoundsError(G,k))
    @inbounds r = G.v[k]
    return r
end

@inline Base.setindex!(G::AbstractGrid,x,k::Integer) = setindex!(G.v,x,k)

"""
    G = FullLatLonGrid{T}

A FullLatLonGrid is a regular latitude-longitude grid that is de

`nlat` equi-spaced latitudes,
and the same `nlon` longitudes for every latitude ring. The grid points are closer in zonal direction
around the poles. The values of all grid points are stored in a vector field `v` that unravels
the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct FullLatLonGrid{T} <: AbstractGrid{T}
    v::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere

    FullLatLonGrid{T}(v,nlat_half) where T = length(v) == 8nlat_half^2 ?
    new(v,nlat_half) : error("$(length(v))-element Vector{$(eltype(v))} cannot be used to create a "*
        "L$nlat_half ($(4nlat_half)x$(2nlat_half)) FullLatLonGrid{$T}.")
end

"""
    G = FullGaussianGrid{T}

A full Gaussian grid is a regular latitude-longitude grid that uses `nlat` Gaussian latitudes,
and the same `nlon` longitudes for every latitude ring. The grid points are closer in zonal direction
around the poles. The values of all grid points are stored in a vector field `v` that unravels
the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct FullGaussianGrid{T} <: AbstractGrid{T}
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
struct OctahedralGaussianGrid{T} <: AbstractGrid{T}
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
struct HEALPixGrid{T} <: AbstractGrid{T}
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
FullLatLonGrid(v::AbstractVector,n::Integer...) = FullLatLonGrid{eltype(v)}(v,n...)
FullGaussianGrid(v::AbstractVector,n::Integer...) = FullGaussianGrid{eltype(v)}(v,n...)
OctahedralGaussianGrid(v::AbstractVector,n::Integer...) = OctahedralGaussianGrid{eltype(v)}(v,n...)
HEALPixGrid(v::AbstractVector,n::Integer...) = HEALPixGrid{eltype(v)}(v,n...)

# infer resolution parameter nlat_half or nside from length of vector
FullLatLonGrid{T}(v::AbstractVector) where T = FullLatLonGrid(v,round(Int,sqrt(length(v)/8)))
FullGaussianGrid{T}(v::AbstractVector) where T = FullGaussianGrid(v,round(Int,sqrt(length(v)/8)))
HEALPixGrid{T}(v::AbstractVector) where T = HEALPixGrid(v,round(Int,sqrt(length(v)/12)))

# for octahedral define the inverse of npoints_octahedral first
nlat_half_octahedral(npoints::Integer) = round(Int,-9/2+sqrt((9/2)^2 + npoints/4))
OctahedralGaussianGrid{T}(v::AbstractVector) where T = OctahedralGaussianGrid(v,nlat_half_octahedral(length(v)))

# MATCHING SPECTRAL TO GRID POINT RESOLUTION
get_truncation(::Type{FullLatLonGrid},nlat_half::Integer) = floor(Int,(4nlat_half-1)/4)
get_truncation(::Type{FullGaussianGrid},nlat_half::Integer) = floor(Int,(4nlat_half-1)/3)
get_truncation(::Type{OctahedralGaussianGrid},nlat_half::Integer) = nlat_half-1
get_truncation(::Type{HEALPixGrid},nside::Integer) = nside_assert(nside) ? 2nside-1 : nothing

get_resolution(::Type{FullLatLonGrid},trunc::Integer) = roundup_fft(ceil(Int,(4*trunc+1)/4))
get_resolution(::Type{FullGaussianGrid},trunc::Integer) = roundup_fft(ceil(Int,(3*trunc+1)/4))
get_resolution(::Type{OctahedralGaussianGrid},trunc::Integer) = roundup_fft(trunc+1)
get_resolution(::Type{HEALPixGrid},trunc::Integer) = roundup_fft(ceil(Int,(trunc+1)/2),small_primes=[2])

# common interface for the resolution parameter nlat_half/nside
get_nresolution(G::FullLatLonGrid) = G.nlat_half
get_nresolution(G::FullGaussianGrid) = G.nlat_half
get_nresolution(G::OctahedralGaussianGrid) = G.nlat_half
get_nresolution(G::HEALPixGrid) = G.nside

# define nlat_half for all grids (HEALPixGrid is different as it doesn't use nlat_half as resolution parameter)
get_nlat_half(::Type{G},nlat_half::Integer) where {G<:AbstractGrid} = nlat_half
get_nlat_half(::Type{HEALPixGrid},nside::Integer) = (nlat_healpix(nside)+1)÷2

# define whether there's an odd number of latitude rings for a grid
nlat_odd(::Type{G}) where {G<:AbstractGrid} = false
nlat_odd(::Type{HEALPixGrid}) = true

# return the maxmimum number of longitude points for a grid and its resolution parameter nlat_half/nside
get_nlon(::Type{FullLatLonGrid},nlat_half::Integer) = 4nlat_half
get_nlon(::Type{FullGaussianGrid},nlat_half::Integer) = 4nlat_half
get_nlon(::Type{OctahedralGaussianGrid},nlat_half::Integer) = nlon_octahedral(nlat_half)
get_nlon(::Type{HEALPixGrid},nside::Integer) = nside_assert(nside) ? nlon_healpix(nside) : nothing

get_nlon_per_ring(::Type{FullLatLonGrid},nlat_half::Integer,iring::Integer) = 4nlat_half
get_nlon_per_ring(::Type{FullGaussianGrid},nlat_half::Integer,iring::Integer) = 4nlat_half
function get_nlon_per_ring(::Type{OctahedralGaussianGrid},nlat_half::Integer,iring::Integer)
    @assert 0 < iring <= 2nlat_half "Ring $iring is outside O$nlat_half grid."
    iring = iring > nlat_half ? 2nlat_half - iring : iring      # flip north south due to symmetry
    return nlon_octahedral(iring)
end
function get_nlon_per_ring(::Type{HEALPixGrid},nside::Integer,iring::Integer)
    nlat = nlat_healpix(nside)
    @assert 0 < iring <= nlat "Ring $iring is outside H$nside grid."
    nlat_half = (nlat+1)÷2
    iring = iring > nlat_half ? nlat - iring : iring      # flip north south due to symmetry
    return nlon_healpix(iring)
end

# total number of grid points per grid type
get_npoints(::Type{FullLatLonGrid},nlat_half::Integer) = 8nlat_half^2
get_npoints(::Type{FullGaussianGrid},nlat_half::Integer) = 8nlat_half^2
get_npoints(::Type{OctahedralGaussianGrid},nlat_half::Integer) = npoints_octahedral(nlat_half)
get_npoints(::Type{HEALPixGrid},nside::Integer) = nside_assert(nside) ? npoints_healpix(nside) : nothing

# colatitude [radians] vectors
get_colat(::Type{FullLatLonGrid},nlat_half::Integer) = [j/(2nlat_half+1)*π for j in 1:nlat_half]
get_colat(::Type{FullGaussianGrid},nlat_half::Integer) =
            π .- acos.(FastGaussQuadrature.gausslegendre(2nlat_half)[1])
get_colat(::Type{OctahedralGaussianGrid},nlat_half::Integer) = get_colat(OctahedralGaussianGrid,nlat_half)
get_colat(::Type{HEALPixGrid},nside::Integer) =
            [acos(Healpix.ring2z(Healpix.Resolution(nside),i)) for i in nlat_healpix(nside)]

# lon [radians] vectors for full grids (empty vectors otherwise)
get_lon(::Type{FullLatLonGrid},nlat_half::Integer) = 
            collect(range(0,2π,step=2π/get_nlon(FullLatLonGrid,nlat_half))[1:end-1])
get_lon(::Type{FullGaussianGrid},nlat_half::Integer) = get_lon(FullLatLonGrid,nlat_half)
get_lon(::Type{OctahedralGaussianGrid},nlat_half::Integer) = Float64[]
get_lon(::Type{HEALPixGrid},nside::Integer) = Float64[]

# get coordinates
function get_colatlons(::Type{G},nlat_half::Integer) where {G<:Union{FullLatLonGrid,FullGaussianGrid}}

    colat = get_colat(G,nlat_half)
    lon = get_lon(G,nlat_half)

    colats = zeros(get_npoints(G,nlat_half))
    lons = zeros(get_npoints(G,nlat_half))

    for j in 1:2nlat_half
        for i in 1:get_nlon(G,nlat_half)
            ij = i + (j-1)*2nlat_half
            colats[ij] = colat[j]
            lons[ij] = lon[i]
        end
    end

    return colats,lons
end

# function get_colatlons(::Type{OctahedralGaussianGrid},nlat_half::Integer)

# QUADRATURE WEIGHTS
# gaussian_weights are exact for Gaussian latitudes when nlat > (2T+1)/2
# clenshaw_curtis_weights are exact for equi-angle latitudes when nlat > 2T+1
# riemann_weights not exact but used for HEALPix
gaussian_weights(nlat_half::Integer) = FastGaussQuadrature.gausslegendre(2nlat_half)[2][1:nlat_half]

function clenshaw_curtis_weights(nlat_half::Integer)
    nlat = 2nlat_half
    θs = get_colat(FullLatLonGrid,nlat_half)
    return [4sin(θj)/(nlat+1)*sum([sin(p*θj)/p for p in 1:2:nlat]) for θj in θs[1:nlat_half]]
end

get_quadrature_weights(::Type{FullLatLonGrid},nlat_half::Integer) = clenshaw_curtis_weights(nlat_half)
get_quadrature_weights(::Type{FullGaussianGrid},nlat_half::Integer) = gaussian_weights(nlat_half)
get_quadrature_weights(::Type{OctahedralGaussianGrid},nlat_half::Integer) = gaussian_weights(nlat_half)
get_quadrature_weights(::Type{HEALPixGrid},nside::Integer) =
                nside_assert(nside) ? 2/12nside^2*[min(4iring,4nside) for iring in 1:2nside+1] : nothing

# generator functions for grid
Base.zeros(::Type{FullLatLonGrid{T}},nlat_half::Integer) where T = 
                FullLatLonGrid(zeros(T,8nlat_half^2),nlat_half)
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

"""
    m = roundup_fft(n::Int;
                    small_primes::Vector{Int}=[2,3,5])

Returns an integer `m >= n` with only small prime factors 2, 3, 5 (default, others can be specified
with the keyword argument `small_primes`) to obtain an efficiently fourier-transformable number of
longitudes, m = 2^i * 3^j * 5^k >= n, with i,j,k >=0.
"""
function roundup_fft(n::Integer;small_primes::Vector{T}=[2,3,5]) where {T<:Integer}
    factors_not_in_small_primes = true      # starting condition for while loop
    n += isodd(n) ? 1 : 0                   # start with an even n
    while factors_not_in_small_primes
        
        factors = Primes.factor(n)          # prime factorization
        all_factors_small = true            # starting condition
        
        for i in 1:length(factors)          # loop over factors and check they are small
            factor = factors.pe[i].first    # extract factor from factors
            all_factors_small &= factor in small_primes
        end
        
        factors_not_in_small_primes = ~all_factors_small    # all factors small will abort while loop
        n += 2                                              # test for next larger even n
    
    end
    return n-2      # subtract unnecessary last += 2 addition
end