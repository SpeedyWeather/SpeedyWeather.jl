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
    @assert is_power_2(nside) "nside=$nside is not a power of two."
    return true
end

npoints_healpix(nside::Integer) = nside_assert(nside) ? 12nside^2 : nothing
nlat_healpix(nside::Integer) = nside_assert(nside) ? 4nside-1 : nothing
nlon_healpix(nside::Integer,ilat::Integer) = nside_assert(nside) ? min(4ilat,4nside) : nothing
nlon_healpix(nside::Integer) = nside_assert(nside) ? 4nside : nothing

# define for all grids that the type T can be infered from the elements in data vector
for grid in [Symbol(grid) for grid in subtypes(AbstractGrid)]
    @eval begin
        # n is the resolution parameter, either nlat_half for Gaussian/LatLon grids, or nside for HEALPix
        $grid(v::AbstractVector{T},n::Integer) where T = $grid{T}(v,n)
        $grid(v::AbstractVector{T}) where T = $grid{T}(v)
    end
end

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
Base.zeros(::Type{G},n::Integer) where {G<:AbstractGrid} = zero(G{Float64},n)

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