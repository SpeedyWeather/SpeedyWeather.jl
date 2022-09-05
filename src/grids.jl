"""
    abstract type AbstractFullGrid{T} <: AbstractGrid{T} end

An `AbstractFullGrid` is a horizontal grid with a constant number of longitude
points across latitude rings. Different latitudes can be used, Gaussian latitudes,
equi-angle latitdes, or others."""
abstract type AbstractFullGrid{T} <: AbstractGrid{T} end

get_nresolution(grid::AbstractFullGrid) = grid.nlat_half
get_nlon(::Type{<:AbstractFullGrid},nlat_half::Integer) = 4nlat_half
get_nlon_per_ring(::Type{<:AbstractFullGrid},nlat_half::Integer,j::Integer) = 4nlat_half
get_lon(Grid::Type{<:AbstractFullGrid},nlat_half::Integer) = 
            collect(range(0,2π,step=2π/get_nlon(Grid,nlat_half))[1:end-1])

function get_colatlons(Grid::Type{<:AbstractFullGrid},nlat_half::Integer)

    colat = get_colat(Grid,nlat_half)
    lon = get_lon(Grid,nlat_half)
    nlon = get_nlon(Grid,nlat_half)

    colats = zeros(get_npoints(Grid,nlat_half))    # preallocate
    lons = zeros(get_npoints(Grid,nlat_half))

    for j in 1:2nlat_half-nlat_odd(Grid)           # populate preallocated colats,lons
        for i in 1:nlon
            ij = i + (j-1)*nlon                    # continuous index ij
            colats[ij] = colat[j]
            lons[ij] = lon[i]
        end
    end

    return colats,lons
end

function each_index_in_ring(Grid::Type{<:AbstractFullGrid},     # function for full grids
                            j::Integer,                         # ring index north to south
                            nlat_half::Integer)                 # resolution param

    @boundscheck 0 < j <= (2nlat_half-nlat_odd(Grid)) || throw(BoundsError)    # valid ring index?
    nlon = 4nlat_half                                       # number of longitudes per ring (const)
    index_1st = (j-1)*nlon + 1                              # first in-ring index i
    index_end = j*nlon                                      # last in-ring index i  
    return index_1st:index_end                              # range of js in ring
end

"""
    abstract type AbstractOctahedralGrid{T} <: AbstractGrid{T} end

An `AbstractOctahedralGrid` is a horizontal grid with 16+4i longitude
points on the latitude ring i starting with i=1 around the pole.
Different latitudes can be used, Gaussian latitudes, equi-angle latitdes, or others."""
abstract type AbstractOctahedralGrid{T} <: AbstractGrid{T} end

get_nresolution(grid::AbstractOctahedralGrid) = grid.nlat_half
get_nlon(::Type{<:AbstractOctahedralGrid},nlat_half::Integer) = nlon_octahedral(nlat_half)

function get_nlon_per_ring(Grid::Type{<:AbstractOctahedralGrid},nlat_half::Integer,j::Integer)
    nlat = 2nlat_half-nlat_odd(Grid)
    @assert 0 < j <= nlat "Ring $j is outside O$nlat_half grid."
    j = j > nlat_half ? nlat - j + 1 : j      # flip north south due to symmetry
    return nlon_octahedral(j)
end

get_lon(::Type{<:AbstractOctahedralGrid},nlat_half::Integer) = Float64[]    # only defined for full grids

function get_colatlons(Grid::Type{<:AbstractOctahedralGrid},nlat_half::Integer)
    
    npoints = get_npoints(Grid,nlat_half)
    colat = get_colat(Grid,nlat_half)

    colats = zeros(npoints)                 # preallocate arrays
    lons = zeros(npoints)

    ij = 1                                  # continuous index
    for j in 1:2nlat_half-nlat_odd(Grid)    # populate arrays ring by ring
        nlon = get_nlon_per_ring(Grid,nlat_half,j)
        lon = collect(0:2π/nlon:2π-π/nlon)

        colats[ij:ij+nlon-1] .= colat[j]
        lons[ij:ij+nlon-1] .= lon

        ij += nlon
    end

    return colats, lons
end

function each_index_in_ring(Grid::Type{<:AbstractOctahedralGrid},
                            j::Integer,                     # ring index north to south
                            nlat_half::Integer)             # resolution param

    @boundscheck 0 < j <= (2nlat_half-nlat_odd(Grid)) || throw(BoundsError)  # ring index valid?
    if j <= nlat_half                                       # northern hemisphere incl Equator
        index_1st = 2j*(j+7) - 15                           # first in-ring index i
        index_end = 2j*(j+9)                                # last in-ring index i
    else                                                    # southern hemisphere excl Equator
        j = 2nlat_half-nlat_odd(Grid) - j + 1               # mirror ring index around Equator
        n = npoints_octahedral(nlat_half,nlat_odd(Grid))+1  # number of grid points + 1
        index_1st = n - 2j*(j+9)                            # count backwards
        index_end = n - (2j*(j+7) - 15)
    end
    return index_1st:index_end                              # range of i's in ring
end

"""
    abstract type AbstractHEALPixGrid{T} <: AbstractGrid{T} end

An `AbstractHEALPixGrid` is a horizontal grid similar to the standard HEALPixGrid,
but different latitudes can be used, the default HEALPix latitudes or others."""
abstract type AbstractHEALPixGrid{T} <: AbstractGrid{T} end

get_nresolution(G::AbstractHEALPixGrid) = G.nside
get_nlat_half(::Type{<:AbstractHEALPixGrid},nside::Integer) = 2nside
nlat_odd(::Type{<:AbstractHEALPixGrid}) = true
get_nlon(::Type{<:AbstractHEALPixGrid},nside::Integer) = nside_assert(nside) ? nlon_healpix(nside) : nothing

function get_nlon_per_ring(::Type{<:AbstractHEALPixGrid},nside::Integer,j::Integer)
    nlat = nlat_healpix(nside)
    @assert 0 < j <= nlat "Ring $j is outside H$nside grid."
    nlat_half = (nlat+1)÷2
    j = j > nlat_half ? nlat - j + 1 : j      # flip north south due to symmetry
    return nlon_healpix(nside,j)
end

get_npoints(::Type{<:AbstractHEALPixGrid},nside::Integer) = npoints_healpix(nside)
get_lon(::Type{<:AbstractHEALPixGrid},nside::Integer) = Float64[]    # only defined for full grids

function get_colatlons(Grid::Type{<:AbstractHEALPixGrid},nside::Integer)
    npoints = get_npoints(Grid,nside)
    colats_lons = [Healpix.pix2angRing(Healpix.Resolution(nside),i) for i in 1:npoints]
    colats = [colat_lon[1] for colat_lon in colats_lons]
    lons = [colat_lon[2] for colat_lon in colats_lons]
    return colats, lons
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

"""
    G = FullClenshawGrid{T}

A FullClenshawGrid is a regular latitude-longitude grid with an odd number of `nlat` equi-spaced
latitudes, the central latitude ring is on the Equator. The same `nlon` longitudes for every latitude ring.
The grid points are closer in zonal direction around the poles. The values of all grid points are stored
in a vector field `data` that unravels the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct FullClenshawGrid{T} <: AbstractFullGrid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere (incl Equator)

    FullClenshawGrid{T}(data,nlat_half) where T = length(data) == npoints_clenshaw(nlat_half) ?
    new(data,nlat_half) : error("$(length(data))-element Vector{$(eltype(data))} cannot be used to create a "*
        "L$nlat_half ($(4nlat_half)x$(2nlat_half - 1)) FullClenshawGrid{$T}.")
end

# subtract the otherwise double-counted 4nlat_half equator points
npoints_clenshaw(nlat_half::Integer) = 8nlat_half^2 - 4nlat_half
nlat_half_clenshaw(npoints::Integer) = round(Int,1/4 + sqrt(1/16 + npoints/8))  # inverse
# infer resolution parameter nlat_half from length of vector
FullClenshawGrid{T}(data::AbstractVector) where T = FullClenshawGrid(data,nlat_half_clenshaw(length(data)))

truncation_order(::Type{<:FullClenshawGrid}) = 3            # cubic
get_truncation(::Type{<:FullClenshawGrid},nlat_half::Integer) = floor(Int,(4nlat_half-1)/4)
get_resolution(::Type{<:FullClenshawGrid},trunc::Integer) = roundup_fft(ceil(Int,(4*trunc+1)/4))
# get_nresolution() is already implemented for AbstractFullGrid
nlat_odd(::Type{<:FullClenshawGrid}) = true
# get_nlon() is already implemented for AbstractFullGrid
# get_nlon_per_ring() is already implemented for AbstractFullGrid
get_npoints(::Type{<:FullClenshawGrid},nlat_half::Integer) = npoints_clenshaw(nlat_half)
get_colat(::Type{<:FullClenshawGrid},nlat_half::Integer) = [j/(2nlat_half)*π for j in 1:2nlat_half-1]
# get_lon() is already implemented for AbstractFullGrid
# get_colatlons() is already implemented for AbstractFullGrid
# each_index_in_ring() is already implemented for AbstractFullGrid
get_quadrature_weights(::Type{<:FullClenshawGrid},nlat_half::Integer) = clenshaw_curtis_weights(nlat_half)

"""
    G = FullGaussianGrid{T}

A full Gaussian grid is a regular latitude-longitude grid that uses `nlat` Gaussian latitudes,
and the same `nlon` longitudes for every latitude ring. The grid points are closer in zonal direction
around the poles. The values of all grid points are stored in a vector field `v` that unravels
the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct FullGaussianGrid{T} <: AbstractFullGrid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere

    FullGaussianGrid{T}(data,nlat_half) where T = length(data) == 8nlat_half^2 ?
    new(data,nlat_half) : error("$(length(data))-element Vector{$(eltype(data))} cannot be used to create a "*
        "F$nlat_half ($(4nlat_half)x$(2nlat_half)) FullGaussianGrid{$T}.")
end

npoints_gaussian(nlat_half::Integer) = 8nlat_half^2
nlat_half_gaussian(npoints::Integer) = round(Int,sqrt(npoints/8))
# infer resolution parameter nlat_half from length of vector
FullGaussianGrid{T}(data::AbstractVector) where T = FullGaussianGrid(data,nlat_half_gaussian(length(data)))

truncation_order(::Type{<:FullGaussianGrid}) = 2            # quadratic
get_truncation(::Type{<:FullGaussianGrid},nlat_half::Integer) = floor(Int,(4nlat_half-1)/3)
get_resolution(::Type{<:FullGaussianGrid},trunc::Integer) = roundup_fft(ceil(Int,(3*trunc+1)/4))
# get_nresolution() is already defined for AbstractFullGrid
nlat_odd(::Type{<:FullGaussianGrid}) = false
# get_nlon() is already implemented for AbstractFullGrid
# get_nlon_per_ring() is already implemented for AbstractFullGrid
get_npoints(::Type{<:FullGaussianGrid},nlat_half::Integer) = npoints_gaussian(nlat_half)
get_colat(::Type{<:FullGaussianGrid},nlat_half::Integer) =
            π .- acos.(FastGaussQuadrature.gausslegendre(2nlat_half)[1])
# get_lon() is already implemented for AbstractFullGrid
# get_colatlons() is already implemented for AbstractFullGrid
# each_index_in_ring() is already implemented for AbstractFullGrid
get_quadrature_weights(::Type{<:FullGaussianGrid},nlat_half::Integer) = gaussian_weights(nlat_half)

"""
    G = OctahedralGaussianGrid{T}

An Octahedral Gaussian grid that uses `nlat` Gaussian latitudes, but a decreasing number of longitude
points per latitude ring towards the poles. Starting with 20 equi-spaced longitude points (starting at 0˚E)
on the rings around the poles, each latitude ring towards the equator has consecuitively 4 more points,
one for each face of the octahedron. E.g. 20,24,28,32,...nlon-4,nlon,nlon,nlon-4,...,32,28,24,20.
The maximum number of longitue points is `nlon`. The values of all grid points are stored in a vector
field `v` that unravels the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct OctahedralGaussianGrid{T} <: AbstractOctahedralGrid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere

    # check that `nlat_half` match the vector `v` length
    OctahedralGaussianGrid{T}(data,nlat_half) where T = length(data) == npoints_octahedral(nlat_half,false) ?
    new(data,nlat_half) : error("$(length(data))-element Vector{$(eltype(data))}"*
    "cannot be used to create a O$(nlat_half) OctahedralGaussianGrid{$T}.")
end

# number of points and longitudes per ring on the octahedral grid
npoints_octahedral(nlat_half::Integer,nlat_oddp::Bool) =
    nlat_oddp ? max(0,4nlat_half^2 + 32nlat_half - 16) : 4nlat_half^2 + 36nlat_half # max(0,...) needed to avoid negative array size when nlat_half==0
nlat_half_octahedral(npoints::Integer,nlat_oddp::Bool) =
    nlat_oddp ? round(Int,-4+sqrt(20 + npoints/4)) : round(Int,-9/2+sqrt((9/2)^2 + npoints/4))  # inverse
nlon_octahedral(j::Integer) = 16+4j
OctahedralGaussianGrid{T}(data::AbstractVector) where T = OctahedralGaussianGrid(data,nlat_half_octahedral(length(data),false))

truncation_order(::Type{<:OctahedralGaussianGrid}) = 3      # cubic
get_truncation(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = nlat_half-1
get_resolution(::Type{<:OctahedralGaussianGrid},trunc::Integer) = roundup_fft(trunc+1)
# get_nresolution() is already implemented for AbstractOctahedralGrid
nlat_odd(::Type{<:OctahedralGaussianGrid}) = false
# get_nlon() is already implemented for AbstractOctahedralGrid
# get_nlon_per_ring() is already implemented for AbstractOctahedralGrid
get_npoints(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = npoints_octahedral(nlat_half,false)
get_colat(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = get_colat(FullGaussianGrid,nlat_half)
# get_lon() is already implemented for AbstractOctaherdalGrid
# get_colatlons() is already implemented for AbstractOctaherdalGrid
# each_index_in_ring() is already implemented for AbstractOctaherdalGrid
get_quadrature_weights(::Type{<:OctahedralGaussianGrid},nlat_half::Integer) = gaussian_weights(nlat_half)

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
    data::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere (incl Equator)

    # check that `nlat_half` match the vector `v` length
    OctahedralClenshawGrid{T}(data,nlat_half) where T = length(data) == npoints_octahedral(nlat_half,true) ?
    new(data,nlat_half) : error("$(length(data))-element Vector{$(eltype(data))}"*
    "cannot be used to create a O$(nlat_half) OctahedralClenshawGrid{$T}.")
end

OctahedralClenshawGrid{T}(data::AbstractVector) where T = OctahedralClenshawGrid(data,nlat_half_octahedral(length(data),true))

truncation_order(::Type{<:OctahedralClenshawGrid}) = 3      # cubic
get_truncation(::Type{<:OctahedralClenshawGrid},nlat_half::Integer) = nlat_half-1
get_resolution(::Type{<:OctahedralClenshawGrid},trunc::Integer) = roundup_fft(trunc+1)
# get_nresolution() is already implemented for AbstractOctahedralGrid
nlat_odd(::Type{<:OctahedralClenshawGrid}) = true
# get_nlon() is already implemented for AbstractOctahedralGrid
# get_nlon_per_ring() is already implemented for AbstractOctahedralGrid
get_npoints(::Type{<:OctahedralClenshawGrid},nlat_half::Integer) = npoints_octahedral(nlat_half,true)
get_colat(::Type{<:OctahedralClenshawGrid},nlat_half::Integer) = get_colat(FullClenshawGrid,nlat_half)
# get_lon() is already implemented for AbstractOctaherdalGrid
# get_colatlons() is already implemented for AbstractOctaherdalGrid
# each_index_in_ring() is already implemented for AbstractOctaherdalGrid
get_quadrature_weights(::Type{<:OctahedralClenshawGrid},nlat_half::Integer) = clenshaw_curtis_weights(nlat_half)

"""
    H = HEALPixGrid{T}

A HEALPix grid with 12 faces, each `nside`x`nside` grid points, each covering the same area.
The values of all grid points are stored in a vector field `v` that unravels the data 0 to 360˚,
then ring by ring, which are sorted north to south."""
struct HEALPixGrid{T} <: AbstractHEALPixGrid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nside::Int      # nside^2 is the number of pixel in each of the 12 base pixel

    HEALPixGrid{T}(data,nside) where T = length(data) == npoints_healpix(nside) ?
    new(data,nside) : error("$(length(data))-element Vector{$(eltype(data))}"*
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

HEALPixGrid{T}(data::AbstractVector) where T = HEALPixGrid(data,nside_healpix(length(data)))

truncation_order(::Type{<:HEALPixGrid}) = 1                 # linear (in longitude)
get_truncation(::Type{<:HEALPixGrid},nside::Integer) = nside_assert(nside) ? 2nside-1 : nothing
get_resolution(::Type{<:HEALPixGrid},trunc::Integer) = roundup_fft(ceil(Int,(trunc+1)/2),small_primes=[2])
get_colat(G::Type{<:HEALPixGrid},nside::Integer) =
            [acos(Healpix.ring2z(Healpix.Resolution(nside),j)) for j in 1:nlat_healpix(nside)]
# get_lon() is already implemented for AbstractHEALPixGrid
# get_colatlons() is already implemented for AbstractHEALPixGrid
# each_index_in_ring() is already implemented for AbstractHEALPixGrid

# define for all grids that the type T can be infered from the elements in data vector
# whether the resolution parameter n is provided or not (hence the ...)
(Grid::Type{<:AbstractGrid})(data::AbstractVector,n::Integer...) = Grid{eltype(data)}(data,n...)

# convert an AbstractMatrix to the full grids, and vice versa
(Grid::Type{<:AbstractFullGrid})(M::AbstractMatrix{T}) where T = Grid{T}(vec(M))
Base.Matrix(G::AbstractFullGrid{T}) where T = Matrix{T}(reshape(G.data,:,2G.nlat_half - nlat_odd(G)))

# QUADRATURE WEIGHTS
# gaussian_weights are exact for Gaussian latitudes when nlat > (2T+1)/2
# clenshaw_curtis_weights are exact for equi-angle latitudes when nlat > 2T+1
gaussian_weights(nlat_half::Integer) = FastGaussQuadrature.gausslegendre(2nlat_half)[2][1:nlat_half]
function clenshaw_curtis_weights(nlat_half::Integer)
    nlat = 2nlat_half - 1
    θs = get_colat(FullClenshawGrid,nlat_half)
    return [4sin(θj)/(nlat+1)*sum([sin(p*θj)/p for p in 1:2:nlat]) for θj in θs[1:nlat_half]]
end

# SOLID ANGLES ΔΩ = sinθ Δθ Δϕ
get_solid_angles(Grid::Type{<:AbstractGrid},nlat_half::Integer) = 
    get_quadrature_weights(Grid,nlat_half) .* (2π./get_nlons(Grid,nlat_half))
get_solid_angles(Grid::Type{<:HEALPixGrid},nside::Integer) =
    4π/get_npoints(Grid,nside)*ones(get_nlat_half(Grid,nside))
