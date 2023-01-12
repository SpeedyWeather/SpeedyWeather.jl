"""
    abstract type AbstractFullGrid{T} <: AbstractGrid{T} end

An `AbstractFullGrid` is a horizontal grid with a constant number of longitude
points across latitude rings. Different latitudes can be used, Gaussian latitudes,
equi-angle latitdes, or others."""
abstract type AbstractFullGrid{T} <: AbstractGrid{T} end

get_nlon(Grid::Type{<:AbstractFullGrid},nlat_half::Integer) = get_nlon_max(Grid,nlat_half)
get_nlon_max(::Type{<:AbstractFullGrid},nlat_half::Integer) = 4nlat_half
get_nlon_per_ring(Grid::Type{<:AbstractFullGrid},nlat_half::Integer,j::Integer) = 
    get_nlon_max(Grid,nlat_half)

function get_lon(Grid::Type{<:AbstractFullGrid},nlat_half::Integer)
    nlon = get_nlon(Grid,nlat_half)
    return collect(range(0,2π-π/nlon,step=2π/nlon))
end

function get_lond(Grid::Type{<:AbstractFullGrid},nlat_half::Integer)
    lon = get_lon(Grid,nlat_half)
    lon .*= 360/2π      # convert to lond in-place
    return lon          # = lond
end

# convert an AbstractMatrix to the full grids, and vice versa
(Grid::Type{<:AbstractFullGrid})(M::AbstractMatrix{T}) where T = Grid{T}(vec(M))
Base.Matrix(grid::AbstractFullGrid{T}) where T = Matrix{T}(reshape(grid.data,:,get_nlat(grid)))
matrix_size(grid::AbstractFullGrid) = (get_nlon_max(grid),get_nlat(grid))
matrix_size(Grid::Type{<:AbstractFullGrid},n::Integer) = (get_nlon_max(Grid,n),get_nlat(Grid,n))

function get_colatlons(Grid::Type{<:AbstractFullGrid},nlat_half::Integer)

    colat = get_colat(Grid,nlat_half)       # vector of colats [0,π]
    lon = get_lon(Grid,nlat_half)           # vector of longitudes [0,2π)
    nlon = get_nlon(Grid,nlat_half)         # number of longitudes
    nlat = get_nlat(Grid,nlat_half)         # number of latitudes

    npoints = get_npoints(Grid,nlat_half)   # total number of grid points
    colats = zeros(npoints)                 # preallocate
    lons = zeros(npoints)

    for j in 1:nlat                         # populate preallocated colats,lons
        for i in 1:nlon
            ij = i + (j-1)*nlon             # continuous index ij
            colats[ij] = colat[j]
            lons[ij] = lon[i]
        end
    end

    return colats,lons
end

function each_index_in_ring(Grid::Type{<:AbstractFullGrid},     # function for full grids
                            j::Integer,                         # ring index north to south
                            nlat_half::Integer)                 # resolution param

    @boundscheck 0 < j <= get_nlat(Grid,nlat_half) || throw(BoundsError)    # valid ring index?
    nlon = 4nlat_half               # number of longitudes per ring (const)
    index_1st = (j-1)*nlon + 1      # first in-ring index i
    index_end = j*nlon              # last in-ring index i  
    return index_1st:index_end      # range of js in ring
end

function each_index_in_ring!(   rings::Vector{<:UnitRange{<:Integer}},
                                Grid::Type{<:AbstractFullGrid},
                                nlat_half::Integer)
    nlat = length(rings)                # number of latitude rings
    @boundscheck nlat == get_nlat(Grid,nlat_half) || throw(BoundsError)

    nlon = get_nlon(Grid,nlat_half)     # number of longitudes
    index_end = 0                       
    @inbounds for j in 1:nlat
        index_1st = index_end + 1       # 1st index is +1 from prev ring's last index
        index_end += nlon               # only calculate last index per ring
        rings[j] = index_1st:index_end  # write UnitRange to rings vector
    end
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

    FullClenshawGrid{T}(data::AbstractVector,nlat_half::Integer) where T = length(data) == npoints_clenshaw(nlat_half) ?
    new(data,nlat_half) : error("$(length(data))-element Vector{$(eltype(data))} cannot be used to create a "*
        "L$nlat_half ($(4nlat_half)x$(2nlat_half - 1)) FullClenshawGrid{$T}.")
end

# subtract the otherwise double-counted 4nlat_half equator points
npoints_clenshaw(nlat_half::Integer) = 8nlat_half^2 - 4nlat_half
nlat_half_clenshaw(npoints::Integer) = round(Int,1/4 + sqrt(1/16 + npoints/8))  # inverse

# infer nlat_half from data vector length, infer parametric type from eltype of data
FullClenshawGrid{T}(data::AbstractVector) where T = FullClenshawGrid{T}(data,nlat_half_clenshaw(length(data)))
FullClenshawGrid(data::AbstractVector,n::Integer...) = FullClenshawGrid{eltype(data)}(data,n...)

truncation_order(::Type{<:FullClenshawGrid}) = 3            # cubic
get_truncation(::Type{<:FullClenshawGrid},nlat_half::Integer) = floor(Int,(4nlat_half-1)/4)
get_resolution(::Type{<:FullClenshawGrid},trunc::Integer) = roundup_fft(ceil(Int,(4*trunc+1)/4))
nlat_odd(::Type{<:FullClenshawGrid}) = true
get_npoints(::Type{<:FullClenshawGrid},nlat_half::Integer) = npoints_clenshaw(nlat_half)
get_colat(::Type{<:FullClenshawGrid},nlat_half::Integer) = [j/(2nlat_half)*π for j in 1:2nlat_half-1]
get_quadrature_weights(::Type{<:FullClenshawGrid},nlat_half::Integer) = clenshaw_curtis_weights(nlat_half)
full_grid(::Type{<:FullClenshawGrid}) = FullClenshawGrid    # the full grid with same latitudes

"""
    G = FullGaussianGrid{T}

A full Gaussian grid is a regular latitude-longitude grid that uses `nlat` Gaussian latitudes,
and the same `nlon` longitudes for every latitude ring. The grid points are closer in zonal direction
around the poles. The values of all grid points are stored in a vector field `v` that unravels
the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct FullGaussianGrid{T} <: AbstractFullGrid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere

    FullGaussianGrid{T}(data::AbstractVector,nlat_half::Integer) where T = length(data) == 8nlat_half^2 ?
    new(data,nlat_half) : error("$(length(data))-element Vector{$(eltype(data))} cannot be used to create a "*
        "F$nlat_half ($(4nlat_half)x$(2nlat_half)) FullGaussianGrid{$T}.")
end

npoints_gaussian(nlat_half::Integer) = 8nlat_half^2
nlat_half_gaussian(npoints::Integer) = round(Int,sqrt(npoints/8))

# infer nlat_half from data vector length, infer parametric type from eltype of data
FullGaussianGrid{T}(data::AbstractVector) where T = FullGaussianGrid{T}(data,nlat_half_gaussian(length(data)))
FullGaussianGrid(data::AbstractVector,n::Integer...) = FullGaussianGrid{eltype(data)}(data,n...)

truncation_order(::Type{<:FullGaussianGrid}) = 2            # quadratic
get_truncation(::Type{<:FullGaussianGrid},nlat_half::Integer) = floor(Int,(4nlat_half-1)/3)
get_resolution(::Type{<:FullGaussianGrid},trunc::Integer) = roundup_fft(ceil(Int,(3*trunc+1)/4))
nlat_odd(::Type{<:FullGaussianGrid}) = false
get_npoints(::Type{<:FullGaussianGrid},nlat_half::Integer) = npoints_gaussian(nlat_half)
get_colat(::Type{<:FullGaussianGrid},nlat_half::Integer) =
            π .- acos.(FastGaussQuadrature.gausslegendre(2nlat_half)[1])
get_quadrature_weights(::Type{<:FullGaussianGrid},nlat_half::Integer) = gaussian_weights(nlat_half)
full_grid(::Type{<:FullGaussianGrid}) = FullGaussianGrid    # the full grid with same latitudes

"""
    G = FullHEALPixGrid{T}

A full HEALPix grid is a regular latitude-longitude grid that uses `nlat` latitudes from the HEALPix grid,
and the same `nlon` longitudes for every latitude ring. The grid points are closer in zonal direction
around the poles. The values of all grid points are stored in a vector field `v` that unravels
the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct FullHEALPixGrid{T} <: AbstractFullGrid{T}
    data::Vector{T}     # data vector, ring by ring, north to south
    nlat_half::Int      # number of latitudes on one hemisphere

    FullHEALPixGrid{T}(data::AbstractVector,nlat_half::Integer) where T = length(data) == npoints_fullhealpix(nlat_half) ?
    new(data,nlat_half) : error("$(length(data))-element Vector{$(eltype(data))} cannot be used to create a "*
        "H$nlat_half ($(4nlat_half)x$(2nlat_half-1)) FullHEALPixGrid{$T}.")
end

npoints_fullhealpix(nlat_half::Integer) = 4nlat_half*(2nlat_half-1)
nlat_half_fullhealpix(npoints::Integer) = round(Int,1/4 + sqrt(1/16 + npoints/8))

# infer nlat_half from data vector length, infer parametric type from eltype of data
FullHEALPixGrid{T}(data::AbstractVector) where T = FullHEALPixGrid{T}(data,nlat_half_fullhealpix(length(data)))
FullHEALPixGrid(data::AbstractVector,n::Integer...) = FullHEALPixGrid{eltype(data)}(data,n...)
truncation_order(::Type{<:FullHEALPixGrid}) = 2            # quadratic
get_truncation(::Type{<:FullHEALPixGrid},nlat_half::Integer) = get_truncation(HEALPixGrid,nlat_half)
get_resolution(::Type{<:FullHEALPixGrid},trunc::Integer) = get_resolution(HEALPixGrid,trunc)
nlat_odd(::Type{<:FullHEALPixGrid}) = true
get_npoints(::Type{<:FullHEALPixGrid},nlat_half::Integer) = npoints_fullhealpix(nlat_half)
get_colat(::Type{<:FullHEALPixGrid},nlat_half::Integer) = get_colat(HEALPixGrid,nlat_half)
full_grid(::Type{<:FullHEALPixGrid}) = FullHEALPixGrid      # the full grid with same latitudes
get_quadrature_weights(::Type{<:FullHEALPixGrid},nlat_half::Integer) = healpix_weights(nlat_half)

"""
    G = FullHEALPix4Grid{T}

A full HEALPix4 grid is a regular latitude-longitude grid that uses `nlat` HEALPix4 latitudes,
and the same `nlon` longitudes for every latitude ring. The grid points are closer in zonal direction
around the poles. The values of all grid points are stored in a vector field `v` that unravels
the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct FullHEALPix4Grid{T} <: AbstractFullGrid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere

    FullHEALPix4Grid{T}(data::AbstractVector,nlat_half::Integer) where T = length(data) == npoints_fullhealpix4(nlat_half) ?
    new(data,nlat_half) : error("$(length(data))-element Vector{$(eltype(data))} cannot be used to create a "*
        "F$nlat_half ($(4nlat_half)x$(2nlat_half - 1)) FullHEALPix4Grid{$T}.")
end
npoints_fullhealpix4(nlat_half::Integer) = 8nlat_half^2 - 4nlat_half
nlat_half_fullhealpix4(npoints::Integer) = round(Int,1/4 + sqrt(1/16 + npoints/8))

# infer nlat_half from data vector length, infer parametric type from eltype of data
FullHEALPix4Grid{T}(data::AbstractVector) where T = FullHEALPix4Grid{T}(data,nlat_half_fullhealpix4(length(data)))
FullHEALPix4Grid(data::AbstractVector,n::Integer...) = FullHEALPix4Grid{eltype(data)}(data,n...)

truncation_order(::Type{<:FullHEALPix4Grid}) = 3            # cubic
get_truncation(::Type{<:FullHEALPix4Grid},nlat_half::Integer) = floor(Int,(4nlat_half-1)/4)
get_resolution(::Type{<:FullHEALPix4Grid},trunc::Integer) = roundup_fft(ceil(Int,(4*trunc+1)/4))
nlat_odd(::Type{<:FullHEALPix4Grid}) = true
get_npoints(::Type{<:FullHEALPix4Grid},nlat_half::Integer) = npoints_fullhealpix4(nlat_half)
get_colat(::Type{<:FullHEALPix4Grid},nlat_half::Integer) = get_colat(HEALPix4Grid,nlat_half)
get_quadrature_weights(::Type{<:FullHEALPix4Grid},nlat_half::Integer) = healpix4_weights(nlat_half)
full_grid(::Type{<:FullHEALPix4Grid}) = FullHEALPix4Grid    # the full grid with same latitudes