"""
    G = OctahedralGaussianGrid{T}

An Octahedral Gaussian grid that uses `nlat` Gaussian latitudes, but a decreasing number of longitude
points per latitude ring towards the poles. Starting with 20 equi-spaced longitude points (starting at 0˚E)
on the rings around the poles, each latitude ring towards the equator has consecuitively 4 more points,
one for each face of the octahedron. E.g. 20, 24, 28, 32, ...nlon-4, nlon, nlon, nlon-4, ..., 32, 28, 24, 20.
The maximum number of longitue points is `nlon`. The values of all grid points are stored in a vector
field `v` that unravels the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct OctahedralGaussianGrid{T} <: AbstractOctahedralGrid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere

    # check that `nlat_half` match the vector `v` length
    OctahedralGaussianGrid{T}(data::AbstractVector, nlat_half::Integer) where T = length(data) == npoints_octahedral(nlat_half, false) ?
    new(data, nlat_half) : error("$(length(data))-element Vector{$(eltype(data))}"*
    "cannot be used to create a O$(nlat_half) OctahedralGaussianGrid{$T}.")
end

nonparametric_type(::Type{<:OctahedralGaussianGrid}) = OctahedralGaussianGrid

# number of points and longitudes per ring on the octahedral grid
npoints_octahedral(nlat_half::Integer, nlat_oddp::Bool) =
    nlat_oddp ? max(0, 4nlat_half^2 + 32nlat_half - 16) : 4nlat_half^2 + 36nlat_half # max(0, ...) needed to avoid negative array size when nlat_half==0
nlat_half_octahedral(npoints::Integer, nlat_oddp::Bool) =
    nlat_oddp ? round(Int, -4+sqrt(20 + npoints/4)) : round(Int, -9/2+sqrt((9/2)^2 + npoints/4))  # inverse
nlon_octahedral(j::Integer) = 16+4j

# infer nside from data vector length, infer parametric type from eltype of data
OctahedralGaussianGrid{T}(data::AbstractVector) where T = OctahedralGaussianGrid{T}(data, nlat_half_octahedral(length(data), false))
OctahedralGaussianGrid(data::AbstractVector, n::Integer...) = OctahedralGaussianGrid{eltype(data)}(data, n...)

nlat_odd(::Type{<:OctahedralGaussianGrid}) = false
get_npoints(::Type{<:OctahedralGaussianGrid}, nlat_half::Integer) = npoints_octahedral(nlat_half, false)
get_colat(::Type{<:OctahedralGaussianGrid}, nlat_half::Integer) = get_colat(FullGaussianGrid, nlat_half)
get_quadrature_weights(::Type{<:OctahedralGaussianGrid}, nlat_half::Integer) = gaussian_weights(nlat_half)
full_grid(::Type{<:OctahedralGaussianGrid}) = FullGaussianGrid    # the full grid with same latitudes
matrix_size(G::OctahedralGaussianGrid) = (2*(4+G.nlat_half), 2*(4+G.nlat_half+1))

"""
    abstract type AbstractOctahedralGrid{T} <: AbstractGrid{T} end

An `AbstractOctahedralGrid` is a horizontal grid with 16+4i longitude
points on the latitude ring i starting with i=1 around the pole.
Different latitudes can be used, Gaussian latitudes, equi-angle latitdes, or others."""
abstract type AbstractOctahedralGrid{T} <: AbstractGrid{T} end

get_nlon_max(::Type{<:AbstractOctahedralGrid}, nlat_half::Integer) = nlon_octahedral(nlat_half)

function get_nlon_per_ring(Grid::Type{<:AbstractOctahedralGrid}, nlat_half::Integer, j::Integer)
    nlat = get_nlat(Grid, nlat_half)
    @assert 0 < j <= nlat "Ring $j is outside O$nlat_half grid."
    j = j > nlat_half ? nlat - j + 1 : j      # flip north south due to symmetry
    return nlon_octahedral(j)
end

function get_colatlons(Grid::Type{<:AbstractOctahedralGrid}, nlat_half::Integer)
    
    colat = get_colat(Grid, nlat_half)
    nlat = get_nlat(Grid, nlat_half)
    
    npoints = get_npoints(Grid, nlat_half)
    colats = zeros(npoints)                 # preallocate arrays
    lons = zeros(npoints)

    ij = 1                                  # continuous index
    for j in 1:nlat                         # populate arrays ring by ring
        nlon = get_nlon_per_ring(Grid, nlat_half, j)
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

    nlat = get_nlat(Grid, nlat_half)
    @boundscheck 0 < j <= nlat || throw(BoundsError)        # ring index valid?
    if j <= nlat_half                                       # northern hemisphere incl Equator
        index_1st = 2j*(j+7) - 15                           # first in-ring index i
        index_end = 2j*(j+9)                                # last in-ring index i
    else                                                    # southern hemisphere excl Equator
        j = nlat - j + 1                                    # mirror ring index around Equator
        n = get_npoints(Grid, nlat_half) + 1                 # number of grid points + 1
        index_1st = n - 2j*(j+9)                            # count backwards
        index_end = n - (2j*(j+7) - 15)
    end
    return index_1st:index_end                              # range of i's in ring
end

function each_index_in_ring!(   rings::Vector{<:UnitRange{<:Integer}},
                                Grid::Type{<:AbstractOctahedralGrid},
                                nlat_half::Integer) # resolution param

    nlat = length(rings)
    @boundscheck nlat == get_nlat(Grid, nlat_half) || throw(BoundsError)

    index_end = 0
    @inbounds for j in 1:nlat_half                  # North incl Eq only
        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += 16 + 4j                        # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end
    @inbounds for (j, j_mir) in zip( nlat_half+1:nlat,       # South only
                                    nlat-nlat_half:-1:1)    # reverse index

        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += 16 + 4j_mir                    # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end
end