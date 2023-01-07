"""
    abstract type AbstractHEALPixGrid{T} <: AbstractGrid{T} end

An `AbstractHEALPixGrid` is a horizontal grid similar to the standard HEALPixGrid,
but different latitudes can be used, the default HEALPix latitudes or others."""
abstract type AbstractHEALPixGrid{T} <: AbstractGrid{T} end

npoints_healpix(nlat_half::Integer) = 3 * nlat_half^2
nside_healpix(nlat_half::Integer) = nlat_half ÷ 2
nlat_half_healpix(npoints::Integer) = round(Int, sqrt(npoints / 3))  # inverse of npoints_healpix
nlon_healpix(nlat_half::Integer, j::Integer) = min(4j, 2nlat_half, 8nlat_half - 4j)
nlon_max_healpix(nlat_half::Integer) = 2nlat_half

nlat_odd(::Type{<:AbstractHEALPixGrid}) = true
function get_nlon_max(::Type{<:AbstractHEALPixGrid}, nlat_half::Integer)
    nlon_max_healpix(nlat_half)
end

function get_nlon_per_ring(G::Type{<:AbstractHEALPixGrid}, nlat_half::Integer, j::Integer)
    nlat = get_nlat(G, nlat_half)
    @assert 0<j<=nlat "Ring $j is outside H$nlat_half grid."
    return nlon_healpix(nlat_half, j)
end

get_npoints(::Type{<:AbstractHEALPixGrid}, nlat_half::Integer) = npoints_healpix(nlat_half)
get_lon(::Type{<:AbstractHEALPixGrid}, nlat_half::Integer) = Float64[]    # only defined for full grids

function get_colatlons(Grid::Type{<:AbstractHEALPixGrid}, nlat_half::Integer)
    nlat = get_nlat(Grid, nlat_half)
    npoints = get_npoints(Grid, nlat_half)
    nside = nside_healpix(nlat_half)
    colat = get_colat(Grid, nlat_half)

    colats = zeros(npoints)
    lons = zeros(npoints)

    ij = 1
    for j in 1:nlat
        nlon = get_nlon_per_ring(Grid, nlat_half, j)

        # s = 1 for polar caps, s=2,1,2,1,... in the equatorial zone
        s = (j < nside) || (j >= 3nside) ? 1 : ((j - nside) % 2 + 1)
        lon = [π / (nlon ÷ 2) * (i - s / 2) for i in 1:nlon]

        colats[ij:(ij + nlon - 1)] .= colat[j]
        lons[ij:(ij + nlon - 1)] .= lon

        ij += nlon
    end

    return colats, lons
end

function each_index_in_ring(::Type{<:AbstractHEALPixGrid},  # function for HEALPix grids
                            j::Integer,                     # ring index north to south
                            nlat_half::Integer)             # resolution param
    nside = nside_healpix(nlat_half)

    @boundscheck 0 < j < 4nside || throw(BoundsError)       # ring index valid?
    if j < nside                                            # northern polar cap
        index_1st = 2j * (j - 1) + 1                            # first in-ring index i
        index_end = 2j * (j + 1)                                # last in-ring index i
    elseif j <= 3nside                                      # equatorial zone with const nlon
        n = 2nside^2 - 2nside                                 # points in polar cap
        nlon = 4nside                                       # points on latitude rings
        j = j - nside + 1                                   # offset ring index into eq zone
        index_1st = n + (j - 1) * nlon + 1                      # add constant nlon per ring
        index_end = n + j * nlon
    else                                                    # southern polar cap
        n = 12nside^2                                       # total number of points
        j = 4nside - j                                      # count ring index from south pole
        index_1st = n - 2j * (j + 1) + 1                        # count backwards
        index_end = n - 2j * (j - 1)
    end
    return index_1st:index_end                              # range of i's in ring
end

function each_index_in_ring!(rings::Vector{<:UnitRange{<:Integer}},
                             Grid::Type{<:AbstractHEALPixGrid},
                             nlat_half::Integer) # resolution param
    nlat = length(rings)
    @boundscheck nlat == get_nlat(Grid, nlat_half) || throw(BoundsError)

    index_end = 0
    nside = nside_healpix(nlat_half)                # side length of a basepixel

    # North polar cap
    @inbounds for j in 1:(nside - 1)
        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += 4j                             # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end

    # Equatorial belt
    nlon_max = get_nlon_max(Grid, nlat_half)         # number of grid points on belt
    @inbounds for j in nside:(3nside)
        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += nlon_max                       # nlon constant in belt
        rings[j] = index_1st:index_end              # turn into UnitRange
    end

    # South polar cap
    @inbounds for (j, j_mir) in zip((3nside + 1):nlat,  # South only
                                    (nside - 1):-1:1)   # mirror index
        index_1st = index_end + 1                   # 1st index is +1 from prev ring's last index
        index_end += 4j_mir                         # add number of grid points per ring
        rings[j] = index_1st:index_end              # turn into UnitRange
    end
end

"""
    H = HEALPixGrid{T}

A HEALPix grid with 12 faces, each `nside`x`nside` grid points, each covering the same area.
The number of latitude rings on one hemisphere (incl Equator) `nlat_half` is used as resolution parameter.
The values of all grid points are stored in a vector field `v` that unravels the data 0 to 360˚,
then ring by ring, which are sorted north to south."""
struct HEALPixGrid{T} <: AbstractHEALPixGrid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int     # number of latitude rings on one hemisphere (Equator included)

    function HEALPixGrid{T}(data, nlat_half) where {T}
        length(data) == npoints_healpix(nlat_half) ?
        new(data, nlat_half) :
        error("$(length(data))-element Vector{$(eltype(data))}" *
              "cannot be used to create an H$nlat_half HEALPixGrid{$T}.")
    end
end

# infer nlat_half from data vector length, infer parametric type from eltype of data
function HEALPixGrid{T}(data::AbstractVector) where {T}
    HEALPixGrid{T}(data, nlat_half_healpix(length(data)))
end
HEALPixGrid(data::AbstractVector, n::Integer...) = HEALPixGrid{eltype(data)}(data, n...)

truncation_order(::Type{<:HEALPixGrid}) = 1                 # linear (in longitude)
get_truncation(::Type{<:HEALPixGrid}, nlat_half::Integer) = nlat_half - 1
get_resolution(::Type{<:HEALPixGrid}, trunc::Integer) = trunc + 1

function get_colat(::Type{<:HEALPixGrid}, nlat_half::Integer)
    nlat = get_nlat(HEALPixGrid, nlat_half)
    nside = nside_healpix(nlat_half)
    colat = zeros(nlat)

    for j in 1:nside
        colat[j] = acos(1 - j^2 / 3nside^2)
    end     # north polar cap
    for j in (nside + 1):(3nside)
        colat[j] = acos(4 / 3 - 2j / 3nside)
    end     # equatorial belt
    for j in (3nside + 1):nlat
        colat[j] = acos((2nlat_half - j)^2 / 3nside^2 - 1)
    end     # south polar cap

    return colat
end

full_grid(::Type{<:HEALPixGrid}) = FullHEALPixGrid    # the full grid with same latitudes
function matrix_size(G::HEALPixGrid)
    nside = nside_healpix(G.nlat_half)
    return (5nside, 5nside)
end
