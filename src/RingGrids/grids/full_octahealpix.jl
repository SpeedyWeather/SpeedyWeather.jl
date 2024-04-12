"""
    G = FullOctaHEALPixGrid{T}

A full OctaHEALPix grid is a regular latitude-longitude grid that uses `nlat` OctaHEALPix latitudes,
and the same `nlon` longitudes for every latitude ring. The grid points are closer in zonal direction
around the poles. The values of all grid points are stored in a vector field `v` that unravels
the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct FullOctaHEALPixGrid{T} <: AbstractFullGrid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere

    FullOctaHEALPixGrid{T}(data::AbstractVector, nlat_half::Integer) where T = length(data) == npoints_fulloctahealpix(nlat_half) ?
    new(data, nlat_half) : error("$(length(data))-element Vector{$(eltype(data))} cannot be used to create a "*
        "F$nlat_half ($(4nlat_half)x$(2nlat_half - 1)) FullOctaHEALPixGrid{$T}.")
end

nonparametric_type(::Type{<:FullOctaHEALPixGrid}) = FullOctaHEALPixGrid

npoints_fulloctahealpix(nlat_half::Integer) = 8nlat_half^2 - 4nlat_half
nlat_half_fulloctahealpix(npoints::Integer) = round(Int, 1/4 + sqrt(1/16 + npoints/8))

# infer nlat_half from data vector length, infer parametric type from eltype of data
FullOctaHEALPixGrid{T}(data::AbstractVector) where T = FullOctaHEALPixGrid{T}(data, nlat_half_fulloctahealpix(length(data)))
FullOctaHEALPixGrid(data::AbstractVector, n::Integer...) = FullOctaHEALPixGrid{eltype(data)}(data, n...)

nlat_odd(::Type{<:FullOctaHEALPixGrid}) = true
get_npoints(::Type{<:FullOctaHEALPixGrid}, nlat_half::Integer) = npoints_fulloctahealpix(nlat_half)
get_colat(::Type{<:FullOctaHEALPixGrid}, nlat_half::Integer) = get_colat(OctaHEALPixGrid, nlat_half)
get_quadrature_weights(::Type{<:FullOctaHEALPixGrid}, nlat_half::Integer) = octahealpix_weights(nlat_half)
full_grid(::Type{<:FullOctaHEALPixGrid}) = FullOctaHEALPixGrid    # the full grid with same latitudes