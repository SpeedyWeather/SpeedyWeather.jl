"""
    G = FullHEALPixGrid{T}

A full HEALPix grid is a regular latitude-longitude grid that uses `nlat` latitudes from the HEALPix grid,
and the same `nlon` longitudes for every latitude ring. The grid points are closer in zonal direction
around the poles. The values of all grid points are stored in a vector field `v` that unravels
the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct FullHEALPixGrid{T} <: AbstractFullGrid{T}
    data::Vector{T}     # data vector, ring by ring, north to south
    nlat_half::Int      # number of latitudes on one hemisphere

    FullHEALPixGrid{T}(data::AbstractVector, nlat_half::Integer) where T = length(data) == npoints_fullhealpix(nlat_half) ?
    new(data, nlat_half) : error("$(length(data))-element Vector{$(eltype(data))} cannot be used to create a "*
        "H$nlat_half ($(4nlat_half)x$(2nlat_half-1)) FullHEALPixGrid{$T}.")
end

nonparametric_type(::Type{<:FullHEALPixGrid}) = FullHEALPixGrid

npoints_fullhealpix(nlat_half::Integer) = 4nlat_half*(2nlat_half-1)
nlat_half_fullhealpix(npoints::Integer) = round(Int, 1/4 + sqrt(1/16 + npoints/8))

# infer nlat_half from data vector length, infer parametric type from eltype of data
FullHEALPixGrid{T}(data::AbstractVector) where T = FullHEALPixGrid{T}(data, nlat_half_fullhealpix(length(data)))
FullHEALPixGrid(data::AbstractVector, n::Integer...) = FullHEALPixGrid{eltype(data)}(data, n...)
nlat_odd(::Type{<:FullHEALPixGrid}) = true
get_npoints(::Type{<:FullHEALPixGrid}, nlat_half::Integer) = npoints_fullhealpix(nlat_half)
get_colat(::Type{<:FullHEALPixGrid}, nlat_half::Integer) = get_colat(HEALPixGrid, nlat_half)
full_grid(::Type{<:FullHEALPixGrid}) = FullHEALPixGrid      # the full grid with same latitudes
get_quadrature_weights(::Type{<:FullHEALPixGrid}, nlat_half::Integer) = healpix_weights(nlat_half)

