"""
    G = FullClenshawGrid{T}

A FullClenshawGrid is a regular latitude-longitude grid with an odd number of `nlat` equi-spaced
latitudes, the central latitude ring is on the Equator. The same `nlon` longitudes for every latitude ring.
The grid points are closer in zonal direction around the poles. The values of all grid points are stored
in a vector field `data` that unravels the data 0 to 360˚, then ring by ring, which are sorted north to south."""
struct FullClenshawGrid{T} <: AbstractFullGrid{T}
    data::Vector{T}    # data vector, ring by ring, north to south
    nlat_half::Int  # number of latitudes on one hemisphere (incl Equator)

    FullClenshawGrid{T}(data::AbstractVector, nlat_half::Integer) where T = length(data) == npoints_clenshaw(nlat_half) ?
    new(data, nlat_half) : error("$(length(data))-element Vector{$(eltype(data))} cannot be used to create a "*
        "L$nlat_half ($(4nlat_half)x$(2nlat_half - 1)) FullClenshawGrid{$T}.")
end

nonparametric_type(::Type{<:FullClenshawGrid}) = FullClenshawGrid

# subtract the otherwise double-counted 4nlat_half equator points
npoints_clenshaw(nlat_half::Integer) = 8nlat_half^2 - 4nlat_half
nlat_half_clenshaw(npoints::Integer) = round(Int, 1/4 + sqrt(1/16 + npoints/8))  # inverse

# infer nlat_half from data vector length, infer parametric type from eltype of data
FullClenshawGrid{T}(data::AbstractVector) where T = FullClenshawGrid{T}(data, nlat_half_clenshaw(length(data)))
FullClenshawGrid(data::AbstractVector, n::Integer...) = FullClenshawGrid{eltype(data)}(data, n...)

nlat_odd(::Type{<:FullClenshawGrid}) = true
get_npoints(::Type{<:FullClenshawGrid}, nlat_half::Integer) = npoints_clenshaw(nlat_half)
get_colat(::Type{<:FullClenshawGrid}, nlat_half::Integer) = [j/(2nlat_half)*π for j in 1:2nlat_half-1]
get_quadrature_weights(::Type{<:FullClenshawGrid}, nlat_half::Integer) = clenshaw_curtis_weights(nlat_half)
full_grid(::Type{<:FullClenshawGrid}) = FullClenshawGrid    # the full grid with same latitudes