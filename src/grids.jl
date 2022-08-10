abstract type AbstractGrid{T} <: AbstractArray{T} end

struct FullGaussianGrid{T} <: AbstractGrid{T}
    v::Vector{T}    # data vector, ring by ring, north to south
    nlon::Int       # number of longitudes
    nlat::Int       # number of latitudes

    FullGaussianGrid{T}(v,nlon,nlat) where T = length(v) == nlon*nlat ?
    new(v,nlon,nlat) : error("$(length(v))-element Vector{$(eltype(v))} cannot be used to create a "*
        "$(m)x$(n) FullGaussianGrid{$T}.")
end

struct OctahedralGaussianGrid{T} <: AbstractGrid{T}
    v::Vector{T}    # data vector, ring by ring, north to south
    nlon::Int       # number of longitudes at Equator (20 at poles, +4 for every ring towards Equator)
    nlat::Int       # number of latitudes
end

struct HEALPixGrid{T} <: AbstractGrid{T}
    v::Vector{T}    # data vector, ring by ring, north to south
    nlon::Int       # number of longitudes at Equator (20 at poles, +4 for every ring towards Equator)
    nlat::Int       # number of latitudes
    nside::Int      # nside^2 is the number of pixel in each of the 12 base pixel
end