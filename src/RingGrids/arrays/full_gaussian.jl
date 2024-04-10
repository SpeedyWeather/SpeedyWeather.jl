struct FullGaussianArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractFullGridArray{T, N, ArrayType}
    data::ArrayType                 # data array, ring by ring, north to south
    nlat_half::Int                  # number of latitudes on one hemisphere
    rings::Vector{UnitRange{Int}}   # TODO make same array type as data?
end

const FullGaussianGrid{T} = FullGaussianArray{T, 1, Vector{T}}



nonparametric_type(::Type{<:FullGaussianGrid}) = FullGaussianGrid

npoints_gaussian(nlat_half::Integer) = 8nlat_half^2
nlat_half_gaussian(npoints::Integer) = round(Int, sqrt(npoints/8))

# infer nlat_half from data vector length, infer parametric type from eltype of data
FullGaussianGrid{T}(data::AbstractVector) where T = FullGaussianGrid{T}(data, nlat_half_gaussian(length(data)))
FullGaussianGrid(data::AbstractVector, n::Integer...) = FullGaussianGrid{eltype(data)}(data, n...)

nlat_odd(::Type{<:FullGaussianGrid}) = false
get_npoints(::Type{<:FullGaussianGrid}, nlat_half::Integer) = npoints_gaussian(nlat_half)
get_colat(::Type{<:FullGaussianGrid}, nlat_half::Integer) =
            π .- acos.(FastGaussQuadrature.gausslegendre(2nlat_half)[1])
get_quadrature_weights(::Type{<:FullGaussianGrid}, nlat_half::Integer) = gaussian_weights(nlat_half)
full_grid(::Type{<:FullGaussianGrid}) = FullGaussianGrid    # the full grid with same latitudes