struct FullGaussianArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractFullGridArray{T, N, ArrayType}
    data::ArrayType                 # data array, ring by ring, north to south
    nlat_half::Int                  # number of latitudes on one hemisphere
    rings::Vector{UnitRange{Int}}   # TODO make same array type as data?

    FullGaussianArray(data::A, nlat_half, rings) where {A <: AbstractArray{T, N}} where {T, N} =
        check_inputs(data, nlat_half, rings, FullGaussianArray) ?
        new{T, N, A}(data, nlat_half, rings) :
        error_message(data, nlat_half, rings, FullGaussianArray, T, N, A)
end

# TYPES
const FullGaussianGrid{T} = FullGaussianArray{T, 1, Vector{T}}
nonparametric_type(::Type{<:FullGaussianArray}) = FullGaussianArray
horizontal_grid_type(::Type{<:FullGaussianArray}) = FullGaussianGrid

# SIZE
nlat_odd(::Type{<:FullGaussianArray}) = false
get_npoints2D(::Type{<:FullGaussianArray}, nlat_half::Integer) = 8 * nlat_half^2
get_nlat_half(::Type{<:FullGaussianArray}, npoints2D::Integer) = round(Int, sqrt(npoints2D/8))
get_nlon(::Type{<:FullGaussianArray}, nlat_half::Integer) = 4nlat_half

## COORDINATES
function get_colat(::Type{<:FullGaussianArray}, nlat_half::Integer)
    return π .- acos.(FastGaussQuadrature.gausslegendre(2nlat_half)[1])
end

function get_lon(::Type{<:FullGaussianArray}, nlat_half::Integer)
    nlat_half == 0 && return Float64[]
    nlon = get_nlon(FullGaussianArray, nlat_half)
    return collect(range(0, 2π-π/nlon, step=2π/nlon))
end

# QUADRATURE
get_quadrature_weights(::Type{<:FullGaussianArray}, nlat_half::Integer) = gaussian_weights(nlat_half)