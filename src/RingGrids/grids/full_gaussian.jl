"""A `FullGaussianGrid` is a discretization of the sphere that uses Gaussian latitudes for each ring.
As a "full" grid it has the same number of longitude points on every latitude ring, i.e. it is 
representable as a matrix, with denser grid points towards the poles. The Gaussian latitudes enable to
use the Gaussian quadrature for the spectral transform, hence the name. The resolution parameter of the
horizontal grid is `nlat_half` (number of latitude rings on one hemisphere, Equator included) and the
ring indices are precomputed in `rings`. Note that a `Grid` does not contain any data, it only describes
the discretization of the space, see `Field` for a data on a `Grid`.
$(TYPEDFIELDS)"""
struct FullGaussianGrid{A, V} <: AbstractFullGrid{A}
    nlat_half::Int              # number of latitudes on one hemisphere
    architecture::A             # information about device, CPU/GPU
    rings::V                    # precomputed ring indices
end

nonparametric_type(::Type{<:FullGaussianGrid}) = FullGaussianGrid

# FIELD
const FullGaussianField{T, N} = Field{T, N, Architecture, Grid} where {Architecture, Grid<:FullGaussianGrid}

# define grid_type (i) without T, N, (ii) with T, (iii) with T, N but not with <:FullGaussianField
# to not have precendence over grid_type(::Type{Field{...})
grid_type(::Type{FullGaussianField}) = FullGaussianGrid
grid_type(::Type{FullGaussianField{T}}) where T = FullGaussianGrid
grid_type(::Type{FullGaussianField{T, N}}) where {T, N} = FullGaussianGrid

Base.show(io::IO, F::Type{<:FullGaussianField{T, N}}) where {T, N} = print(io, "FullGaussianField{$T, $N}")

# SIZE
nlat_odd(::Type{<:FullGaussianGrid}) = false        # Gaussian latitudes always even
get_npoints(::Type{<:FullGaussianGrid}, nlat_half::Integer) = 8 * nlat_half^2
get_nlat_half(::Type{<:FullGaussianGrid}, npoints::Integer) = round(Int, sqrt(npoints/8))
get_nlon(::Type{<:FullGaussianGrid}, nlat_half::Integer) = 4nlat_half

## COORDINATES
function get_latd(::Type{<:FullGaussianGrid}, nlat_half::Integer)
    return acosd.(FastGaussQuadrature.gausslegendre(2nlat_half)[1]) .- 90
end

function get_lond(::Type{<:FullGaussianGrid}, nlat_half::Integer)
    nlat_half == 0 && return Float64[]      # necessary to avoid error from /0 below
    nlon = get_nlon(FullGaussianGrid, nlat_half)
    return collect(range(0, 360 - 180/nlon, step=360/nlon))
end

# QUADRATURE
get_quadrature_weights(::Type{<:FullGaussianGrid}, nlat_half::Integer) = gaussian_weights(nlat_half)