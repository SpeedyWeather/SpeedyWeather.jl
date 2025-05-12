abstract type AbstractGrid end
abstract type AbstractFullGrid <: AbstractGrid end
abstract type AbstractReducedGrid <: AbstractGrid end

struct FullGaussianGrid{V} <: AbstractFullGrid
    nlat_half::Int
    rings::V
end

struct OctahedralGaussianGrid{V} <: AbstractReducedGrid
    nlat_half::Int
    rings::V
end

abstract type AbstractField{T, N, ArrayType, Grid} <: AbstractArray{T, N} end

struct Field{T, N, ArrayType <: AbstractArray{T, N}, Grid <: AbstractGrid} <: AbstractField{T, N, ArrayType, Grid}
    data::ArrayType
    grid::Grid
end

isreduced(::Type{<:AbstractField}) = true
isreduced(::Type{<:AbstractField{T, N, A, <:AbstractFullGrid}}) where {T, N, A} = false
isreduced(field::Field) = isreduced(typeof(field))
isreduced(::Type{<:AbstractGrid}) = true
isreduced(::Type{<:AbstractFullGrid}) = false
isreduced(grid::AbstractGrid) = isreduced(typeof(grid))
isfull(F) = ~isreduced(F)

Base.size(f::Field, args...) = size(f.data, args...)
Base.@propagate_inbounds Base.getindex(f::Field, args...) = getindex(f.data, args...)