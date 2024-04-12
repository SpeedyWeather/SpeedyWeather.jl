abstract type AbstractReducedGridArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractGridArray{T, N, ArrayType} end
const AbstractReducedGrid{T} = AbstractReducedGridArray{T, 1, Vector{T}}

full_grid(G::Type{<:AbstractReducedGridArray}) = @warn "Please define full_grid(::$(nonparametric_type(G))"
