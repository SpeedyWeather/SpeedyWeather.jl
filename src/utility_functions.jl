"""
    true/false = isincreasing(v::Vector)

Check whether elements of a vector `v` are strictly increasing."""
function isincreasing(x::Vector)
    is_increasing = true
    for i in 2:length(x)
        is_increasing &= x[i-1] < x[i] ? true : false
    end
    return is_increasing
end

"""
    true/false = isdecreasing(v::Vector)

Check whether elements of a vector `v` are strictly decreasing."""
function isdecreasing(x::Vector)
    is_decreasing = true
    for i in 2:length(x)
        is_decreasing &= x[i-1] > x[i] ? true : false
    end
    return is_decreasing
end

"""
    clip_negatives!(A::AbstractArray)

Set all negative entries `a` in `A` to zero."""
function clip_negatives!(A::AbstractArray{T}) where T
    @inbounds for i in eachindex(A)
        A[i] = max(A[i],zero(T))
    end
end

"""
    underflow!(A::AbstractArray,ϵ::Real)

Underflows element `a` in `A` to zero if `abs(a) < ϵ`."""
function underflow!(A::AbstractArray{T},ϵ::Real) where T
    ϵT = convert(T,abs(ϵ))
    @inbounds for i in eachindex(A)
        A[i] = abs(A[i]) < ϵT ? zero(T) : A[i]
    end
end

"""
    flipgsign!(A::AbstractArray)

Like `-A` but in-place."""
function flipsign!(A::AbstractArray)
    @inbounds for i in eachindex(A)
        A[i] = -A[i]
    end
    A
end

# define as will only become availale in Julia 1.9
pkgversion(m::Module) = VersionNumber(TOML.parsefile(joinpath(
    dirname(string(first(methods(m.eval)).file)), "..", "Project.toml"))["version"])

# NAN initialisation
"""
    A = nans(T,dims...)

Allocate array A with NaNs of type T. Similar to zeros(T,dims...)."""
function nans(::Type{T},dims...) where T
    return fill(convert(T,NaN),dims...)
end

"""
    A = nans(dims...)

Allocate A::Array{Float64} with NaNs."""
nans(dims...) = nans(Float64,dims...)