"""
$(TYPEDSIGNATURES)
Check whether elements of a vector `v` are strictly increasing."""
function isincreasing(x::AbstractVector)
    is_increasing = true
    for i in 2:length(x)
        is_increasing &= x[i-1] < x[i] ? true : false
    end
    return is_increasing
end

"""
$(TYPEDSIGNATURES)

Check whether elements of a vector `v` are strictly decreasing."""
function isdecreasing(x::AbstractVector)
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

"""
$(TYPEDSIGNATURES)
Prints to `io` all fields of a struct `A` identified by their
`keys`."""
function print_fields(io::IO,A,keys;arrays::Bool=false)
    keys_filtered = arrays ? keys : filter(key -> ~(getfield(A,key) isa AbstractArray),keys)
    n = length(keys_filtered)
    filtered = n < length(keys)
    for (i,key) in enumerate(keys_filtered)
        last = (i == n) & ~filtered
        key = keys[i]
        val = getfield(A,key)
        ~last ? println(io,"├ $key::$(typeof(val)) = $val") :
                print(io,  "└ $key::$(typeof(val)) = $val")
    end
    if filtered                 # add the names of arrays
        s = "└── arrays: "
        for key in keys
            if ~(key in keys_filtered)
                s *= "$key, "
            end
        end
        print(io,s[1:end-2])    # remove last ", "
    end
end