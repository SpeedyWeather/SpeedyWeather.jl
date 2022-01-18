"""Check whether elements of a vector are increasing."""
function isincreasing(x::Vector)
    is_increasing = true
    for i in 2:length(x)
        is_increasing &= x[i] > x[i-1] ? true : false
    end
    return is_increasing
end

"""Set all negative entries in an array to zero."""
function clip_negatives!(A::AbstractArray{T}) where T
    @inbounds for i in eachindex(A)
        A[i] = max(A[i],zero(T))
    end
end