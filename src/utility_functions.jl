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

"""
    time_string = readable_secs(secs::Real)

Returns a human readable string representing seconds in terms of days, hours, minutes or seconds,
rounding to either (days, hours), (hours, minutes), (minutes, seconds), or seconds with 1 decimal
place accuracy for >10s and two for less.
E.g. 
```julia
julia> readable_secs(12345)
"3h, 25min"
```
"""
function readable_secs(secs::Real)
    days = floor(Int,secs/3600/24)
    hours = floor(Int,(secs/3600) % 24)
    minutes = floor(Int,(secs/60) % 60)
    seconds = floor(Int,secs%3600%60)
    secs1f = @sprintf "%.1fs" secs%3600%60
    secs2f = @sprintf "%.2fs" secs%3600%60

    days > 0 && return "$(days)d, $(hours)h"
    hours > 0 && return "$(hours)h, $(minutes)min"
    minutes > 0 && return "$(minutes)min, $(seconds)s"
    seconds > 10 && return secs1f
    return secs2f
end