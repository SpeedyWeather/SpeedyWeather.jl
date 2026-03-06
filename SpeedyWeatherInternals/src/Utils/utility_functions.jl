"""$(TYPEDSIGNATURES)
Check whether elements of a vector `v` are strictly increasing."""
function isincreasing(v::AbstractVector)
    is_increasing = true
    for i in 2:length(v)
        is_increasing &= v[i - 1] < v[i] ? true : false
    end
    return is_increasing
end

"""$(TYPEDSIGNATURES)
Check whether elements of a vector `v` are strictly decreasing."""
function isdecreasing(v::AbstractVector)
    is_decreasing = true
    for i in 2:length(v)
        is_decreasing &= v[i - 1] > v[i] ? true : false
    end
    return is_decreasing
end

"""
$(TYPEDSIGNATURES)
Prints to `io` all fields of a struct `A` identified by their
`keys`."""
function print_fields(io::IO, A, keys; arrays::Bool = false, values::Bool = true)
    keys_filtered = arrays ? keys : filter(key -> ~(getfield(A, key) isa AbstractArray), keys)
    n = length(keys_filtered)
    filtered = n < length(keys)
    for (i, key) in enumerate(keys_filtered)
        last = (i == n) & ~filtered
        key = keys_filtered[i]
        val = getfield(A, key)
        val_str = values ? val : ""
        equal_sign = values ? " = " : ""
        s = styled"{info:$key}{note:::$(typeof(val))}$equal_sign$val_str"
        ~last ? println(io, "├ " * s) : print(io, "└ " *  s)
    end
    if filtered                 # add the names of arrays
        s = styled"└ arrays: "
        for key in keys
            if ~(key in keys_filtered)
                s *= styled"{info:$key}, "
            end
        end

        # remove last ", "
        s_without_comma = s[1:prevind(s, findlast(==(','), s))]
        print(io, s_without_comma)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Returns `Dates.CompoundPeriod` rounding to either (days, hours), (hours, minutes), (minutes,
seconds), or seconds with 1 decimal place accuracy for >10s and two for less.
E.g.
```@example
julia> using SpeedyWeather: readable_secs

julia> readable_secs(12345)
```
"""
function readable_secs(secs::Real)
    millisecs = Dates.Millisecond(round(secs * 1000))
    if millisecs >= Dates.Day(1)
        return Dates.canonicalize(round(millisecs, Dates.Hour))
    elseif millisecs >= Dates.Hour(1)
        return Dates.canonicalize(round(millisecs, Dates.Minute))
    elseif millisecs >= Dates.Minute(1)
        return Dates.canonicalize(round(millisecs, Dates.Second))
    elseif millisecs >= Dates.Second(10)
        return Dates.canonicalize(round(millisecs, Dates.Millisecond(100)))
    end
    return Dates.canonicalize(round(millisecs, Dates.Millisecond(10)))
end
