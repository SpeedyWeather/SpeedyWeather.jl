"""$(TYPEDSIGNATURES)
Check whether elements of a vector `v` are strictly increasing."""
function isincreasing(v::AbstractVector)
    return all(diff(v) .> 0)
end

"""$(TYPEDSIGNATURES)
Check whether elements of a vector `v` are strictly decreasing."""
function isdecreasing(v::AbstractVector)
    return all(diff(v) .< 0)
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
        s_short = textwidth(s) > 75 ? first(s, 75) * "..." : s
        ~last ? println(io, "├ " * s_short) : print(io, "└ " * s_short)
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
Returns a human-readable string for a duration given in seconds, rounding to
(days, hours), (hours, minutes), (minutes, seconds), or seconds with
1-2 decimal places for short durations. E.g.
```@example
julia> using SpeedyWeather: readable_secs

julia> readable_secs(12345)
```
"""
function readable_secs(secs::Real)
    millisecs = Dates.Millisecond(round(Int, secs * 1000))
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

"""
    set_scalar!(x::AbstractArray, i, v)

Sets the value of `x[i]` to `v` in a GPU compatible way. This uses `@allowscalar` with
Reactant and `x[i:i] .= v` as a workaround for Enzyme compatibility otherwise.
"""
set_scalar!(x::AbstractArray, i, v) = _set_scalar!(architecture(x), x, i, v)

function _set_scalar!(arch::AbstractArchitecture, x::AbstractArray, i::Integer, v::Number)
    x[i:i] .= v
    return x
end

"""
    componentwise_ifelse(condition::Bool, x, y)

`ifelse(condition, x, y)` that selects real and imaginary part of complex numbers separately.
In kernels, `ifelse` on a `Complex` compiles to a select of an LLVM aggregate (`[2 x float]`)
which Reactant can't raise (`raise = true`), two selects of the components it can.
For other numbers the same as `ifelse`.
TODO: Might be sovled upstream in Reactant in the future 
"""
@inline componentwise_ifelse(condition::Bool, x, y) = ifelse(condition, x, y)
@inline componentwise_ifelse(condition::Bool, x::Complex, y::Complex) =
    complex(ifelse(condition, real(x), real(y)), ifelse(condition, imag(x), imag(y)))

"""
$(TYPEDSIGNATURES)
Fallback for `@maybe_jit` when Reactant is not available. Just calls `f(args...; kwargs...)`."""
_jit(::AbstractArchitecture, f, args...; kwargs...) = f(args...; kwargs...)

"""Optimization patterns of Reactant excluded from every compilation (`@maybe_jit`, the compiled time
loop, see `reactant_compile_options`), to work around Reactant limitations and bugs: add a pattern's
name here, remove it once fixed upstream. None needed with Reactant v0.2.294 (it fixed the
`concat_view_cse` loop of v0.2.293)."""
const REACTANT_EXCLUDED_PASSES = String[]

"""
    reactant_compile_options(; kwargs...)

Compile options (a `Reactant.CompileOptions`) for all compilations with Reactant: the KernelAbstractions
kernels are raised to StableHLO (`raise = true`) and the patterns in `REACTANT_EXCLUDED_PASSES` are
excluded; `kwargs` are passed on to `Reactant.CompileOptions` and can override both. Defined in the
Reactant extension."""
function reactant_compile_options end

"""
    @maybe_jit arch expr

Macro that conditionally applies `Reactant.@jit` based on the architecture.
For `ReactantDevice`, the extension overloads `_jit` to use `@jit`.
For other architectures, just executes the expression directly.

Usage: `@maybe_jit model.architecture initialize!(model.geometry, model)`
"""
macro maybe_jit(arch, expr)
    if expr.head == :call
        f = expr.args[1]
        rest = expr.args[2:end]
        # keyword arguments appear as Expr(:parameters, ...) at the front of rest
        if !isempty(rest) && rest[1] isa Expr && rest[1].head == :parameters
            kwargs = rest[1]
            args = rest[2:end]
            return :(_jit($(esc(arch)), $(esc(f)), $(esc.(args)...); $(esc.(kwargs.args)...)))
        else
            return :(_jit($(esc(arch)), $(esc(f)), $(esc.(rest)...)))
        end
    else
        return esc(expr)
    end
end
