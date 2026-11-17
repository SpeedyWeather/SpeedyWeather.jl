"""
    ProgressMeters

A lightweight progress meter, a progress line that is redrawn in place while a loop runs. The
line is built from *elements*, each printing one value (e.g. the percentage or the estimated remaining
time), joined by a separator, e.g.

```julia
p = Progress(100; description = "Computing:", elements = (Percentage(), Bar(), ETA(), Counter()))
for i in 1:100
    work()
    next!(p)
end
```

prints `Computing:  42% ━━━━━━━━━━━━              ETA: 0:00:07, 42`. Define a custom element by
subtyping [`AbstractProgressElement`](@ref) and extending [`print_element`](@ref).
"""
module ProgressMeters

using DocStringExtensions: TYPEDSIGNATURES, TYPEDFIELDS
import Printf: @sprintf

export Progress, next!, finish!,
    AbstractProgressElement, print_element, Percentage, Bar, ETA, ElapsedTime, Counter

"""Supertype of the elements of a progress line, see [`print_element`](@ref)."""
abstract type AbstractProgressElement end

"""
Progress meter for a loop of `n` iterations, create with
`Progress(n; description, elements, separator, color, enabled, dt, output)`, advance with
[`next!`](@ref), end with [`finish!`](@ref). Fields are
$(TYPEDFIELDS)"""
mutable struct Progress
    "Total number of iterations"
    n::Int

    "Iterations so far"
    counter::Int

    "Printed at the start of the line (followed by a space) unless empty"
    description::String

    "Elements of the line, see `print_element`"
    elements::Tuple

    "Separator between elements, a single space is used next to a `Bar`"
    separator::String

    "Color of the line, a key of `Base.text_colors`, only used if `output` has `:color`"
    color::Symbol

    "Draw the line? Iterations are still counted otherwise"
    enabled::Bool

    "Minimum wallclock time [s] between redraws"
    dt::Float64

    "Stream the line is printed to"
    output::IO

    "time() of creation"
    tinit::Float64

    "time() after the last redraw"
    tlast::Float64

    "time() of the current redraw, for the elements"
    tcurrent::Float64

    "Has the line been drawn at least once?"
    printed::Bool

    "Has the final line been drawn (or the meter finished without drawing)?"
    finished::Bool
end

"""$(TYPEDSIGNATURES)
Progress meter for `n` iterations. The line is drawn to `output` at most every `dt` seconds by
[`next!`](@ref) and a last time when `n` is reached (only if drawn before, so loops faster than `dt`
print nothing). `elements` default to [`default_elements`](@ref)."""
function Progress(
        n::Integer;
        description::AbstractString = "",
        elements = default_elements(),
        separator::AbstractString = ", ",
        color::Symbol = :blue,
        enabled::Bool = true,
        dt::Real = 0.1,
        output::IO = stderr,
    )
    t = time()
    return Progress(
        n, 0, description, Tuple(elements), separator, color, enabled, dt, output,
        t, t, t, false, false
    )
end

function Base.show(io::IO, p::Progress)
    return print(io, "Progress($(p.counter)/$(p.n), elements = $(p.elements))")
end

"""$(TYPEDSIGNATURES)
Advance the progress meter `p` by `step` iterations and redraw its line if more than `p.dt` seconds
have passed since the last redraw, or draw the final line once `p.n` is reached. Does not allocate
unless the line is drawn."""
function next!(p::Progress; step::Integer = 1)
    p.counter += step
    p.enabled || return nothing
    p.finished && return nothing
    if p.counter >= p.n
        p.printed ? draw!(p, true) : (p.finished = true)
    elseif time() > p.tlast + p.dt
        draw!(p, false)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Set the progress meter `p` to its last iteration and draw the final line if the line has been drawn before."""
function finish!(p::Progress)
    p.counter = max(p.counter, p.n)
    p.enabled || return nothing
    p.finished && return nothing
    p.printed ? draw!(p, true) : (p.finished = true)
    return nothing
end

# draw the line over the current one (or the final line, followed by a newline), written at once
function draw!(p::Progress, finished::Bool)
    p.tcurrent = time()
    p.finished = finished
    buffer = IOBuffer()
    print(buffer, '\r')
    colored = get(p.output, :color, false)::Bool
    colored && print(buffer, get(Base.text_colors, p.color, ""))
    print(buffer, render_line(p))
    colored && print(buffer, Base.text_colors[:normal])
    print(buffer, "\e[K")       # clear the rest of a previous, longer line
    finished && print(buffer, '\n')
    write(p.output, take!(buffer))
    flush(p.output)
    p.tlast = time()
    p.printed = true
    return nothing
end

"""
    print_element(element, p::Progress) -> String

Text of `element` in the line of the progress meter `p`, only its value, the separators are added
by the meter. Called for every element on each redraw, so this is where an element computes what it
shows. Useful fields of `p` are `counter`, `n`, `tinit` (`time()` of creation), `tcurrent` (`time()`
of this redraw) and `finished` (whether this is the final line). Strings print themselves, an empty
string leaves the element out of the line."""
function print_element end
print_element(s::AbstractString, p::Progress) = s

"""$(TYPEDSIGNATURES)
The line of the progress meter `p` for its `elements`: the description, then the non-empty texts of the
elements joined by `p.separator`, with a single space next to a bar instead."""
function render_line(p::Progress)
    (; elements, separator, description) = p
    texts = String[isbar(element) ? "" : print_element(element, p) for element in elements]

    # bars without a length share the terminal width the rest of the line leaves
    flexible = [isbar(element) && isnothing(element.length) for element in elements]
    bar_width = 0
    if any(flexible)
        visible = [!isempty(text) || isbar(element) for (element, text) in zip(elements, texts)]
        width = textwidth(description) + !isempty(description) +
            sum(textwidth, texts; init = 0) + separators_width(elements, visible, separator) +
            sum(element -> isbar(element) ? something(element.length, 0) : 0, elements; init = 0)
        bar_width = max(0, (displaysize(p.output)[2] - width - 1) ÷ count(flexible))
    end

    for (i, element) in enumerate(elements)
        if isbar(element)
            texts[i] = bar_string(element, p, something(element.length, bar_width))
        end
    end

    line = IOBuffer()
    print(line, description)
    previous = nothing      # previous element that was printed
    for (element, text) in zip(elements, texts)
        isempty(text) && continue
        if !isnothing(previous)
            print(line, isbar(previous) || isbar(element) ? " " : separator)
        elseif !isempty(description)
            print(line, " ")
        end
        print(line, text)
        previous = element
    end
    return String(take!(line))
end

# width of all separators if only the elements with `visible` are printed
function separators_width(elements, visible, separator)
    width = 0
    previous = nothing
    for (element, isvisible) in zip(elements, visible)
        isvisible || continue
        if !isnothing(previous)
            width += isbar(previous) || isbar(element) ? 1 : textwidth(separator)
        end
        previous = element
    end
    return width
end

"""$(TYPEDSIGNATURES)
Default elements of a `Progress`: percentage, bar and estimated remaining time."""
default_elements() = (Percentage(), Bar(), ETA())

"""Percentage of completed iterations, e.g. ` 42%`, padded to 3 digits, 100% only once finished."""
struct Percentage <: AbstractProgressElement end
function print_element(::Percentage, p::Progress)
    percentage = p.finished ? 100 : min(99, floor(Int, 100 * p.counter / max(1, p.n)))
    return @sprintf "%3d%%" percentage
end

"""
    Bar(; length = nothing, fill = '━', empty = ' ')

Progress bar of `length` characters, `fill` for the completed and `empty` for the remaining part.
`length = nothing` fills the width of the terminal the other elements leave (shared between several
such bars)."""
struct Bar <: AbstractProgressElement
    length::Union{Int, Nothing}
    fill::Char
    empty::Char
end
Bar(; length::Union{Integer, Nothing} = nothing, fill::Char = '━', empty::Char = ' ') = Bar(length, fill, empty)

isbar(element) = false
isbar(::Bar) = true

# a bar printed on its own (outside of `render_line`) has 20 characters if not given a length
print_element(bar::Bar, p::Progress) = bar_string(bar, p, something(bar.length, 20))

function bar_string(bar::Bar, p::Progress, width::Integer)
    width <= 0 && return ""
    nfill = p.finished ? width : clamp(floor(Int, width * p.counter / max(1, p.n)), 0, width)
    return repeat(bar.fill, nfill) * repeat(bar.empty, width - nfill)
end

"""Estimated remaining time, e.g. `ETA: 0:01:23`, the total time `Time: 0:02:00` once finished."""
struct ETA <: AbstractProgressElement end
function print_element(::ETA, p::Progress)
    p.finished && return "Time: " * durationstring(p.tcurrent - p.tinit)
    return "ETA: " * durationstring(eta(p, p.tcurrent))
end

"""$(TYPEDSIGNATURES)
Estimated remaining time [s] of the progress meter `p` at time `t`, extrapolating the
average time per iteration so far. `NaN` before the first iteration."""
function eta(p::Progress, t::Real = time())
    p.counter <= 0 && return NaN
    elapsed = t - p.tinit
    return elapsed * (p.n - p.counter) / p.counter
end

"""Time elapsed since the progress meter was created, e.g. `0:00:42`."""
struct ElapsedTime <: AbstractProgressElement end
print_element(::ElapsedTime, p::Progress) = durationstring(p.tcurrent - p.tinit)

"""Number of iterations so far."""
struct Counter <: AbstractProgressElement end
print_element(::Counter, p::Progress) = string(p.counter)

"""$(TYPEDSIGNATURES)
Duration of `seconds` as `h:mm:ss`, with days for 1 to 9 days (`2 days, 3:04:05`) and in days
beyond (`12.50 days`), `N/A` if not finite or negative."""
function durationstring(seconds::Real)
    (isfinite(seconds) && seconds >= 0) || return "N/A"
    seconds_per_day = 24 * 60 * 60
    seconds >= 10 * seconds_per_day && return @sprintf "%.2f days" seconds / seconds_per_day
    s = floor(Int, seconds)
    days, s = divrem(s, seconds_per_day)
    hours, s = divrem(s, 60 * 60)
    minutes, s = divrem(s, 60)
    hhmmss = @sprintf "%d:%02d:%02d" hours minutes s
    return days > 0 ? "$days days, $hhmmss" : hhmmss
end

end # module ProgressMeters
