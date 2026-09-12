export Radiation

"""Radiation scheme bundling a `shortwave` and a `longwave` scheme into one model component
`model.radiation`. Inside the fused column kernel the shortwave scheme is called first, then the
longwave scheme, exactly as the two previously separate components were. Either can be `nothing`
to disable that stream. A radiation scheme that computes both streams at once (e.g. a
correlated-k scheme sharing gas optics) subtypes `AbstractRadiation` directly instead, declares
its variables (e.g. `variables(::MyRadiation) = (variables(OneBandShortwave(SG))..., variables(OneBandLongwave(SG))...)`
for the standard diagnostics) and is passed as `radiation = MyRadiation(...)`. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct Radiation{SW, LW} <: AbstractRadiation
    "[OPTION] Shortwave radiation scheme, `<: AbstractShortwave` or `nothing`"
    @component shortwave::SW

    "[OPTION] Longwave radiation scheme, `<: AbstractLongwave` or `nothing`"
    @component longwave::LW
end

Adapt.@adapt_structure Radiation

"""$(TYPEDSIGNATURES) Bundle `shortwave` and `longwave` into one `Radiation` component,
defaulting to `OneBandShortwave` and `OneBandLongwave` (the `PrimitiveWetModel` defaults).
For a dry model use `shortwave = OneBandGreyShortwave(spectral_grid)` and
`longwave = OneBandGreyLongwave(spectral_grid)`."""
function Radiation(
        SG::SpectralGrid;
        shortwave = OneBandShortwave(SG),
        longwave = OneBandLongwave(SG),
    )
    return Radiation(shortwave, longwave)
end

Base.show(io::IO, radiation::Radiation) = show(io, radiation, values = false)

function initialize!(radiation::Radiation, model::PrimitiveEquation)
    initialize!(radiation.shortwave, model)
    initialize!(radiation.longwave, model)
    return nothing
end

# variables of both sub-schemes; duplicates are removed by identifier in filter_variables
variables(radiation::Radiation) = (variables(radiation.shortwave)..., variables(radiation.longwave)...)
variables(radiation::Radiation, model::AbstractModel) =
    (variables(radiation.shortwave, model)..., variables(radiation.longwave, model)...)

@propagate_inbounds function parameterization!(ij, vars, radiation::Radiation, model)
    parameterization!(ij, vars, radiation.shortwave, model)     # shortwave first,
    parameterization!(ij, vars, radiation.longwave, model)      # then longwave, as before
    return nothing
end

# DEPRECATION of the shortwave_radiation/longwave_radiation model fields (v0.23)
"""$(TYPEDSIGNATURES)
Translate the deprecated model keywords `shortwave_radiation` and `longwave_radiation` into a
single `radiation = Radiation(shortwave, longwave)` (missing stream taken from `default`), and
replace `:shortwave_radiation`, `:longwave_radiation` in a user-provided `parameterizations` tuple
by `:radiation`. Returns the keyword arguments unchanged if none of these are present."""
function deprecated_radiation_kwargs(default::Radiation, kwargs)
    has_shortwave = haskey(kwargs, :shortwave_radiation)
    has_longwave = haskey(kwargs, :longwave_radiation)
    has_old_parameterizations = haskey(kwargs, :parameterizations) &&
        any(p -> p in (:shortwave_radiation, :longwave_radiation), kwargs[:parameterizations])
    (has_shortwave || has_longwave || has_old_parameterizations) || return kwargs

    new_kwargs = Dict{Symbol, Any}(pairs(kwargs))

    if has_shortwave || has_longwave
        haskey(kwargs, :radiation) && throw(
            ArgumentError(
                "Pass either `radiation = Radiation(shortwave, longwave)` or the deprecated " *
                    "`shortwave_radiation`/`longwave_radiation` keywords, not both."
            )
        )
        shortwave = pop!(new_kwargs, :shortwave_radiation, default.shortwave)
        longwave = pop!(new_kwargs, :longwave_radiation, default.longwave)
        new_kwargs[:radiation] = Radiation(shortwave, longwave)
        Base.depwarn(
            "The model keywords `shortwave_radiation` and `longwave_radiation` are deprecated, " *
                "use `radiation = Radiation(spectral_grid; shortwave, longwave)` instead.",
            :deprecated_radiation_kwargs, force = true
        )
    end

    if has_old_parameterizations
        old = Tuple(kwargs[:parameterizations])
        replaced = false
        new = Symbol[]
        for p in old
            if p in (:shortwave_radiation, :longwave_radiation)
                replaced || push!(new, :radiation)      # :radiation at the position of the first
                replaced = true
            else
                push!(new, p)
            end
        end
        new_kwargs[:parameterizations] = Tuple(new)
        Base.depwarn(
            "`:shortwave_radiation` and `:longwave_radiation` in the `parameterizations` tuple are " *
                "deprecated, use a single `:radiation` instead.",
            :deprecated_radiation_kwargs, force = true
        )
    end

    return new_kwargs
end
