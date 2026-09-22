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

variables(radiation::Radiation, model::AbstractModel) =
    (variables(radiation.shortwave, model)..., variables(radiation.longwave, model)...)

@propagate_inbounds function parameterization!(ij, vars, radiation::Radiation, model)
    parameterization!(ij, vars, radiation.shortwave, model)     # shortwave first,
    parameterization!(ij, vars, radiation.longwave, model)      # then longwave, as before
    return nothing
end
