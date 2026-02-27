abstract type AbstractGreenhouseGas <: AbstractModelComponent end
abstract type AbstractCO2 <: AbstractGreenhouseGas end
# variables(::AbstractCO2) = PrognosticVariable(:co2, ...)

const DEFAULT_CO2 = 280     # preindustrial CO2 [ppm]

export CO2
"""$(TYPEDSIGNATURES) CO2 concentration [ppm] that can be a constant, e.g. `CO2(420)`
or a function of `t::DateTime`, e.g. `CO2(t -> 280 + year(t) - 1850)` for a linear
(but stepwise) increase by 1ppm per year."""
@kwdef struct CO2{F} <: AbstractCO2
    "[OPTION] CO2 concentration"
    concentration::F
end

CO2(SG::SpectralGrid, concentration::Real = DEFAULT_CO2) = CO2{SG.NF}(concentration)
CO2(SG::SpectralGrid, concentration::Function) = CO2(concentration)

(C::CO2)(::DateTime) = C.concentration
(C::CO2{<:Function})(t::DateTime) = C.concentration(t)

export TwoTimesCO2, FourTimesCO2
"$(TYPEDSIGNATURES) 2xCO2 scenario: Doubles preindustrial `DEFAULT_CO2` in the year 2000."
TwoTimesCO2(t0 = DateTime(2000)) = CO2(t -> t > t0 ? 2DEFAULT_CO2 : DEFAULT_CO2)

"$(TYPEDSIGNATURES) 4xCO2 scenario: Quadruples preindustrial `DEFAULT_CO2` in the year 2000."
FourTimesCO2(t0 = DateTime(2000)) = CO2(t -> t > t0 ? 4DEFAULT_CO2 : DEFAULT_CO2)

export ExponentialCO2

"""$(TYPEDSIGNATURES) Exponential CO2 concentration with defaults fitted
to the Keeling curve, ignoring seasonal variation."""
@kwdef struct ExponentialCO2{NF} <: AbstractCO2
    "[OPTION] preindustrial CO2 concentration [ppm]"
    base_concentration::NF = 280

    "[OPTION] time of preindustrial CO2"
    start::DateTime = DateTime(1850)

    "[OPTION] Scaling `a` in `c + a*exp((t-t0)*b)` [ppm]"
    a::NF = 3.85

    "[OPTION] e-folding frequency b in `c + a*exp((t-t0)*b)` [1/year]"
    b::NF = 1/48
end

ExponentialCO2(SG::SpectralGrid; kwargs...) = ExponentialCO2{SG.NF}(; kwargs...)
function (C::ExponentialCO2{NF})(t::DateTime) where NF

    # time since start in seconds, convert from Float64 to NF
    Δt = convert(NF, Dates.datetime2unix(t) - Dates.datetime2unix(C.start))
    Δt /= (365.25f0 * 24 * 3600)    # time in years
    Δt = convert(NF, Δt)            # to compute the exponential in NF
    
    return C.base_concentration + C.a * exp(C.b * Δt)
end

# for gas in keys(model.greenhouse_gases)
#     simulation.variables[gas][] = model.greenhouse_gases[gas](time)
# end