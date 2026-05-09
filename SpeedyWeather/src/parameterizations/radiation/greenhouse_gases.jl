abstract type AbstractGreenhouseGas <: AbstractModelComponent end
abstract type AbstractCO2 <: AbstractGreenhouseGas end

# greenhouse gases are implemented as named tuple so to identify them by symbol
# also pass on the name as the key of the named tuple element
variables(name_ghg::Pair{<:Symbol, <:AbstractGreenhouseGas}, model::AbstractModel) =
    variables(name_ghg.second, name_ghg.first, model)

variables(::AbstractCO2, name::Symbol, args...) = (
    PrognosticVariable(name, ScalarDim(), units="ppm", desc="Carbon dioxide", namespace=:greenhouse_gases),
)

"""$(TYPEDSIGNATURES)
Time step greenhouse gas concentration as function of time only."""
function greenhouse_gases_time_step!(vars::Variables, model::PrimitiveEquation)
    isnothing(model.greenhouse_gases) && return nothing     # if greenhouse_gases is nothing, then skip this step
    t = vars.prognostic.clock.time
    for gas in keys(model.greenhouse_gases)
        # greenhouse gases are functions of time only, overwrite current value
        vars.prognostic.greenhouse_gases[gas][] = model.greenhouse_gases[gas](t)
    end
end

# fallback: for models without greenhouse gases, do nothing
greenhouse_gases_time_step!(vars::Variables, model::AbstractModel) = nothing   

"""$(TYPEDSIGNATURES) CO2 uses unit of ppm, so *1e-6 to convert to kg/kg."""
@inline unit(::AbstractCO2) = 1f-6

"""$(TYPEDSIGNATURES) The molar mass ratio CO2 / molar mass of dry air"""
@inline molar_mass_ratio(::AbstractCO2) = 44f0 / 29f0

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

@inline (C::CO2)(::DateTime) = C.concentration
@inline (C::CO2{<:Function})(t::DateTime) = C.concentration(t)

export TwoTimesCO2, FourTimesCO2
"$(TYPEDSIGNATURES) 2xCO2 scenario: Doubles preindustrial `DEFAULT_CO2` in the year 2000."
TwoTimesCO2(t0 = DateTime(2000)) = CO2(t -> t > t0 ? 2DEFAULT_CO2 : DEFAULT_CO2)

"$(TYPEDSIGNATURES) 4xCO2 scenario: Quadruples preindustrial `DEFAULT_CO2` in the year 2000."
FourTimesCO2(t0 = DateTime(2000)) = CO2(t -> t > t0 ? 4DEFAULT_CO2 : DEFAULT_CO2)

export ExponentialCO2

"""$(TYPEDSIGNATURES) Exponential CO2 concentration with defaults fitted
to the Keeling curve, ignoring seasonal variation."""
@kwdef struct ExponentialCO2{NF, DT} <: AbstractCO2
    "[OPTION] preindustrial CO2 concentration [ppm]"
    base_concentration::NF = 280

    "[OPTION] time of preindustrial CO2"
    start::DT = DateTime(1850)

    "[OPTION] Scaling `a` in `c + a*exp((t-t0)*b)` [ppm]"
    a::NF = 3.85

    "[OPTION] e-folding frequency b in `c + a*exp((t-t0)*b)` [1/year]"
    b::NF = 1 / 48
end

ExponentialCO2(SG::SpectralGrid; kwargs...) = ExponentialCO2{SG.NF, DateTime}(; kwargs...)
function (C::ExponentialCO2{NF})(t) where {NF}

    # time since start in seconds, convert from Float64 to NF
    Δt = convert(NF, Dates.datetime2unix(t) - Dates.datetime2unix(C.start))
    Δt *= convert(NF, 1 / (365.25f0 * 24 * 3600))    # time in years
    return C.base_concentration + C.a * exp(C.b * Δt)
end