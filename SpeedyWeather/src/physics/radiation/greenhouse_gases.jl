abstract type AbstractGreenhouseGas <: AbstractModelComponent end

DEFAULT_CO2 = 280 # ppm
export CO2
@kwdef struct CO2{F} <: AbstractGreenhouseGas
    concentration::F
end

CO2(SG::SpectralGrid, concentration::Real = DEFAULT_CO2) = CO2{SG.NF}(concentration)
CO2(SG::SpectralGrid, concentration::Function) = CO2(concentration)

(C::CO2)(::DateTime) = C.concentration
(C::CO2{<:Function})(t::DateTime) = C.concentration(t)

export TwoTimesCO2, FourTimesCO2
TwoTimesCO2() = CO2(t -> t > DateTime(2000) ? 2DEFAULT_CO2 : DEFAULT_CO2)
FourTimesCO2() = CO2(t -> t > DateTime(2000) ? 4DEFAULT_CO2 : DEFAULT_CO2)

export ExponentialCO2
@kwdef struct ExponentialCO2{NF} <: AbstractGreenhouseGas
    base_concentration::NF = 280
    start::DateTime = DateTime(1850)
    a::NF = 3.6
    b::NF = 1/47
end

ExponentialCO2(SG::SpectralGrid; kwargs...) = ExponentialCO2{SG.NF}(; kwargs...)
function (C::ExponentialCO2{NF})(t::DateTime) where NF

    # time since start in seconds
    Δt = Dates.datetime2unix(t) - Dates.datetime2unix(C.start)
    Δt /= (365.25f0 * 24 * 3600) # time in years
    
    return C.base_concentration + C.a * exp(C.b * Δt)
end