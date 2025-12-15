abstract type AbstractClausiusClapeyron <: AbstractParameterization end

export ClausiusClapeyron

"""
Parameters for computing saturation vapor pressure of water using the Clausis-Clapeyron equation,

    e(T) = e₀ * exp( -Lᵥ/Rᵥ * (1/T - 1/T₀)),

where T is in Kelvin, Lᵥ the the latent heat of condensation and Rᵥ the gas constant of water vapor
$(TYPEDFIELDS)"""
@kwdef struct ClausiusClapeyron{NF<:AbstractFloat} <: AbstractClausiusClapeyron
    "Saturation water vapor pressure at 0°C [Pa]"
    e₀::NF = 610.78

    "0°C in Kelvin"
    T₀::NF = 273.15

    "Latent heat of condensation/vaporization of water [J/kg]"
    latent_heat_condensation::NF = 2.5e6

    "Latent heat of freezing/fusion of ice [J/kg]"
    latent_heat_fusion::NF = 3.3e5

    "Latent heat of freezing/sublimation of ice [J/kg]"
    latent_heat_sublimation::NF = 2.8e6

    "Specific heat of air at constant pressure [J/K/kg]"
    heat_capacity::NF = 1004.64

    "Gas constant of water vapor [J/kg/K]"
    R_vapor::NF = 461.5

    "Gas constant of dry air [J/kg/K]"
    R_dry::NF = 287.04

    "[DERIVED] Latent heat of condensation divided by gas constant of water vapor [K]"    
    Lᵥ_Rᵥ::NF = latent_heat_condensation/R_vapor

    "[DERIVED] Inverse of T₀, one over 0°C in Kelvin"
    T₀⁻¹::NF = inv(T₀)

    "[DERIVED] Ratio of molecular masses [1] of water vapor over dry air (=R_dry/R_vapor)."
    mol_ratio::NF = R_dry/R_vapor
end

Adapt.@adapt_structure ClausiusClapeyron

# generator function
function ClausiusClapeyron(SG::SpectralGrid, atm::AbstractAtmosphere; kwargs...)
    (; R_dry, R_vapor, latent_heat_condensation, heat_capacity) = atm
    return ClausiusClapeyron{SG.NF}(; latent_heat_condensation, R_dry, R_vapor, heat_capacity, kwargs...)
end

"""
$(TYPEDSIGNATURES)
Saturation water vapor pressure as a function of temperature using the
Clausius-Clapeyron equation,

    e(T) = e₀ * exp( -Lᵥ/Rᵥ * (1/T - 1/T₀)),

where T is in Kelvin, Lᵥ the the latent heat of vaporization and Rᵥ the gas constant
of water vapor, T₀ is 0˚C in Kelvin."""
@inline function saturation_vapor_pressure(CC::ClausiusClapeyron{NF}, temp_kelvin::NF) where NF
    (; e₀, T₀⁻¹, Lᵥ_Rᵥ) = CC
    return e₀ * exp(Lᵥ_Rᵥ*(T₀⁻¹ - inv(temp_kelvin)))
end

# convert to number format of struct
@inline function saturation_vapor_pressure(CC::ClausiusClapeyron{NF}, temp_kelvin) where NF
    saturation_vapor_pressure(CC, convert(NF, temp_kelvin))
end

# Keep functor for backwards compatibility
@inline (CC::ClausiusClapeyron)(temp_kelvin) = saturation_vapor_pressure(CC, temp_kelvin)

"""
$(TYPEDSIGNATURES)
Gradient of Clausius-Clapeyron wrt to temperature, evaluated at `temp_kelvin`."""
@inline function grad(CC::ClausiusClapeyron{NF}, temp_kelvin::NF) where NF
    e = saturation_vapor_pressure(CC, temp_kelvin)
    return e*CC.Lᵥ_Rᵥ/temp_kelvin^2
end

# convert to input argument to number format from struct
@inline grad(CC::ClausiusClapeyron{NF}, temp_kelvin) where NF = grad(CC, convert(NF, temp_kelvin))

"""
$(TYPEDSIGNATURES)
Saturation humidity from saturation vapor pressure and pressure via

    qsat = mol_ratio*sat_vap_pres/pres

with both pressures in same units and qsat in kg/kg."""
@inline function saturation_humidity(
    sat_vap_pres::NF,                   # saturation vapor pressure [Pa]
    pres::NF;                           # pressure [Pa]
    mol_ratio::NF = NF(287.04/461.5)    # ratio of mol masses dry air - water vapor
) where NF
    # simpler version given that pres >> sat_vap_pres for most temperatures T<40˚C
    return mol_ratio*sat_vap_pres / pres
    # more accurate version, with a more complicated derivative though
    # return mol_ratio*sat_vap_pres / (pres - (1-mol_ratio)*sat_vap_pres)
end

"""
$(TYPEDSIGNATURES)
Saturation humidity [kg/kg] from temperature [K], pressure [Pa] via

    sat_vap_pres = saturation_vapor_pressure(clausius_clapeyron, temperature)
    saturation humidity = mol_ratio * sat_vap_pres / pressure"""
@inline function saturation_humidity(
    temp_kelvin::NF,
    pres::NF,
    clausius_clapeyron::AbstractClausiusClapeyron,
) where NF
    sat_vap_pres = saturation_vapor_pressure(clausius_clapeyron, temp_kelvin)
    return saturation_humidity(sat_vap_pres, pres; mol_ratio=clausius_clapeyron.mol_ratio)
end

"""
$(TYPEDSIGNATURES)
Gradient of Clausius-Clapeyron wrt to temperature, evaluated at `temp_kelvin`."""
@inline function grad_saturation_humidity(CC::ClausiusClapeyron{NF}, temp_kelvin::NF, pres::NF) where NF
    qsat = saturation_humidity(temp_kelvin, pres, CC)
    return qsat*CC.Lᵥ_Rᵥ/temp_kelvin^2
end

export TetensEquation

"""
Parameters for computing saturation vapor pressure of water using the Tetens' equation,

    eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T + Tᵢ)),

where T is in Kelvin and i = 1, 2 for saturation above and below freezing,
respectively. From Tetens (1930), and Murray (1967) for below freezing.
$(TYPEDFIELDS)"""
@kwdef struct TetensEquation{NF<:AbstractFloat} <: AbstractClausiusClapeyron
    "Saturation water vapor pressure at 0°C [Pa]"
    e₀::NF = 610.78

    "0°C in Kelvin"
    T₀::NF = 273.15

    "Tetens denominator (water) [˚C]"
    T₁::NF = 237.3

    "Tetens denominator following Murray (1967, below freezing) [˚C]"
    T₂::NF = 265.5

    "Tetens numerator scaling [1], above freezing"
    C₁::NF = 17.27

    "Tetens numerator scaling [1], below freezing"
    C₂::NF = 21.875
end

"""
$(TYPEDSIGNATURES)
Saturation water vapor pressure as a function of temperature using the
Tetens equation,

    eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T - Tᵢ)),

where T is in Kelvin and i = 1, 2 for saturation above and below freezing,
respectively."""
@inline function saturation_vapor_pressure(TetensCoefficients::TetensEquation{NF}, temp_kelvin::NF) where NF
    (; e₀, T₀, C₁, C₂, T₁, T₂) = TetensCoefficients
    C, T = temp_kelvin > T₀ ? (C₁, T₁) : (C₂, T₂)      # change coefficients above/below freezing
    temp_celsius = temp_kelvin - T₀
    return e₀ * exp(C * temp_celsius / (temp_celsius + T))
end

# Keep functor for backwards compatibility
@inline (TE::TetensEquation)(temp_kelvin) = saturation_vapor_pressure(TE, temp_kelvin)

"""
$(TYPEDSIGNATURES)
Gradient of the Tetens equation wrt to temperature, evaluated at `temp_kelvin`."""
@inline function grad(TetensCoefficients::TetensEquation{NF}, temp_kelvin::NF) where NF
    (; T₀, C₁, C₂, T₁, T₂) = TetensCoefficients
    e = saturation_vapor_pressure(TetensCoefficients, temp_kelvin)  # saturation vapor pressure
    C, T = temp_kelvin > T₀ ? (C₁, T₁) : (C₂, T₂)   # change coefficients above/below freezing
    return e*C*T/(temp_kelvin - T₀ - T)^2           # chain rule: times derivative of inner function
end
