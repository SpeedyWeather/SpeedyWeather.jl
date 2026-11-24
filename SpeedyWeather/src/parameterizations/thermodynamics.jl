abstract type AbstractSaturationVaporPressure end

export ClausiusClapeyronConstantL, ClausiusClapeyronLinearL, IFSMixedPhase, MagnusBolton

# Interface: every formulation implements
#   saturation_vapor_pressure(svp, T, atmosphere)  -> eₛ(T) [Pa], T in Kelvin
#   grad(svp, T, atmosphere)                       -> deₛ/dT [Pa/K]
# The atmosphere argument lets a formulation use constants consistent with the rest
# of the model (e.g. latent heat, gas constant), formulations with own constants ignore it.
# All functions are branch-free (ifelse) so that they can be traced by Reactant.

"""$(TYPEDSIGNATURES)
Generator function for saturation vapor pressure formulations, converting
to the number format `NF`, e.g. `ClausiusClapeyronConstantL(Float32)` or
`ClausiusClapeyronConstantL(spectral_grid)`."""
function (::Type{SVP})(::Type{NF}; kwargs...) where {SVP <: AbstractSaturationVaporPressure, NF}
    return SVP{NF}(; kwargs...)
end
(::Type{SVP})(SG::SpectralGrid; kwargs...) where {SVP <: AbstractSaturationVaporPressure} = SVP(SG.NF; kwargs...)

"""
Saturation vapor pressure over liquid water from the Clausius-Clapeyron equation
integrated with a constant latent heat of condensation Lᵥ (and ideal gas, condensed volume neglected),

    eₛ(T) = e₀ exp(Lᵥ/Rᵥ (1/T₀ - 1/T))

with Lᵥ = `atmosphere.latent_heat_condensation` and Rᵥ = `atmosphere.R_vapor` taken from
the atmosphere. Overestimates eₛ by ≈2.4% at 30˚C and ≈4% at 40˚C as Lᵥ decreases with
temperature. Fields are
$(TYPEDFIELDS)"""
@kwdef struct ClausiusClapeyronConstantL{NF} <: AbstractSaturationVaporPressure
    "Saturation vapor pressure at T₀ [Pa]"
    e₀::NF = 610.78

    "Reference temperature [K]"
    T₀::NF = 273.15
end

@inline function saturation_vapor_pressure(svp::ClausiusClapeyronConstantL, T, atmosphere)
    (; e₀, T₀) = svp
    Lᵥ = atmosphere.latent_heat_condensation
    Rᵥ = atmosphere.R_vapor
    T = convert(typeof(T₀), T)
    return e₀ * exp(Lᵥ / Rᵥ * (inv(T₀) - inv(T)))
end

@inline function grad(svp::ClausiusClapeyronConstantL, T, atmosphere)
    T = convert(typeof(svp.T₀), T)
    eₛ = saturation_vapor_pressure(svp, T, atmosphere)
    return eₛ * atmosphere.latent_heat_condensation / (atmosphere.R_vapor * T^2)
end

"""
Saturation vapor pressure over liquid water from the Clausius-Clapeyron equation
integrated with a latent heat that varies linearly with temperature (constant heat capacities),
following Romps (2017, J. Atmos. Sci., eq. 6 therein),

    eₛ(T) = e₀ (T/T₀)^a exp(b (1/T₀ - 1/T)),
    a = (c_pv - c_vl)/Rᵥ,   b = (E₀ᵥ - (c_vv - c_vl) T₀)/Rᵥ

with the reference point at the triple point. Only assumes an ideal gas and neglects
the condensed volume, so it is thermodynamically exact for the given constants
(accurate to ≈0.1% between -40˚C and 50˚C). Uses its own constants (not the atmosphere's)
for consistency with Romps (2017). Fields are
$(TYPEDFIELDS)"""
@kwdef struct ClausiusClapeyronLinearL{NF} <: AbstractSaturationVaporPressure
    "Triple point pressure [Pa]"
    e₀::NF = 611.65

    "Triple point temperature [K]"
    T₀::NF = 273.16

    "Difference in specific internal energy between vapor and liquid at the triple point [J/kg]"
    E₀ᵥ::NF = 2.374e6

    "Gas constant of water vapor [J/kg/K]"
    Rᵥ::NF = 461.0

    "Heat capacity of water vapor at constant volume [J/kg/K]"
    c_vv::NF = 1418.0

    "Heat capacity of liquid water at constant volume [J/kg/K]"
    c_vl::NF = 4119.0
end

@inline function saturation_vapor_pressure(svp::ClausiusClapeyronLinearL, T, atmosphere)
    (; e₀, T₀, E₀ᵥ, Rᵥ, c_vv, c_vl) = svp
    T = convert(typeof(T₀), T)
    a = (c_vv + Rᵥ - c_vl) / Rᵥ          # c_pv = c_vv + Rᵥ
    b = (E₀ᵥ - (c_vv - c_vl) * T₀) / Rᵥ
    return e₀ * (T / T₀)^a * exp(b * (inv(T₀) - inv(T)))
end

@inline function grad(svp::ClausiusClapeyronLinearL, T, atmosphere)
    (; T₀, E₀ᵥ, Rᵥ, c_vv, c_vl) = svp
    T = convert(typeof(T₀), T)
    a = (c_vv + Rᵥ - c_vl) / Rᵥ
    b = (E₀ᵥ - (c_vv - c_vl) * T₀) / Rᵥ
    eₛ = saturation_vapor_pressure(svp, T, atmosphere)
    return eₛ * (a / T + b / T^2)       # d ln eₛ/dT = a/T + b/T² = L(T)/(Rᵥ T²)
end

"""
Saturation vapor pressure as in the ECMWF IFS (IFS documentation, Part IV: Physical processes), Tetens' formula
over liquid water (subscript l) and ice (subscript i)

    eₛ,ₓ(T) = e₀ exp(aₓ (T - T₀)/(T - bₓ)),

blended in a mixed-phase range T_ice < T < T₀ with a liquid fraction α = ((T - T_ice)/(T₀ - T_ice))²

    eₛ(T) = α eₛ,ₗ(T) + (1 - α) eₛ,ᵢ(T).

Pure ice below T_ice (-23˚C), pure liquid above T₀ (0˚C), continuous in between. Fields are
$(TYPEDFIELDS)"""
@kwdef struct IFSMixedPhase{NF} <: AbstractSaturationVaporPressure
    "Saturation vapor pressure at T₀ [Pa]"
    e₀::NF = 611.21

    "Triple point temperature, upper end of mixed phase [K]"
    T₀::NF = 273.16

    "Lower end of mixed phase, pure ice below [K]"
    T_ice::NF = 250.16

    "Tetens numerator over liquid water [1]"
    a_liquid::NF = 17.502

    "Tetens denominator temperature over liquid water [K]"
    b_liquid::NF = 32.19

    "Tetens numerator over ice [1]"
    a_ice::NF = 22.587

    "Tetens denominator temperature over ice [K]"
    b_ice::NF = -0.7
end

# Tetens formula and its derivative, shared by liquid and ice
@inline tetens(e₀, T₀, a, b, T) = e₀ * exp(a * (T - T₀) / (T - b))
@inline grad_tetens(e₀, T₀, a, b, T) = tetens(e₀, T₀, a, b, T) * a * (T₀ - b) / (T - b)^2

# liquid fraction α ∈ [0, 1] of the mixed phase, quadratic in temperature
@inline function liquid_fraction(svp::IFSMixedPhase, T)
    (; T₀, T_ice) = svp
    return clamp((T - T_ice) / (T₀ - T_ice), 0, 1)^2
end

@inline function saturation_vapor_pressure(svp::IFSMixedPhase, T, atmosphere)
    (; e₀, T₀, a_liquid, b_liquid, a_ice, b_ice) = svp
    T = convert(typeof(T₀), T)
    α = liquid_fraction(svp, T)
    e_liquid = tetens(e₀, T₀, a_liquid, b_liquid, T)
    e_ice = tetens(e₀, T₀, a_ice, b_ice, T)
    return α * e_liquid + (1 - α) * e_ice
end

@inline function grad(svp::IFSMixedPhase, T, atmosphere)
    (; e₀, T₀, T_ice, a_liquid, b_liquid, a_ice, b_ice) = svp
    T = convert(typeof(T₀), T)
    α = liquid_fraction(svp, T)
    mixed = (T > T_ice) & (T < T₀)
    dα_dT = ifelse(mixed, 2 * (T - T_ice) / (T₀ - T_ice)^2, zero(T))
    e_liquid = tetens(e₀, T₀, a_liquid, b_liquid, T)
    e_ice = tetens(e₀, T₀, a_ice, b_ice, T)
    de_liquid = grad_tetens(e₀, T₀, a_liquid, b_liquid, T)
    de_ice = grad_tetens(e₀, T₀, a_ice, b_ice, T)
    return α * de_liquid + (1 - α) * de_ice + dα_dT * (e_liquid - e_ice)
end

"""
Saturation vapor pressure over liquid water following the Magnus-type formula of Bolton (1980),

    eₛ(T) = e₀ exp(a (T - T₀)/(T - T₀ + b)),

i.e. with temperature in ˚C in the exponent. Accurate to ≈0.1% between -35˚C and 35˚C. Fields are
$(TYPEDFIELDS)"""
@kwdef struct MagnusBolton{NF} <: AbstractSaturationVaporPressure
    "Saturation vapor pressure at T₀ [Pa]"
    e₀::NF = 611.2

    "0˚C in Kelvin [K]"
    T₀::NF = 273.15

    "Magnus numerator [1]"
    a::NF = 17.67

    "Magnus denominator [˚C]"
    b::NF = 243.5
end

@inline function saturation_vapor_pressure(svp::MagnusBolton, T, atmosphere)
    (; e₀, T₀, a, b) = svp
    T = convert(typeof(T₀), T)
    return tetens(e₀, T₀, a, T₀ - b, T)
end

@inline function grad(svp::MagnusBolton, T, atmosphere)
    (; e₀, T₀, a, b) = svp
    T = convert(typeof(T₀), T)
    return grad_tetens(e₀, T₀, a, T₀ - b, T)
end
