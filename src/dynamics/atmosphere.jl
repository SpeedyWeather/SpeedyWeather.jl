abstract type AbstractAtmosphere <: AbstractModelComponent end
abstract type AbstractWetAtmosphere <: AbstractAtmosphere end
abstract type AbstractDryAtmosphere <: AbstractAtmosphere end

export EarthAtmosphere

"""
$(TYPEDSIGNATURES)
Create a struct `EarthAtmosphere <: AbstractAtmosphere`, with the following physical/chemical
characteristics. Keyword arguments are
$(TYPEDFIELDS)"""
@parameterized @kwdef struct EarthAtmosphere{NF <: AbstractFloat} <: AbstractWetAtmosphere
    "[OPTION] Molar mass of dry air [g/mol]"
    @param mol_mass_dry_air::NF = 28.9649 (bounds = Positive,)

    "[OPTION] Molar mass of water vapor [g/mol]"
    @param mol_mass_vapor::NF = 18.0153 (bounds = Positive,)

    "[OPTION] Specific heat at constant pressure cₚ [J/K/kg]"
    @param heat_capacity::NF = 1004.0 (bounds = Nonnegative,)

    "[OPTION] Universal gas constant [J/K/mol]"
    R_gas::NF = 8.3145

    "[OPTION] Specific gas constant for dry air [J/kg/K]"
    R_dry::NF = 1000 * R_gas / mol_mass_dry_air

    "[OPTION] Specific gas constant for water vapor [J/kg/K]"
    @param R_vapor::NF = 1000 * R_gas / mol_mass_vapor (bounds = Positive,)

    "[OPTION] Ratio of gas constants: dry air / water vapor, often called ε [1]"
    @param mol_ratio::NF = R_dry / R_vapor (bounds = Positive,)

    "[OPTION] Virtual temperature Tᵥ calculation, Tᵥ = T(1 + μ*q), specific humidity q [kg/kg], absolute tempereature T [K]"
    @param μ_virt_temp::NF = (1 - mol_ratio) / mol_ratio (bounds = Positive,)

    "[OPTION] κ = R_dry/cₚ, gas constant for dry air over heat capacity [1]"
    @param κ::NF = R_dry / heat_capacity

    "[OPTION] Water density [kg/m³]"
    @param water_density::NF = 1000.0 (bounds = Positive,)

    "[OPTION] Latent heat of condensation [J/kg]"
    @param latent_heat_condensation::NF = 2501.0e3 (bounds = Nonnegative,)

    "[OPTION] Latent heat of sublimation [J/kg]"
    @param latent_heat_sublimation::NF = 2801.0e3 (bounds = Nonnegative,)

    "[OPTION] Latent heat of freezing/fusion of ice [J/kg]"
    @param latent_heat_fusion::NF = 3.3e5 (bounds = Nonnegative,)

    "[OPTION] Stefan-Boltzmann constant [W/m²/K⁴]"
    stefan_boltzmann::NF = 5.67e-8

    "[OPTION] Surface reference pressure [Pa]"
    @param pressure_reference::NF = 1.0e5 (bounds = Positive,)

    "[OPTION] Saturation vapor pressure at freezing point (0°C) [Pa]"
    @param saturation_vapor_pressure::NF = 610.78 (bounds = Positive,)

    "[OPTION] Surface reference temperature [K]"
    @param temperature_reference::NF = 288.0 (bounds = Nonnegative,)

    "[OPTION] Temperature of freezing point of water [K]"
    @param temperature_freezing::NF = 273.15 (bounds = Nonnegative,)

    "[OPTION] Reference moist-adiabatic temperature lapse rate [K/m]"
    @param moist_lapse_rate::NF = 5 / 1000

    "[OPTION] Reference dry-adiabatic temperature lapse rate [K/m]"
    @param dry_lapse_rate::NF = 9.8 / 1000
end

Adapt.@adapt_structure EarthAtmosphere
EarthAtmosphere(SG::SpectralGrid; kwargs...) = EarthAtmosphere{SG.NF}(; kwargs...)
EarthAtmosphere(::Type{NF}; kwargs...) where {NF} = EarthAtmosphere{NF}(; kwargs...)
Base.eltype(::EarthAtmosphere{NF}) where {NF} = NF

export EarthDryAtmosphere
@parameterized @kwdef struct EarthDryAtmosphere{NF <: AbstractFloat} <: AbstractDryAtmosphere
    "[OPTION] Molar mass of dry air [g/mol]"
    @param mol_mass_dry_air::NF = 28.9649 (bounds = Positive,)

    "[OPTION] Specific heat at constant pressure cₚ [J/K/kg]"
    @param heat_capacity::NF = 1004.0 (bounds = Nonnegative,)

    "[OPTION] Universal gas constant [J/K/mol]"
    R_gas::NF = 8.3145

    "[OPTION] Specific gas constant for dry air [J/kg/K]"
    R_dry::NF = 1000 * R_gas / mol_mass_dry_air

    "[OPTION] κ = R_dry/cₚ, gas constant for dry air over heat capacity"
    @param κ::NF = R_dry / heat_capacity

    "[OPTION] Stefan-Boltzmann constant [W/m²/K⁴]"
    stefan_boltzmann::NF = 5.67e-8

    "[OPTION] Surface reference pressure [Pa]"
    @param pressure_reference::NF = 1.0e5 (bounds = Positive,)

    "[OPTION] Surface reference temperature [K]"
    @param temperature_reference::NF = 288.0 (bounds = Nonnegative,)

    "[OPTION] Reference dry-adiabatic temperature lapse rate [K/m]"
    @param dry_lapse_rate::NF = 9.8 / 1000

    "[OPTION] Layer thickness for the shallow water model [m]"
    @param layer_thickness::NF = 8500 (bounds = Positive,)
end

Adapt.@adapt_structure EarthDryAtmosphere
EarthDryAtmosphere(SG::SpectralGrid; kwargs...) = EarthDryAtmosphere{SG.NF}(; kwargs...)
EarthDryAtmosphere(::Type{NF}; kwargs...) where {NF} = EarthDryAtmosphere{NF}(; kwargs...)
Base.eltype(::EarthDryAtmosphere{NF}) where {NF} = NF

## FUNCTIONS WITH ATMOSPHERE AS ARGUMENT

lapse_rate(A::AbstractWetAtmosphere) = A.moist_lapse_rate
lapse_rate(A::AbstractDryAtmosphere) = A.dry_lapse_rate

latent_heat_condensation(A::AbstractWetAtmosphere) = A.latent_heat_condensation
latent_heat_condensation(A::AbstractDryAtmosphere) = zero(eltype(A))

latent_heat_sublimation(A::AbstractWetAtmosphere) = A.latent_heat_sublimation
latent_heat_sublimation(A::AbstractDryAtmosphere) = zero(eltype(A))

saturation_vapor_pressure(T, A::AbstractDryAtmosphere) = zero(T)

"""$(TYPEDSIGNATURES)
Saturation water vapor pressure as a function of temperature using the
Clausius-Clapeyron equation, 

    e(T) = e₀ * exp( -Lᵥ/Rᵥ * (1/T - 1/T₀)),

where T is in Kelvin, Lᵥ the the latent heat of condensation and Rᵥ the gas constant of water vapor"""
function saturation_vapor_pressure(T, A::AbstractWetAtmosphere)
    e₀ = A.saturation_vapor_pressure
    Lᵥ = A.latent_heat_condensation
    (; R_vapor) = A
    T₀ = A.temperature_freezing
    return e₀ * exp(Lᵥ/R_vapor*(inv(T₀) - inv(T)))
end

saturation_humidity(T, p, A::AbstractDryAtmosphere) = zero(T)

"""$(TYPEDSIGNATURES)
Saturation specific humidity as a function of temperature [K] and pressure [Pa], defined as
    qₛ = ε * eₛ(T) / p,
with saturation vapor pressure eₛ(T), pressure p [Pa] and ratio of gas constants ε = R_dry/R_vapor."""
function saturation_humidity(T, p, A::AbstractWetAtmosphere)
    e_s = saturation_vapor_pressure(T, A)
    ϵ = A.mol_ratio
    return ϵ * e_s / p
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
