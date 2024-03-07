abstract type AbstractAtmosphere <: AbstractModelComponent end
export EarthAtmosphere

"""
$(TYPEDSIGNATURES)
Create a struct `EarthAtmosphere <: AbstractAtmosphere`, with the following physical/chemical
characteristics. Keyword arguments are
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct EarthAtmosphere{NF<:AbstractFloat} <: AbstractAtmosphere
    "molar mass of dry air [g/mol]"
    mol_mass_dry_air::NF = 28.9649

    "molar mass of water vapour [g/mol]"
    mol_mass_vapour::NF = 18.0153

    "specific heat at constant pressure cₚ [J/K/kg]"
    heat_capacity::NF = 1004

    "universal gas constant [J/K/mol]"
    R_gas::NF = 8.3145

    "specific gas constant for dry air [J/kg/K]"
    R_dry::NF = 1000*R_gas/mol_mass_dry_air

    "specific gas constant for water vapour [J/kg/K]"
    R_vapour::NF = 1000*R_gas/mol_mass_vapour

    "Ratio of gas constants: dry air / water vapour, often called ε [1]"
    mol_ratio::NF = R_dry/R_vapour

    "Virtual temperature Tᵥ calculation, Tᵥ = T(1 + μ*q), humidity q, absolute tempereature T"
    μ_virt_temp::NF = (1-mol_ratio)/mol_ratio

    "= R_dry/cₚ, gas const for air over heat capacity"
    κ::NF = R_dry/heat_capacity

    "water density [kg/m³]"
    water_density::NF = 1000

    "latent heat of condensation [J/kg] for consistency with specific humidity [kg/kg]"
    latent_heat_condensation::NF = 2501e3

    "latent heat of sublimation [J/kg]"
    latent_heat_sublimation::NF = 2801e3

    "stefan-Boltzmann constant [W/m²/K⁴]"
    stefan_boltzmann::NF = 5.67e-8

    "surface reference pressure [Pa]"
    pres_ref::NF = 1e5

    "surface reference temperature [K]"
    temp_ref::NF = 288

    "reference moist-adiabatic temperature lapse rate [K/m]"
    lapse_rate::NF = 5/1000

    "layer thickness for the shallow water model [m]"
    layer_thickness::NF = 8500
end

EarthAtmosphere(SG::SpectralGrid; kwargs...) = EarthAtmosphere{SG.NF}(; kwargs...)
EarthAtmosphere(::Type{NF}; kwargs...) where NF = EarthAtmosphere{NF}(; kwargs...)
