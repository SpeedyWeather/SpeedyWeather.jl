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

    "latent heat of condensation [J/kg]"
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
    moist_lapse_rate::NF = 5/1000
    
    "reference dry-adiabatic temperature lapse rate [K/m]"
    dry_lapse_rate::NF = 9.8/1000

    "layer thickness for the shallow water model [m]"
    layer_thickness::NF = 8500
end

EarthAtmosphere(SG::SpectralGrid; kwargs...) = EarthAtmosphere{SG.NF}(; kwargs...)
EarthAtmosphere(::Type{NF}; kwargs...) where NF = EarthAtmosphere{NF}(; kwargs...)
Base.eltype(::EarthAtmosphere{NF}) where NF = NF

parameters(atm::EarthAtmosphere; kwargs...) =
    (
        mol_mass_dry_air = parameters(atm.mol_mass_dry_air; desc="molar mass of dry air [g/mol]", kwargs...),
        mol_mass_vapour = parameters(atm.mol_mass_vapour; desc="molar mass of water vapour [g/mol]", kwargs...),
        heat_capacity = parameters(atm.heat_capacity; desc="specific heat at constant pressure cₚ [J/K/kg]", kwargs...),
        R_gas = parameters(atm.R_gas; desc="universal gas constant [J/K/mol]", kwargs...),
        R_dry = parameters(atm.R_dry; desc="specific gas constant for dry air [J/kg/K]", kwargs...),
        R_vapour = parameters(atm.R_vapour; desc="specific gas constant for water vapour [J/kg/K]", kwargs...),
        mol_ratio = parameters(atm.mol_ratio; desc="Ratio of gas constants: dry air / water vapour, often called ε [1]", kwargs...),
        μ_virt_temp = parameters(atm.μ_virt_temp; desc="Virtual temperature Tᵥ calculation factor", kwargs...),
        κ = parameters(atm.κ; desc="= R_dry/cₚ, gas const for air over heat capacity", kwargs...),
        water_density = parameters(atm.water_density; desc="water density [kg/m³]", kwargs...),
        latent_heat_condensation = parameters(atm.latent_heat_condensation; desc="latent heat of condensation [J/kg]", kwargs...),
        latent_heat_sublimation = parameters(atm.latent_heat_sublimation; desc="latent heat of sublimation [J/kg]", kwargs...),
        stefan_boltzmann = parameters(atm.stefan_boltzmann; desc="stefan-Boltzmann constant [W/m²/K⁴]", kwargs...),
        pres_ref = parameters(atm.pres_ref; desc="surface reference pressure [Pa]", kwargs...),
        temp_ref = parameters(atm.temp_ref; desc="surface reference temperature [K]", kwargs...),
        moist_lapse_rate = parameters(atm.moist_lapse_rate; desc="reference moist-adiabatic temperature lapse rate [K/m]", kwargs...),
        dry_lapse_rate = parameters(atm.dry_lapse_rate; desc="reference dry-adiabatic temperature lapse rate [K/m]", kwargs...),
        layer_thickness = parameters(atm.layer_thickness; desc="layer thickness for the shallow water model [m]", kwargs...),
    )
