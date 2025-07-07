abstract type AbstractAtmosphere <: AbstractModelComponent end
export EarthAtmosphere

"""
$(TYPEDSIGNATURES)
Create a struct `EarthAtmosphere <: AbstractAtmosphere`, with the following physical/chemical
characteristics. Keyword arguments are
$(TYPEDFIELDS)"""
@parameterized Base.@kwdef mutable struct EarthAtmosphere{NF<:AbstractFloat} <: AbstractAtmosphere
    "molar mass of dry air [g/mol]"
    @param mol_mass_dry_air::NF = 28.9649 (bounds=Positive,)

    "molar mass of water vapour [g/mol]"
    @param mol_mass_vapour::NF = 18.0153 (bounds=Positive,)

    "specific heat at constant pressure cₚ [J/K/kg]"
    @param heat_capacity::NF = 1004 (bounds=Nonnegative,)

    "universal gas constant [J/K/mol]"
    R_gas::NF = 8.3145

    "specific gas constant for dry air [J/kg/K]"
    @param R_dry::NF = 1000*R_gas/mol_mass_dry_air (bounds=Positive,)

    "specific gas constant for water vapour [J/kg/K]"
    @param R_vapour::NF = 1000*R_gas/mol_mass_vapour (bounds=Positive,)

    "Ratio of gas constants: dry air / water vapour, often called ε [1]"
    @param mol_ratio::NF = R_dry/R_vapour (bounds=Positive,)

    "Virtual temperature Tᵥ calculation, Tᵥ = T(1 + μ*q), humidity q, absolute tempereature T"
    @param μ_virt_temp::NF = (1-mol_ratio)/mol_ratio (bounds=Positive,)

    "= R_dry/cₚ, gas const for air over heat capacity"
    @param κ::NF = R_dry/heat_capacity

    "water density [kg/m³]"
    @param water_density::NF = 1000 (bounds=Positive,)

    "latent heat of condensation [J/kg]"
    @param latent_heat_condensation::NF = 2501e3 (bounds=Nonnegative,)

    "latent heat of sublimation [J/kg]"
    @param latent_heat_sublimation::NF = 2801e3 (bounds=Nonnegative,)

    "stefan-Boltzmann constant [W/m²/K⁴]"
    stefan_boltzmann::NF = 5.67e-8

    "surface reference pressure [Pa]"
    @param pres_ref::NF = 1e5 (bounds=Positive,)

    "surface reference temperature [K]"
    @param temp_ref::NF = 288 (bounds=Nonnegative,)

    "reference moist-adiabatic temperature lapse rate [K/m]"
    @param moist_lapse_rate::NF = 5/1000
    
    "reference dry-adiabatic temperature lapse rate [K/m]"
    @param dry_lapse_rate::NF = 9.8/1000

    "layer thickness for the shallow water model [m]"
    @param layer_thickness::NF = 8500 (bounds=Positive,)
end

EarthAtmosphere(SG::SpectralGrid; kwargs...) = EarthAtmosphere{SG.NF}(; kwargs...)
EarthAtmosphere(::Type{NF}; kwargs...) where NF = EarthAtmosphere{NF}(; kwargs...)
Base.eltype(::EarthAtmosphere{NF}) where NF = NF
