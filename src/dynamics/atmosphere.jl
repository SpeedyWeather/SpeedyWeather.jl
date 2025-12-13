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

    "[OPTION] Molar mass of water vapour [g/mol]"
    @param mol_mass_vapour::NF = 18.0153 (bounds = Positive,)

    "[OPTION] Specific heat at constant pressure cₚ [J/K/kg]"
    @param heat_capacity::NF = 1004.0 (bounds = Nonnegative,)

    "[OPTION] Universal gas constant [J/K/mol]"
    R_gas::NF = 8.3145

    "[OPTION] Specific gas constant for dry air [J/kg/K]"
    R_dry::NF = 1000 * R_gas / mol_mass_dry_air

    "[OPTION] Specific gas constant for water vapour [J/kg/K]"
    @param R_vapour::NF = 1000 * R_gas / mol_mass_vapour (bounds = Positive,)

    "[OPTION] Ratio of gas constants: dry air / water vapour, often called ε [1]"
    @param mol_ratio::NF = R_dry / R_vapour (bounds = Positive,)

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

    "[OPTION] Stefan-Boltzmann constant [W/m²/K⁴]"
    stefan_boltzmann::NF = 5.67e-8

    "[OPTION] Surface reference pressure [Pa]"
    @param pres_ref::NF = 1.0e5 (bounds = Positive,)

    "[OPTION] Surface reference temperature [K]"
    @param temp_ref::NF = 288.0 (bounds = Nonnegative,)

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
    @param pres_ref::NF = 1.0e5 (bounds = Positive,)

    "[OPTION] Surface reference temperature [K]"
    @param temp_ref::NF = 288.0 (bounds = Nonnegative,)

    "[OPTION] Reference dry-adiabatic temperature lapse rate [K/m]"
    @param dry_lapse_rate::NF = 9.8 / 1000

    "[OPTION] Layer thickness for the shallow water model [m]"
    @param layer_thickness::NF = 8500 (bounds = Positive,)
end

Adapt.@adapt_structure EarthDryAtmosphere
EarthDryAtmosphere(SG::SpectralGrid; kwargs...) = EarthDryAtmosphere{SG.NF}(; kwargs...)
EarthDryAtmosphere(::Type{NF}; kwargs...) where {NF} = EarthDryAtmosphere{NF}(; kwargs...)
Base.eltype(::EarthDryAtmosphere{NF}) where {NF} = NF

lapse_rate(A::AbstractWetAtmosphere) = A.moist_lapse_rate
lapse_rate(A::AbstractDryAtmosphere) = A.dry_lapse_rate

latent_heat_condensation(A::AbstractWetAtmosphere) = A.latent_heat_condensation
latent_heat_condensation(A::AbstractDryAtmosphere) = zero(eltype(A))

latent_heat_sublimation(A::AbstractWetAtmosphere) = A.latent_heat_sublimation
latent_heat_sublimation(A::AbstractDryAtmosphere) = zero(eltype(A))