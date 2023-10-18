"""
Parameters for computing saturation vapour pressure using the August-Roche-Magnus formula,

    eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T - Tᵢ)),

where T is in Kelvin and i = 1,2 for saturation with respect to water and ice,
respectively.
$(TYPEDFIELDS)"""
Base.@kwdef struct MagnusCoefs{NF<:AbstractFloat}
    "Saturation vapour pressure at 0°C [Pa]"
    e₀::NF = 610.8

    "0°C in Kelvin"
    T₀::NF = 273.16
    T₁::NF = 35.86
    T₂::NF = 7.66
    C₁::NF = 17.269
    C₂::NF = 21.875
end

Base.@kwdef struct Thermodynamics{NF} <: AbstractThermodynamics{NF}
    magnus_coefs::MagnusCoefs{NF} = MagnusCoefs{NF}()
    mol_ratio::NF 
    latent_heat_condensation::NF
    latent_heat_sublimation::NF
end

function Thermodynamics(SG::SpectralGrid,atm::AbstractAtmosphere;kwargs...)
    (;latent_heat_condensation, latent_heat_sublimation) = atm
    mol_ratio = atm.mol_mass_vapour/atm.mol_mass_dry_air
    return Thermodynamics{SG.NF}(;mol_ratio,latent_heat_condensation,latent_heat_sublimation,kwargs...)
end

"""
$(TYPEDSIGNATURES)
Calculate the dry static energy for the primitive dry model."""
function get_thermodynamics!(column::ColumnVariables,model::PrimitiveDry)
    dry_static_energy!(column, model.constants)
end

"""
$(TYPEDSIGNATURES)
Calculate thermodynamic quantities like saturation vapour pressure,
saturation specific humidity, dry static energy, moist static energy
and saturation moist static energy from the prognostic column variables."""
function get_thermodynamics!(column::ColumnVariables,model::PrimitiveWet)
    dry_static_energy!(column, model.constants)
    saturation_humidity!(column, model.thermodynamics)
    moist_static_energy!(column, model.thermodynamics)

    # Interpolate certain variables to half-levels
    # interpolate!(column, model)

    return nothing
end

"""
$(TYPEDSIGNATURES)
Compute the dry static energy SE = cₚT + Φ (latent heat times temperature plus geopotential)
for the column."""
function dry_static_energy!(column::ColumnVariables,constants::DynamicsConstants)

    (;cₚ) = constants
    (;dry_static_energy, geopot, temp) = column

    @inbounds for k in eachlayer(column)
        dry_static_energy[k] = cₚ * temp[k] + geopot[k]
    end

    return nothing
end

"""$(TYPEDSIGNATURES)
Compute (1) the saturation vapour pressure as a function of temperature using the
August-Roche-Magnus formula,

    eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T - Tᵢ)),

where T is in Kelvin and i = 1,2 for saturation with respect to water and ice,
respectively. And (2) the saturation specific humidity according to the formula,

    0.622 * e / (p - (1 - 0.622) * e),

where `e` is the saturation vapour pressure, `p` is the pressure, and 0.622 is the ratio of
the molecular weight of water to dry air."""
function saturation_humidity!(
    column::ColumnVariables,
    thermodynamics::Thermodynamics,
)
    (;sat_humid, sat_vap_pres, pres, temp) = column
    # (;e₀, T₀, C₁, C₂, T₁, T₂) = thermodynamics.magnus_coefs
    (;mol_ratio) = thermodynamics      # = mol_mass_vapour/mol_mass_dry_air = 0.622

    for k in eachlayer(column)
        # change coefficients for water (temp > T₀) or ice (else)
        # C, T = temp[k] > T₀ ? (C₁, T₁) : (C₂, T₂)
        # sat_vap_pres[k] = e₀ * exp(C * (temp[k] - T₀) / (temp[k] - T))
        sat_vap_pres[k] = saturation_vapour_pressure(temp[k],thermodynamics.magnus_coefs)
        # sat_humid[k] = mol_ratio*sat_vap_pres[k] / (pres[k] - (1-mol_ratio)*sat_vap_pres[k])
        sat_humid[k] = saturation_humidity(sat_vap_pres[k],pres[k];mol_ratio)
    end
end

function saturation_vapour_pressure(temp::NF,magnus_coefs::MagnusCoefs{NF}) where NF
    (;e₀, T₀, C₁, C₂, T₁, T₂) = magnus_coefs
    C, T = temp > T₀ ? (C₁, T₁) : (C₂, T₂)
    return e₀ * exp(C * (temp - T₀) / (temp - T))
end

function saturation_humidity(sat_vap_pres::NF,pres::NF;mol_ratio::NF = NF(0.622)) where NF
    return mol_ratio*sat_vap_pres / (pres - (1-mol_ratio)*sat_vap_pres)
end

function saturation_humidity(temp::NF,pres::NF,thermodynamics::Thermodynamics{NF}) where NF
    (;mol_ratio) = thermodynamics
    sat_vap_pres = saturation_vapour_pressure(temp,thermodynamics.magnus_coefs)
    return saturation_humidity(sat_vap_pres,pres;mol_ratio)
end

"""$(TYPEDSIGNATURES)
Compute the moist static energy

    MSE = SE + Lc*Q = cₚT + Φ + Lc*Q

with the static energy `SE`, the latent heat of condensation `Lc`,
the geopotential `Φ`. As well as the saturation moist static energy
which replaces Q with Q_sat"""
function moist_static_energy!(column::ColumnVariables,thermodynamics::Thermodynamics)
    (;latent_heat_condensation) = thermodynamics
    (;sat_moist_static_energy, moist_static_energy, dry_static_energy) = column
    (;humid, sat_humid) = column

    for k in eachlayer(column)
        moist_static_energy[k] = dry_static_energy[k] + latent_heat_condensation * humid[k]
        sat_moist_static_energy[k] = dry_static_energy[k] + latent_heat_condensation * sat_humid[k]
    end
end
