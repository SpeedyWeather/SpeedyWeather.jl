"""
Parameters for computing saturation vapour pressure of water using the Tetens' equation,

    eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T + Tᵢ)),

where T is in Kelvin and i = 1,2 for saturation above and below freezing,
respectively. From Tetens (1930), and Murray (1967) for below freezing.
$(TYPEDFIELDS)"""
Base.@kwdef struct TetensCoefs{NF<:AbstractFloat}
    "Saturation water vapour pressure at 0°C [Pa]"
    e₀::NF = 610.78

    "0°C in Kelvin"
    T₀::NF = 273.15

    "Tetens denominator (water) [˚C]"
    T₁::NF = 237.3

    "Tetens denominator following Murray (1967, ice) [˚C]"
    T₂::NF = 265.5

    "Tetens numerator scaling [1], above freezing"
    C₁::NF = 17.27

    "Tetens numerator scaling [1], below freezing"
    C₂::NF = 21.875
end

Base.@kwdef struct Thermodynamics{NF} <: AbstractThermodynamics{NF}
    tetens_coefs::TetensCoefs{NF} = TetensCoefs{NF}()
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
    geopotential!(column.geopot,column.temp,model.constants)
    dry_static_energy!(column, model.constants)
end

"""
$(TYPEDSIGNATURES)
Calculate thermodynamic quantities like saturation vapour pressure,
saturation specific humidity, dry static energy, moist static energy
and saturation moist static energy from the prognostic column variables."""
function get_thermodynamics!(column::ColumnVariables,model::PrimitiveWet)
    geopotential!(column.geopot,column.temp,model.constants)
    dry_static_energy!(column, model.constants)
    saturation_humidity!(column, model.thermodynamics)
    moist_static_energy!(column, model.thermodynamics)

    # Interpolate certain variables to half-levels
    vertical_interpolate!(column, model)

    return nothing
end

"""
$(TYPEDSIGNATURES)
Full to half-level interpolation for humidity, saturation humidity,
dry static energy and saturation moist static energy.
"""
function vertical_interpolate!(
    column::ColumnVariables,
    model::PrimitiveEquation,
)

    for (full, half) in (
        (column.humid,                      column.humid_half),
        (column.sat_humid,                  column.sat_humid_half),
        (column.dry_static_energy,          column.dry_static_energy_half),
        (column.sat_moist_static_energy,    column.sat_moist_static_energy_half),
    )
        vertical_interpolate!(half, full, model.geometry)
    end
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
Compute (1) the saturation water vapour pressure as a function of temperature using the
Tetens equation,

    eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T - Tᵢ)),

where T is in Kelvin and i = 1,2 for saturation above and below freezing,
respectively. And (2) the saturation specific humidity according to the formula,

    0.622 * e / (p - (1 - 0.622) * e),

where `e` is the saturation vapour pressure, `p` is the pressure, and 0.622 is the ratio of
the molecular weight of water to dry air."""
function saturation_humidity!(
    column::ColumnVariables,
    thermodynamics::Thermodynamics,
)
    (;sat_humid, rel_humid, sat_vap_pres, pres, temp, humid) = column
    (;mol_ratio) = thermodynamics      # = mol_mass_vapour/mol_mass_dry_air = 0.622

    for k in eachlayer(column)
        sat_vap_pres[k] = saturation_vapour_pressure(temp[k],thermodynamics.tetens_coefs)
        sat_humid[k] = saturation_humidity(sat_vap_pres[k],pres[k];mol_ratio)
        rel_humid[k] = humid[k]/sat_humid[k]
        # rel_humid[k] > 2 && @info "$(column.ij), $k"
    end
end

# integer version to convert to Float64
saturation_vapour_pressure(temp_kelvin::Integer) = saturation_vapour_pressure(Float64(temp_kelvin))

# version that pulls default Tetens coefficients
function saturation_vapour_pressure(temp_kelvin::NF) where NF
    return saturation_vapour_pressure(temp_kelvin,TetensCoefs{NF}())
end

function relative_humidity(temp_kelvin,humid,pres)
    sat_vap_pres = saturation_vapour_pressure(temp_kelvin)
    sat_vap_pres, pres = promote(sat_vap_pres, pres)
    sat_humid = saturation_humidity(sat_vap_pres,pres)
    return humid/sat_humid  # = relative humidity [1]
end

"""
$(TYPEDSIGNATURES)
Saturation water vapour pressure as a function of temperature using the
Tetens equation,

    eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T - Tᵢ)),

where T is in Kelvin and i = 1,2 for saturation above and below freezing,
respectively."""
function saturation_vapour_pressure(temp_kelvin::NF,tetens_coefs::TetensCoefs{NF}) where NF
    (;e₀, T₀, C₁, C₂, T₁, T₂) = tetens_coefs
    C, T = temp_kelvin > T₀ ? (C₁, T₁) : (C₂, T₂)      # change coefficients above/below freezing
    temp_celsius = temp_kelvin - T₀
    return e₀ * exp(C * temp_celsius / (temp_celsius + T))
end

function saturation_humidity(sat_vap_pres::NF,pres::NF;mol_ratio::NF = NF(0.622)) where NF
    return mol_ratio*sat_vap_pres / (pres - (1-mol_ratio)*sat_vap_pres)
end

function saturation_humidity(temp::NF,pres::NF,thermodynamics::Thermodynamics{NF}) where NF
    (;mol_ratio) = thermodynamics
    sat_vap_pres = saturation_vapour_pressure(temp,thermodynamics.tetens_coefs)
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
