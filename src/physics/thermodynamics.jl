abstract type AbstractClausiusClapeyron <: AbstractParameterization end

export ClausiusClapeyron

"""
Parameters for computing saturation vapour pressure of water using the Clausis-Clapeyron equation,

    e(T) = e₀ * exp( -Lᵥ/Rᵥ * (1/T - 1/T₀)),

where T is in Kelvin, Lᵥ the the latent heat of condensation and Rᵥ the gas constant of water vapour
$(TYPEDFIELDS)"""
Base.@kwdef struct ClausiusClapeyron{NF<:AbstractFloat} <: AbstractClausiusClapeyron
    "Saturation water vapour pressure at 0°C [Pa]"
    e₀::NF = 610.78

    "0°C in Kelvin"
    T₀::NF = 273.16

    "Latent heat of condensation/vaporization of water [J/kg]"
    Lᵥ::NF = 2.5e6

    "Specific heat at constant pressure [J/K/kg]"
    cₚ::NF = 1004.64

    "Gas constant of water vapour [J/kg/K]"
    R_vapour::NF = 461.5

    "Gas constant of dry air [J/kg/K]"
    R_dry::NF = 287.04

    "Latent heat of condensation divided by gas constant of water vapour [K]"    
    Lᵥ_Rᵥ::NF = Lᵥ/R_vapour

    "Inverse of T₀, one over 0°C in Kelvin"
    T₀⁻¹::NF = inv(T₀)

    "Ratio of molecular masses [1] of water vapour over dry air (=R_dry/R_vapour)."
    mol_ratio::NF = R_dry/R_vapour
end

# generator function
function ClausiusClapeyron(SG::SpectralGrid, atm::AbstractAtmosphere; kwargs...)
    (;R_dry, R_vapour, latent_heat_condensation, heat_capacity) = atm
    return ClausiusClapeyron{SG.NF}(;Lᵥ=latent_heat_condensation,R_dry,R_vapour,cₚ=heat_capacity,kwargs...)
end

"""
$(TYPEDSIGNATURES)
Functor: Saturation water vapour pressure as a function of temperature using the
Clausius-Clapeyron equation,

    e(T) = e₀ * exp( -Lᵥ/Rᵥ * (1/T - 1/T₀)),

where T is in Kelvin, Lᵥ the the latent heat of vaporization and Rᵥ the gas constant
of water vapour, T₀ is 0˚C in Kelvin."""
function (CC::ClausiusClapeyron{NF})(temp_kelvin::NF) where NF
    (;e₀, T₀⁻¹, Lᵥ_Rᵥ) = CC
    return e₀ * exp(Lᵥ_Rᵥ*(T₀⁻¹ - inv(temp_kelvin)))
end

# convert to number format of struct
function (CC::ClausiusClapeyron{NF})(temp_kelvin) where NF
    CC(convert(NF,temp_kelvin))
end

"""
$(TYPEDSIGNATURES)
Gradient of Clausius-Clapeyron wrt to temperature, evaluated at `temp_kelvin`."""
function grad(CC::ClausiusClapeyron{NF},temp_kelvin::NF) where NF
    e = CC(temp_kelvin)
    return e*CC.Lᵥ_Rᵥ/temp_kelvin^2
end

# convert to input argument to number format from struct
grad(CC::ClausiusClapeyron{NF},temp_kelvin) where NF = grad(CC,convert(NF,temp_kelvin))

"""
$(TYPEDSIGNATURES)
Saturation humidity from saturation vapour pressure and pressure via

    qsat = mol_ratio*sat_vap_pres/pres

with both pressures in same units and qsat in kg/kg."""
function saturation_humidity(
    sat_vap_pres::NF,                   # saturation vapour pressure [Pa]
    pres::NF;                           # pressure [Pa]
    mol_ratio::NF = NF(287.04/461.5)    # ratio of mol masses dry air - water vapour
) where NF
    # simpler version given that pres >> sat_vap_pres for most temperatures T<40˚C
    return mol_ratio*sat_vap_pres / pres
    # more accurate version, with a more complicated derivative though
    # return mol_ratio*sat_vap_pres / (pres - (1-mol_ratio)*sat_vap_pres)
end

"""
$(TYPEDSIGNATURES)
Saturation humidity [kg/kg] from temperature [K], pressure [Pa] via

    sat_vap_pres = clausius_clapeyron(temperature)
    saturation humidity = mol_ratio * sat_vap_pres / pressure"""
function saturation_humidity(
    temp_kelvin::NF,
    pres::NF,
    clausius_clapeyron::AbstractClausiusClapeyron,
) where NF
    sat_vap_pres = clausius_clapeyron(temp_kelvin)
    return saturation_humidity(sat_vap_pres,pres;mol_ratio=clausius_clapeyron.mol_ratio)
end

"""
$(TYPEDSIGNATURES)
Gradient of Clausius-Clapeyron wrt to temperature, evaluated at `temp_kelvin`."""
function grad_saturation_humidity(CC::ClausiusClapeyron{NF},temp_kelvin::NF,pres::NF) where NF
    qsat = saturation_humidity(temp_kelvin,pres,CC)
    return qsat*CC.Lᵥ_Rᵥ/temp_kelvin^2
end

"""
$(TYPEDSIGNATURES)
Calculate geopotentiala and dry static energy for the primitive equation model."""
function get_thermodynamics!(column::ColumnVariables,model::PrimitiveEquation)
    geopotential!(column.geopot, column.temp, model.geopotential, column.surface_geopotential)
    dry_static_energy!(column, model.atmosphere)
end

"""
$(TYPEDSIGNATURES)
Compute the dry static energy SE = cₚT + Φ (latent heat times temperature plus geopotential)
for the column."""
function dry_static_energy!(
    column::ColumnVariables,
    atmosphere::AbstractAtmosphere
)
    cₚ = atmosphere.heat_capacity
    (;dry_static_energy, geopot, temp) = column

    @inbounds for k in eachlayer(column)
        dry_static_energy[k] = cₚ * temp[k] + geopot[k]
    end

    return nothing
end

function bulk_richardson!(
    column::ColumnVariables,
    atmosphere::AbstractAtmosphere,
)
    cₚ = atmosphere.heat_capacity
    (;u, v, geopot, temp_virt, nlev, bulk_richardson) = column

    V² = u[nlev]^2 + v[nlev]^2
    Θ₀ = cₚ*temp_virt[nlev]
    Θ₁ = Θ₀ + geopot[nlev]
    bulk_richardson[nlev] = geopot[nlev]*(Θ₁ - Θ₀)/Θ₀/V²

    @inbounds for k in nlev-1:-1:1
        V² = u[k]^2 + v[k]^2
        virtual_dry_static_energy = cₚ * temp_virt[k] + geopot[k]
        bulk_richardson[k] = geopot[k]*(virtual_dry_static_energy - Θ₁)/Θ₁/V²
    end

    return nothing
end

"""$(TYPEDSIGNATURES)
Compute the saturation water vapour pressure [Pa], the saturation humidity [kg/kg]
and the relative humidity following `clausius_clapeyron`."""
function saturation_humidity!(
    column::ColumnVariables,
    clausius_clapeyron::AbstractClausiusClapeyron,
)
    (;sat_humid, pres, temp) = column

    for k in eachlayer(column)
        sat_humid[k] = saturation_humidity(temp[k],pres[k],clausius_clapeyron)
    end
end

"""$(TYPEDSIGNATURES)
Compute the moist static energy

    MSE = SE + Lc*Q = cₚT + Φ + Lc*Q

with the static energy `SE`, the latent heat of condensation `Lc`,
the geopotential `Φ`. As well as the saturation moist static energy
which replaces Q with Q_sat"""
function moist_static_energy!(
    column::ColumnVariables,
    clausius_clapeyron::AbstractClausiusClapeyron
)
    (;Lᵥ) = clausius_clapeyron      # latent heat of vaporization
    (;sat_moist_static_energy, moist_static_energy, dry_static_energy) = column
    (;humid, sat_humid) = column

    for k in eachlayer(column)
        moist_static_energy[k] = dry_static_energy[k] + Lᵥ * humid[k]
        sat_moist_static_energy[k] = dry_static_energy[k] + Lᵥ * sat_humid[k]
    end
end

export TetensEquation

"""
Parameters for computing saturation vapour pressure of water using the Tetens' equation,

    eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T + Tᵢ)),

where T is in Kelvin and i = 1,2 for saturation above and below freezing,
respectively. From Tetens (1930), and Murray (1967) for below freezing.
$(TYPEDFIELDS)"""
Base.@kwdef struct TetensEquation{NF<:AbstractFloat} <: AbstractClausiusClapeyron
    "Saturation water vapour pressure at 0°C [Pa]"
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
Functor: Saturation water vapour pressure as a function of temperature using the
Tetens equation,

    eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T - Tᵢ)),

where T is in Kelvin and i = 1,2 for saturation above and below freezing,
respectively."""
function (TetensCoefficients::TetensEquation{NF})(temp_kelvin::NF) where NF
    (;e₀, T₀, C₁, C₂, T₁, T₂) = TetensCoefficients
    C, T = temp_kelvin > T₀ ? (C₁, T₁) : (C₂, T₂)      # change coefficients above/below freezing
    temp_celsius = temp_kelvin - T₀
    return e₀ * exp(C * temp_celsius / (temp_celsius + T))
end

"""
$(TYPEDSIGNATURES)
Gradient of the Tetens equation wrt to temperature, evaluated at `temp_kelvin`."""
function grad(TetensCoefficients::TetensEquation{NF},temp_kelvin::NF) where NF
    (;T₀, C₁, C₂, T₁, T₂) = TetensCoefficients
    e = TetensCoefficients(temp_kelvin)             # saturation vapour pressure
    C, T = temp_kelvin > T₀ ? (C₁, T₁) : (C₂, T₂)   # change coefficients above/below freezing
    return e*C*T/(temp_kelvin - T₀ - T)^2           # chain rule: times derivative of inner function
end