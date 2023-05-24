"""
Parameters for computing saturation vapour pressure using the August-Roche-Magnus formula,

    eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T - Tᵢ)),

where T is in Kelvin and i = 1,2 for saturation with respect to water and ice,
respectively.
"""
Base.@kwdef struct MagnusCoefs{NF<:Real} <: Coefficients
    e₀::NF = 6.108   # Saturation vapour pressure at 0°C
    T₀::NF = 273.16  # 0°C in Kelvin
    T₁::NF = 35.86
    T₂::NF = 7.66
    C₁::NF = 17.269
    C₂::NF = 21.875
end

"""
    get_thermodynamics!(column::ColumnVariables,model::PrimitiveWetCore)

Calculate thermodynamic quantities like saturation vapour pressure,
saturation specific humidity, dry static energy, moist static energy
and saturation moist static energy from the prognostic column variables."""
function get_thermodynamics!(   column::ColumnVariables,
                                model::PrimitiveEquation)

    # Calculate thermodynamic quantities at full levels
    dry_static_energy!(column, model)

    if model isa PrimitiveWetCore
        saturation_vapour_pressure!(column, model)
        saturation_specific_humidity!(column, model)
        moist_static_energy!(column, model)
        saturation_moist_static_energy!(column, model)
    end

    # Interpolate certain variables to half-levels
    # interpolate!(column, model)

    return nothing
end

"""
    interpolate!(
        column::ColumnVariables{NF},
        model::PrimitiveEquation,
    )
"""
function interpolate!(
    column::ColumnVariables{NF},
    model::PrimitiveEquation,
) where {NF<:AbstractFloat}
    (; humid, humid_half ) = column
    (; sat_humid, sat_humid_half ) = column
    (; dry_static_energy, dry_static_energy_half ) = column
    (; sat_moist_static_energy, sat_moist_static_energy_half ) = column

    for (full, half) in (
        (humid, humid_half),
        (sat_humid, sat_humid_half),
        (dry_static_energy, dry_static_energy_half),
        (sat_moist_static_energy, sat_moist_static_energy_half),
    )
        interpolate!(full, half, column, model)
    end

    return nothing
end

"""
    interpolate!(
        A_full_level::Vector{NF},
        A_half_level::Vector{NF},
        column::ColumnVariables{NF},
        model::PrimitiveEquation,
    )

Given some generic column variable A defined at full levels, do a linear interpolation in
log(σ) to calculate its values at half-levels.
"""
function interpolate!(
    A_full_level::Vector{NF},
    A_half_level::Vector{NF},
    column::ColumnVariables{NF},
    model::PrimitiveEquation,
) where {NF<:AbstractFloat}
    (; nlev ) = column
    (; σ_levels_full, σ_levels_half ) = model.geometry

    # For A at each full level k, compute A at the half-level below, i.e. at the boundary
    # between the full levels k and k+1.
    for k = 1:nlev-1
        A_half_level[k] =
            A_full_level[k] +
            (A_full_level[k+1] - A_full_level[k]) *
            (log(σ_levels_half[k+1]) - log(σ_levels_full[k])) /
            (log(σ_levels_full[k+1]) - log(σ_levels_full[k]))
    end

    # Compute the values at the surface separately
    A_half_level[nlev] =
        A_full_level[nlev] +
        (A_full_level[nlev] - A_full_level[nlev-1]) *
        (log(NF(0.99)) - log(σ_levels_full[nlev])) /
        (log(σ_levels_full[nlev]) - log(σ_levels_full[nlev-1]))

    return nothing
end


"""
    saturation_vapour_pressure!(
        column::ColumnVariables{NF},
        model::PrimitiveEquation,
    )

Compute the saturation vapour pressure as a function of temperature using the
August-Roche-Magnus formula,

eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T - Tᵢ)),

where T is in Kelvin and i = 1,2 for saturation with respect to water and ice,
respectively.
"""
function saturation_vapour_pressure!(
    column::ColumnVariables{NF},
    model::PrimitiveEquation,
) where {NF<:AbstractFloat}

    (; sat_vap_pres, temp ) = column
    (; e₀, T₀, C₁, C₂, T₁, T₂ ) = model.parameters.magnus_coefs

    for k in eachlayer(column)
        # change coefficients for water (temp > T₀) or ice (else)
        C, T = temp[k] > T₀ ? (C₁, T₁) : (C₂, T₂)
        sat_vap_pres[k] = e₀ * exp(C * (temp[k] - T₀) / (temp[k] - T))
    end

    return nothing
end

"""
    saturation_specific_humidity!(  column::ColumnVariables{NF},
                                    model::PrimitiveEquation                                    ) where {NF<:AbstractFloat}

Compute the saturation specific humidity according to the formula,

0.622 * e / (p - (1 - 0.622) * e),

where e is the saturation vapour pressure, p is the pressure, and 0.622 is the ratio of
the molecular weight of water to dry air.
"""
function saturation_specific_humidity!(
    column::ColumnVariables{NF},
    model::PrimitiveEquation,
) where {NF<:AbstractFloat}

    (; sat_humid, sat_vap_pres, pres ) = column
    (; mol_mass_vapour, mol_mass_dry_air ) = model.parameters

    mol_ratio = convert(NF, mol_mass_vapour/mol_mass_dry_air)

    @inbounds for k in eachlayer(column)
        sat_humid[k] = mol_ratio*sat_vap_pres[k] / (pres[k] - (1-mol_ratio)*sat_vap_pres[k])
    end

    return nothing
end

"""
    dry_static_energy!(column::ColumnVariables,model::PrimitiveEquation)

Compute the dry static energy SE = cₚT + Φ (latent heat times temperature plus geopotential)
for the column."""
function dry_static_energy!(column::ColumnVariables{NF},
                            model::PrimitiveEquation) where NF

    cₚ = convert(NF,model.parameters.cₚ)
    (;dry_static_energy, geopot, temp) = column

    @inbounds for k in eachlayer(column)
        dry_static_energy[k] = cₚ * temp[k] + geopot[k]
    end

    return nothing
end

"""
    moist_static_energy!(
        column::ColumnVariables{NF},
        model::PrimitiveEquation,
    )
"""
function moist_static_energy!(
    column::ColumnVariables{NF},
    model::PrimitiveEquation,
) where {NF<:AbstractFloat}
    (; alhc ) = model.parameters
    (; moist_static_energy, dry_static_energy, humid ) = column

    for k in eachlayer(column)
        moist_static_energy[k] = dry_static_energy[k] + alhc * humid[k]
    end

    return nothing
end

"""
    saturation_moist_static_energy!(
        column::ColumnVariables{NF},
        model::PrimitiveEquation,
    )
"""
function saturation_moist_static_energy!(
    column::ColumnVariables{NF},
    model::PrimitiveEquation,
) where {NF<:AbstractFloat}
    (; alhc ) = model.parameters
    (; sat_moist_static_energy,
    sat_moist_static_energy_half,
    dry_static_energy,
    sat_humid) = column

    for k in eachlayer(column)
        sat_moist_static_energy[k] = dry_static_energy[k] + alhc * sat_humid[k]
    end

    return nothing
end
