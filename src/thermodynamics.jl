"""
    get_thermodynamics!(
        column::ColumnVariables{NF},
        model::PrimitiveEquationModel{NF},
    )
"""
function get_thermodynamics!(column, model)
    saturation_vapour_pressure!(column, model)
    saturation_specific_humidity!(column, model)
    dry_static_energy!(column, model)
    moist_static_energy!(column, model)
    saturation_moist_static_energy!(column, model)
    return nothing
end

"""
saturation_vapour_pressure!(column::ColumnVariables,
                            model::PrimitiveEquationModel)

Compute the saturation vapour pressure as a function of temperature using the
August-Roche-Magnus formula,

eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T - Tᵢ)),

where T is in Kelvin and i = 1,2 for saturation with respect to water and ice,
respectively.
"""
function saturation_vapour_pressure!(
    column::ColumnVariables{NF},
    model::PrimitiveEquationModel{NF},
) where {NF<:AbstractFloat}

    @unpack sat_vap_pres, temp = column
    @unpack e₀, T₀, C₁, C₂, T₁, T₂ = model.parameters.magnus_coefs

    for k in eachlayer(column)
        # change coefficients for water (temp > T₀) or ice (else)
        C, T = temp[k] > T₀ ? (C₁, T₁) : (C₂, T₂)
        sat_vap_pres[k] = e₀ * exp(C * (temp[k] - T₀) / (temp[k] - T))
    end

    return nothing
end

"""
saturation_specific_humidity!(  column::ColumnVariables{NF},
                                model::PrimitiveEquationModel{NF}
                                ) where {NF<:AbstractFloat}

Compute the saturation specific humidity according to the formula,

0.622 * e / (p - (1 - 0.622) * e),

where e is the saturation vapour pressure, p is the pressure, and 0.622 is the ratio of
the molecular weight of water to dry air.
"""
function saturation_specific_humidity!(
    column::ColumnVariables{NF},
    model::PrimitiveEquationModel{NF},
) where {NF<:AbstractFloat}

    @unpack σ_levels_full = model.geometry
    @unpack sat_humid, sat_vap_pres, pres = column

    mol_ratio = convert(NF, 0.622)
    one_minus_mol_ratio = convert(NF, 1 - mol_ratio)

    for k in eachlayer(column)
        pₖ = pres * σ_levels_full[k]       # pressure in layer k
        sat_humid[k] =
            mol_ratio * sat_vap_pres[k] / (pₖ - one_minus_mol_ratio * sat_vap_pres[k])
    end

    return nothing
end

"""
    dry_moist_static_energy!(
        column::ColumnVariables{NF},
        model::PrimitiveEquationModel{NF},
    )
"""
function dry_static_energy!(
    column::ColumnVariables{NF},
    model::PrimitiveEquationModel{NF},
) where {NF<:AbstractFloat}
    @unpack cp = model.parameters
    @unpack dry_static_energy, geopot, temp = column

    for k in eachlayer(column)
        dry_static_energy[k] = cp * temp[k] + geopot[k]
    end

    return nothing
end

"""
    moist_static_energy!(
        column::ColumnVariables{NF},
        model::PrimitiveEquationModel{NF},
    )
"""
function moist_static_energy!(
    column::ColumnVariables{NF},
    model::PrimitiveEquationModel{NF},
) where {NF<:AbstractFloat}
    @unpack alhc = model.parameters
    @unpack moist_static_energy, dry_static_energy, humid = column

    for k in eachlayer(column)
        moist_static_energy[k] = dry_static_energy[k] + alhc * humid[k]
    end

    return nothing
end

"""
    saturation_moist_static_energy!(
        column::ColumnVariables{NF},
        model::PrimitiveEquationModel{NF},
    )
"""
function saturation_moist_static_energy!(
    column::ColumnVariables{NF},
    model::PrimitiveEquationModel{NF},
) where {NF<:AbstractFloat}
    @unpack alhc = model.parameters
    @unpack sat_moist_static_energy, sat_moist_static_energy_half, dry_static_energy, sat_humid = column

    # Full levels
    for k in eachlayer(column)
        sat_moist_static_energy[k] = dry_static_energy[k] + alhc * sat_humid[k]
    end

    # Half-levels
    interpolate!(sat_moist_static_energy, sat_moist_static_energy_half, column, model)

    return nothing
end

"""
    interpolate!(A_full_level, A_half_level, column, model)

Given some generic column variable A defined at full levels, do a linear interpolation in
log(σ) to calculate its values at half-levels.
"""
function interpolate!(
    A_full_level::Vector{NF},
    A_half_level::Vector{NF},
    column::ColumnVariables{NF},
    model::PrimitiveEquationModel,
) where {NF<:AbstractFloat}
    @unpack nlev = column
    @unpack σ_levels_full, σ_levels_half = model.geometry

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
        (log(0.99) - log(σ_levels_full[nlev])) /
        (log(σ_levels_full[nlev]) - log(σ_levels_full[nlev-1]))

    return nothing
end
