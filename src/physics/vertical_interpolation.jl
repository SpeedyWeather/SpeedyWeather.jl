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