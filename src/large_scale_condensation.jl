function get_large_scale_condensation_tendencies!(
    Diag::DiagnosticVariables{NF},
    M::ModelSetup{NF},
) where {NF<:AbstractFloat}
    @unpack pres_ref, cp, alhc, alhs, gravity, τ, RH¹, ΔRH, rhb = M.constants  # Let's give some of these constants better names?
    @unpack nlon, nlat, nlev, σ_levels_half, σ_levels_full, σ_levels_thick =
        M.geospectral.geometry
    @unpack temp_grid, humid_grid, pres_surf_grid = Diag.grid_variables
    @unpack humid_saturation, cloud_top, precipitation_ls = Diag.parametrization_variables
    @unpack temp_tend, humid_tend = Diag.tendencies

    # Use mathematical notation for consistency with the documentation
    Q = humid_grid
    T = temp_grid
    ∂Q = humid_tend
    ∂T = temp_tend
    Qsat = humid_saturation

    p = exp.(pres_surf_grid)  # Normalised surface pressure - TODO(alistair): check pressure units throughout

    get_saturation_specific_humidity!(T, p, M)

    # 1. Tendencies of humidity and temperature due to large-scale condensation
    for k = 2:nlev
        σₖ = σ_levels_full[k]
        RH_threshold = RH¹ + ΔRH * (σₖ^2 - 1)  # Relative humidity threshold for condensation (Formula 24)
        if k == nlev
            RH_threshold = max(RH_threshold, rhb)
        end

        # Impose a maximum heating rate to avoid grid-point storm instability
        ∂Q_max = 10σₖ^2 / 3600τ  # This formula does not appear in the documentation

        for j = 1:nlat, i = 1:nlon
            Q_threshold = RH_threshold * Qsat[i, j, k]  # Specific humidity threshold for condensation
            if Q[i, j, k] > Q_threshold
                ∂Q[i, j, k] += -(Q[i, j, k] - Q_threshold) / τ  # Formula 22
                ∂T[i, j, k] += -alhc / cp * min(∂Q[i, j, k], ∂Q_max * p[i, j])  # Formula 23
                cloud_top[i, j] = min(k, cloud_top[i, j])  # Page 7 (last sentence)
            end
        end
    end

    # 2. Precipitation due to large-scale condensation#
    for k = 2:nlev
        Δpₖ = p * σ_levels_thick[k]  # Formula 4
        @. precipitation = -1 / grav * Δpₖ * ∂Q[:, :, k]  # Formula 25
    end
end
