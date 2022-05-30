function get_large_scale_condensation_tendencies!(
    Diag::DiagnosticVariables{NF},
    M   ::ModelSetup,
) where {NF<:AbstractFloat}
    @unpack gravity, humid_relax_time, RH_thresh_max, ΔRH, RH_thresh_boundary = M.constants  # Let's give some of these constants better names?
    @unpack cp, alhc = M.parameters
    @unpack nlon, nlat, nlev, σ_levels_full, σ_levels_thick = M.geospectral.geometry
    @unpack temp_grid, humid_grid, pres_surf_grid = Diag.grid_variables
    @unpack temp_tend, humid_tend, sat_spec_humidity, sat_vap_pressure, cloud_top, precip_large_scale = Diag.parametrization_variables

    pres = exp.(pres_surf_grid)  # Normalised surface pressure - TODO(alistair): check pressure units throughout

    get_saturation_specific_humidity!(sat_spec_humidity, sat_vap_pressure, temp_grid, pres, M)

    # 1. Tendencies of humidity and temperature due to large-scale condensation
    for k = 2:nlev
        σₖ = σ_levels_full[k]
        RH_threshold = RH_thresh_max + ΔRH * (σₖ^2 - 1)  # Relative humidity threshold for condensation (Formula 24)
        if k == nlev
            RH_threshold = max(RH_threshold, RH_thresh_boundary)
        end

        # Impose a maximum heating rate to avoid grid-point storm instability
        humid_tend_max = 10σₖ^2 / 3600humid_relax_time  # This formula does not appear in the documentation

        for j = 1:nlat, i = 1:nlon
            humid_threshold = RH_threshold * sat_spec_humidity[i, j, k]  # Specific humidity threshold for condensation
            if humid_grid[i, j, k] > humid_threshold
                humid_tend[i, j, k] += -(humid_grid[i, j, k] - humid_threshold) / humid_relax_time  # Formula 22
                temp_tend[i, j, k] += -alhc / cp * min(humid_tend[i, j, k], humid_tend_max * pres[i, j])  # Formula 23
                cloud_top[i, j] = min(k, cloud_top[i, j])  # Page 7 (last sentence)
            end
        end
    end

    # 2. Precipitation due to large-scale condensation#
    for k = 2:nlev
        Δpₖ = pres * σ_levels_thick[k]  # Formula 4
        for j = 1:nlat, i = 1:nlon
            precip_large_scale[i, j] += -1 / gravity * Δpₖ[i, j] * humid_tend[i, j, k]  # Formula 25
        end
    end
end
