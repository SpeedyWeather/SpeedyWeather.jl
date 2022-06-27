function get_large_scale_condensation_tendencies!(
    Diag::DiagnosticVariables{NF},
    M   ::ModelSetup,
) where {NF<:AbstractFloat}
    @unpack gravity, RH_thresh_max, RH_thresh_range, RH_thresh_boundary, humid_relax_time = M.constants
    @unpack cp, alhc, k_lsc = M.parameters
    @unpack nlon, nlat, nlev, σ_levels_full, σ_levels_thick = M.geospectral.geometry
    @unpack temp_grid, humid_grid, pres_grid = Diag.grid_variables
    @unpack temp_tend_lsc, humid_tend_lsc, sat_spec_humidity, sat_vap_pressure, cloud_top, precip_large_scale = Diag.parametrization_variables

    pres = exp.(pres_grid)  # Normalised surface pressure - TODO(alistair): check pressure units throughout

    get_saturation_specific_humidity!(sat_spec_humidity, sat_vap_pressure, temp_grid, pres, M)

    # 1. Tendencies of humidity and temperature due to large-scale condensation
    for k = k_lsc:nlev  # Used to be 2:nlev in original Speedy
        σₖ = σ_levels_full[k]
        RH_threshold = RH_thresh_max + RH_thresh_range * (σₖ^2 - 1)  # Relative humidity threshold for condensation (Formula 24)
        if k == nlev
            RH_threshold = max(RH_threshold, RH_thresh_boundary)
        end

        # Impose a maximum heating rate to avoid grid-point storm instability
        # This formula does not appear in the original Speedy documentation
        humid_tend_max = 10σₖ^2 / 3600humid_relax_time

        for j = 1:nlat, i = 1:nlon
            humid_threshold = RH_threshold * sat_spec_humidity[i, j, k]  # Specific humidity threshold for condensation
            if humid_grid[i, j, k] > humid_threshold
                humid_tend_lsc[i, j, k] = -(humid_grid[i, j, k] - humid_threshold) / humid_relax_time            # Formula 22
                temp_tend_lsc[i, j, k] = -alhc / cp * min(humid_tend_lsc[i, j, k], humid_tend_max * pres[i, j])  # Formula 23
                cloud_top[i, j] = min(cloud_top[i, j], k)                                                        # Page 7 (last sentence)
            else
                humid_tend_lsc[i, j, k] = zero(NF)
                temp_tend_lsc[i, j, k] = zero(NF)
            end
        end
    end

    # 2. Precipitation due to large-scale condensation
    fill!(precip_large_scale, zero(NF))
    for k = k_lsc:nlev
        Δpₖ = pres * σ_levels_thick[k]  # Formula 4
        for j = 1:nlat, i = 1:nlon
            precip_large_scale[i, j] += -1 / gravity * Δpₖ[i, j] * humid_tend_lsc[i, j, k]  # Formula 25
        end
    end
end
