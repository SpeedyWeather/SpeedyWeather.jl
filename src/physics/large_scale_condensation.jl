function large_scale_condensation!( column::ColumnVariables{NF},
                                    model::PrimitiveEquation,
                                    ) where {NF<:AbstractFloat}

    (; gravity, RH_thresh_max_lsc, RH_thresh_range_lsc, RH_thresh_pbl_lsc, humid_relax_time_lsc ) = model.constants
    (; σ_levels_full, σ_levels_thick, n_stratosphere_levels ) = model.geometry
    (; cₚ, alhc ) = model.parameters

    (; humid, pres ) = column            # prognostic variables: specific humidity, surface pressure
    (; temp_tend, humid_tend ) = column  # tendencies to write into
    (; sat_humid ) = column              # intermediate variable
    (; nlev ) = column

    # 1. Tendencies of humidity and temperature due to large-scale condensation
    for k in eachlayer(column)[n_stratosphere_levels:end]           # top to bottom, skip stratospheric levels

        # Relative humidity threshold for condensation (Formula 24)
        σₖ² = σ_levels_full[k]^2
        RH_threshold = RH_thresh_max_lsc + RH_thresh_range_lsc * (σₖ² - 1)
        if k == nlev
            RH_threshold = max(RH_threshold, RH_thresh_pbl_lsc)
        end

        # Impose a maximum heating rate to avoid grid-point storm instability
        # This formula does not appear in the original Speedy documentation
        humid_tend_max = 10σₖ² / 3600humid_relax_time_lsc
        humid_threshold = RH_threshold * sat_humid[k]               # Specific humidity threshold for condensation

        if humid[k] > humid_threshold
            # accumulate in tendencies (nothing is added if humidity not above threshold)
            humid_tend[k] += -(humid[k] - humid_threshold) / humid_relax_time_lsc   # Formula 22
            temp_tend[k] += -alhc / cₚ * min(humid_tend[k], humid_tend_max*pres)    # Formula 23

            # If there is large-scale condensation at a level higher (i.e. smaller k) than
            # the cloud-top previously diagnosed due to convection, then increase the
            # cloud-top.
            column.cloud_top = min(column.cloud_top, k)                             # Page 7 (last sentence)
        end
    end

    # 2. Precipitation due to large-scale condensation
    pₛ = pres[end]                                                          # surface pressure
    for k in eachlayer(column)[n_stratosphere_levels:end]                   # top to bottom, skip stratosphere
        Δpₖ = pₛ*σ_levels_thick[k]                                          # Formula 4  # Correct index?
        column.precip_large_scale += -1 / gravity * Δpₖ * humid_tend[k]     # Formula 25
    end
end
