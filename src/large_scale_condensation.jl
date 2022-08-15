function large_scale_condensation!( column::ColumnVariables{NF},
                                    model::PrimitiveEquationModel,
                                    ) where {NF<:AbstractFloat}

    @unpack gravity, RH_thresh_max, RH_thresh_range, RH_thresh_boundary, humid_relax_time = model.constants
    @unpack cp, alhc, k_lsc = model.parameters
    @unpack nlon, nlat, nlev, σ_levels_full, σ_levels_thick = M.geometry

    get_saturation_vapour_pressure!(column, model)
    get_saturation_specific_humidity!(column, model)

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
                humid_tend_lsc[i, j, k] = -(humid_grid[i, j, k] - humid_threshold) / humid_relax_time                 # Formula 22
                temp_tend_lsc[i, j, k] = -alhc / cp * min(humid_tend_lsc[i, j, k], humid_tend_max * pres_grid[i, j])  # Formula 23
                cloud_top[i, j] = min(cloud_top[i, j], k)                                                             # Page 7 (last sentence)
            else
                humid_tend_lsc[i, j, k] = zero(NF)
                temp_tend_lsc[i, j, k] = zero(NF)
            end
        end
    end

    # 2. Precipitation due to large-scale condensation
    fill!(precip_large_scale, zero(NF))
    for k = k_lsc:nlev
        Δpₖ = pres_grid * σ_levels_thick[k]  # Formula 4
        for j = 1:nlat, i = 1:nlon
            precip_large_scale[i, j] += -1 / gravity * Δpₖ[i, j] * humid_tend_lsc[i, j, k]  # Formula 25
        end
    end
end

"""
    get_saturation_vapour_pressure!(column::ColumnVariables,
                                    model::PrimitiveEquationModel)

Compute the saturation vapour pressure as a function of temperature using the
August-Roche-Magnus formula,

    eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T - Tᵢ)),

where T is in Kelvin and i = 1,2 for saturation with respect to water and ice,
respectively.
"""
function get_saturation_vapour_pressure!(   column::ColumnVariables{NF},
                                            model::PrimitiveEquationModel{NF}
                                            ) where {NF<:AbstractFloat}

    @unpack sat_vap_pres, temp = column
    @unpack e₀, T₀, C₁, C₂, T₁, T₂ = model.parameters.magnus_coefs

    for k in eachlayer(column)
        # change coefficients for water (temp > T₀) or ice (else)
        C₁or₂, T₁or₂ = temp[k] > T₀ ? (C₁,T₁) : (C₂,T₂)
        sat_vap_pres[k] = e₀ * exp(C₁or₂ * (temp[k] - T₀) / (temp[k] - T₁or₂))
    end
end

"""
    get_saturation_specific_humidity!(  column::ColumnVariables{NF},
                                        model::PrimitiveEquationModel{NF}
                                        ) where {NF<:AbstractFloat}

Compute the saturation specific humidity according to the formula,

    0.622 * e / (p - (1 - 0.622) * e),

where e is the saturation vapour pressure, p is the pressure, and 0.622 is the ratio of
the molecular weight of water to dry air.
"""
function get_saturation_specific_humidity!( column::ColumnVariables{NF},
                                            model::PrimitiveEquationModel{NF}
                                            ) where {NF<:AbstractFloat}

    @unpack σ_levels_full = model.geometry
    @unpack sat_spec_humid, sat_vap_pres, pres = column

    mol_ratio = convert(NF, 0.622)
    one_minus_mol_ratio = convert(NF, 1 - mol_ratio)

    for k in eachlayer(column)
        pres_at_k = pres*σ_levels_full[k]       # pressure in layer k
        sat_spec_humid[k] = mol_ratio*sat_vap_pres[k] / (pres_at_k - one_minus_mol_ratio*sat_vap_pres[k])
    end
end