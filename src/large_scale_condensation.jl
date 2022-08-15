function large_scale_condensation!( column::ColumnVariables{NF},
                                    model::PrimitiveEquationModel,
                                    ) where {NF<:AbstractFloat}

    get_saturation_vapour_pressure!(column, model)
    get_saturation_specific_humidity!(column, model)
    
    @unpack gravity, RH_thresh_max, RH_thresh_range, RH_thresh_boundary, humid_relax_time = model.constants
    @unpack cp, alhc, k_lsc = model.parameters
    @unpack σ_levels_full, σ_levels_thick = model.geometry

    @unpack humid, pres = column            # prognostic variables: specific humidity, surface pressure
    @unpack temp_tend, humid_tend = column  # tendencies to write into

    # 1. Tendencies of humidity and temperature due to large-scale condensation
    # skip the planetary boundary layer, k_lsc is lowest level for precip 
    for k in eachlayer(column)[k_lsc:end]
        
        # Relative humidity threshold for condensation (Formula 24)
        σₖ² = σ_levels_full[k]^2
        RH_threshold = RH_thresh_max + RH_thresh_range * (σₖ² - 1)
        if k == nlev
            RH_threshold = max(RH_threshold, RH_thresh_boundary)
        end

        # Impose a maximum heating rate to avoid grid-point storm instability
        # This formula does not appear in the original Speedy documentation
        humid_tend_max = 10σₖ² / 3600humid_relax_time
        humid_threshold = RH_threshold * sat_spec_humidity[k]  # Specific humidity threshold for condensation

        if humid[k] > humid_threshold
            # accumulate in tendencies (nothing is added if humidity not above threshold)
            humid_tend[k] += -(humid[k] - humid_threshold) / humid_relax_time       # Formula 22
            temp_tend[k] += -alhc / cp * min(humid_tend[k], humid_tend_max*pres)    # Formula 23
            column.cloud_top = min(column.cloud_top, k)                             # Page 7 (last sentence)
        end
    end

    # 2. Precipitation due to large-scale condensation
    column.precip_large_scale = 0
    for k in eachlayer(column)[k_lsc:end]
        Δpₖ = pres*σ_levels_thick[k]                                        # Formula 4
        column.precip_large_scale += -1 / gravity * Δpₖ * humid_tend[k]     # Formula 25
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