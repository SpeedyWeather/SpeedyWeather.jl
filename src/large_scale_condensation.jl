function large_scale_condensation!( column::ColumnVariables{NF},
                                    model::PrimitiveEquationModel,
                                    ) where {NF<:AbstractFloat}

    get_saturation_vapour_pressure!(column, model)
    get_saturation_specific_humidity!(column, model)

    @unpack gravity, RH_thresh_max_lsc, RH_thresh_range_lsc, RH_thresh_PBL_lsc, humid_relax_time_lsc = model.constants
    @unpack cp, alhc, n_stratosphere_levels = model.parameters
    @unpack σ_levels_full, σ_levels_thick = model.geometry

    @unpack humid, pres = column            # prognostic variables: specific humidity, surface pressure
    @unpack temp_tend, humid_tend = column  # tendencies to write into
    @unpack sat_humid = column              # intermediate variable
    @unpack nlev = column

    # 1. Tendencies of humidity and temperature due to large-scale condensation
    for k in eachlayer(column)[n_stratosphere_levels:end]           # top to bottom, skip stratospheric levels

        # Relative humidity threshold for condensation (Formula 24)
        σₖ² = σ_levels_full[k]^2
        RH_threshold = RH_thresh_max_lsc + RH_thresh_range_lsc * (σₖ² - 1)
        if k == nlev
            RH_threshold = max(RH_threshold, RH_thresh_PBL_lsc)
        end

        # Impose a maximum heating rate to avoid grid-point storm instability
        # This formula does not appear in the original Speedy documentation
        humid_tend_max = 10σₖ² / 3600humid_relax_time_lsc
        humid_threshold = RH_threshold * sat_humid[k]               # Specific humidity threshold for condensation

        if humid[k] > humid_threshold
            # accumulate in tendencies (nothing is added if humidity not above threshold)
            humid_tend[k] += -(humid[k] - humid_threshold) / humid_relax_time_lsc   # Formula 22
            temp_tend[k] += -alhc / cp * min(humid_tend[k], humid_tend_max*pres)    # Formula 23

            # highest model level where condensation occurs (=cloud top), initialised with nlev+1 for min here
            column.cloud_top = min(column.cloud_top, k)                             # Page 7 (last sentence)
        end
    end

    # 2. Precipitation due to large-scale condensation
    for k in eachlayer(column)[n_stratosphere_levels:end]                   # top to bottom, skip stratosphere
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
        C, T = temp[k] > T₀ ? (C₁,T₁) : (C₂,T₂)
        sat_vap_pres[k] = e₀ * exp(C * (temp[k] - T₀) / (temp[k] - T))
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
    @unpack sat_humid, sat_vap_pres, pres = column

    mol_ratio = convert(NF, 0.622)
    one_minus_mol_ratio = convert(NF, 1-mol_ratio)

    for k in eachlayer(column)
        pₖ = pres*σ_levels_full[k]       # pressure in layer k
        sat_humid[k] = mol_ratio*sat_vap_pres[k] / (pₖ - one_minus_mol_ratio*sat_vap_pres[k])
    end
end
