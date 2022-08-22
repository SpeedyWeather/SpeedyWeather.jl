"""
    convection!(
        column::ColumnVariables{NF},
        model::PrimitiveEquationModel,
    )
"""
function convection!(
    column::ColumnVariables{NF},
    model::PrimitiveEquationModel,
) where {NF<:AbstractFloat}
    conditional_instability!(column, model)  # Diagnose convection

    if column.cloud_top == column.nlev + 1
        return nothing  # No convection
    end

    @unpack alhc, g = model.constants
    @unpack pres_ref = model.parameters
    @unpack RH_thresh_pbl_cnv, RH_thresh_layers_cnv, pres_thresh_cnv, humid_relax_time_cnv,
    max_entrainment, ratio_secondary_mass_flux = model.constants  # Constants for convection
    @unpack σ_levels_full, σ_levels_thick = model.geometry
    @unpack pres, humid, humid_half, sat_humid, sat_humid_half, dry_static_energy,
    dry_static_energy_half, entrainment_coefficients, cloud_top, nlev = column
    @unpack cloud_base_mass_flux, flux_humid, flux_dry_static_energy, precip_cnv = column

    # Compute the entrainment coefficients for the column. Mass entrainment in layer k is
    # equal to the entrainment coefficient in layer k * the mass flux at the lower boundary.
    # You should pre-compute and reuse these.
    for k = 2:nlev-1
        entrainment_coefficients[k] = max(0.0, σ_levels_full[k] - 0.5)  # Linear or quadratic??
    end
    entrainment_coefficients /= sum(entrainment_coefficients)  # Normalise
    entrainment_coefficients *= max_entrainment                # Scale by maximum entrainment (as a fraction of cloud-base mass flux)

    # 1. Fluxes for the planetary boundary layer
    # Dry static energy and humidity at the upper boundary of the PBL
    dry_static_energy_top_of_pbl = dry_static_energy_half[nlev-1]
    humid_top_of_pbl = min(humid_half[nlev-1], humid[nlev])  # Why?

    # Maximum specific humidity in the PBL
    humid_max_pbl = max(1.01 * humid[nlev], sat_humid[nlev])  # Why?

    # Cloud-base mass flux
    Δp = pres * σ_levels_thick[nlev-1]  # Pressure difference between bottom and top of PBL  # Check index
    # This is the formula as it appears in Fortran Speedy
    mass_flux = pres_ref * Δp / (g * 3600humid_relax_time_cnv) *            # Why?
        min(1.0, 2.0 * (pres - pres_thresh_cnv) / (1 - pres_thresh_cnv)) *  # Why?
        min(5.0, excess_humidity / (humid_max_pbl - humid_top_of_pbl))      # Why?

    # Commented version below is the formula as it appears in the documentation.
    # humid_threshold_pbl = RH_thresh_pbl_cnv * sat_humid[nlev]  # Humidity threshold for convection in the PBL
    # mass_flux_pbl = Δp / (g * 3600humid_relax_time_cnv) * (humid[nlev] - humid_threshold_pbl) /
    #     (sat_humid[nlev] - humid_half[nlev-1])
    column.cloud_base_mass_flux = mass_flux

    # Upward fluxes at upper boundary
    flux_up_humid = mass_flux * humid_max_pbl
    flux_up_dry_static_energy = mass_flux * dry_static_energy[nlev]

    # Downward fluxes at upper boundary
    flux_down_humid = mass_flux * humid_top_of_pbl
    flux_down_dry_static_energy = mass_flux * dry_static_energy_top_of_pbl

    # Net flux
    flux_dry_static_energy[nlev] = flux_down_dry_static_energy - flux_up_dry_static_energy
    flux_humid[nlev] = flux_down_humid - flux_up_humid

    # 2. Fluxes for intermediate layers
    for k in (nlev-1):-1:(cloud_top+1)
        # Fluxes at lower boundary
        flux_dry_static_energy[k] = flux_up_dry_static_energy - flux_down_dry_static_energy  # Change of sign??
        flux_humid[k] = flux_up_humid - flux_down_humid

        # Mass entrainment
        mass_entrainment = entrainment_coefficients[k] * pres * cloud_base_mass_flux  # Why multiply by pres?
        mass_flux += mass_entrainment

        # Upward fluxes at upper boundary
        flux_up_dry_static_energy += mass_entrainment * dry_static_energy[k]
        flux_up_humid += mass_entrainment * humid[k]

        # Downward fluxes at upper boundary
        flux_down_dry_static_energy = mass_flux * dry_static_energy_half[k]
        flux_down_humid = mass_flux * humid_half[k]

        # Net flux of dry static energy and moisture
        flux_dry_static_energy[k] += flux_down_dry_static_energy - flux_up_dry_static_energy
        flux_humid[k] = flux_down_humid - flux_up_humid

        # Secondary moisture flux
        Δhumid = RH_thresh_layers_cnv * sat_humid[k] - humid[k]
        if Δhumid > 0.0
            Δflux_humid = ratio_secondary_mass_flux * cloud_base_mass_flux * Δhumid
            flux_humid[k] += Δflux_humid
            flux_humid[nlev] -= Δflux_humid
        end
    end

    # 3. Fluxes for top-of-convection layer
    # Flux of convective precipitation
    column.precip_cnv = max(flux_up_humid - mass_flux * sat_humid_half[cloud_top], 0.0)

    # Net flux of dry static energy and moisture
    flux_dry_static_energy[cloud_top] = flux_up_dry_static_energy - flux_down_dry_static_energy + alhc * precip_cnv
    flux_humid[cloud_top] = flux_up_humid - flux_down_humid - precip_cnv

    return nothing
end

"""
    conditional_instability!(
        column::ColumnVariables{NF},
        model::PrimitiveEquationModel,
    )
"""
function conditional_instability!(
    column::ColumnVariables{NF},
    model::PrimitiveEquationModel,
) where {NF<:AbstractFloat}
    @unpack pres_thresh_cnv, RH_thresh_PBL_cnv = model.constants
    @unpack alhc = model.parameters

    @unpack nlev = column
    @unpack humid, pres = column
    @unpack sat_humid,
    dry_static_energy,
    moist_static_energy,
    sat_moist_static_energy,
    sat_moist_static_energy_half = column

    if pres > pres_thresh_cnv
        # Saturation (or super-saturated) moist static energy in the PBL
        sat_moist_static_energy_pbl =
            max(moist_static_energy[nlev], sat_moist_static_energy[nlev])

        # Humidity threshold for convection, defined in the PBL and one level above
        humid_threshold_pbl = RH_thresh_PBL_cnv * sat_humid[nlev]
        humid_threshold_above_pbl = RH_thresh_PBL_cnv * sat_humid[nlev-1]

        for k = 3:nlev-3
            # Condition 1: Conditional instability (MSS in PBL < MSS at this half-level)
            if sat_moist_static_energy_pbl > sat_moist_static_energy_half[k]
                # Condition 2a: Gradient of actual moist static energy between lower and upper troposphere
                if moist_static_energy[nlev-1] > sat_moist_static_energy_half[k]
                    column.cloud_top = k
                    column.excess_humidity = max(
                        humid[nlev] - humid_threshold_pbl,
                        (moist_static_energy[nlev] - sat_moist_static_energy_half[k]) /
                        alhc,
                    )
                    break
                    # Condition 2b: Humidity exceeds threshold in both PBL and one layer above
                elseif (humid[nlev] > humid_threshold_pbl) &&
                       (humid[nlev-1] > humid_threshold_above_pbl)
                    column.cloud_top = k
                    column.excess_humidity = humid[nlev] - humid_threshold_pbl
                    break
                end
            end
        end
    end
end
