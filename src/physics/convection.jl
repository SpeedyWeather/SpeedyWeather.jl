"""
    diagnose_convection!(
        column::ColumnVariables{NF},
        model::PrimitiveEquation,
    )

Check whether the convection scheme should be activated in the given atmospheric column.

1. A conditional instability exists when the saturation moist energy (MSS) decreases with
height, that is, there exists an atmospheric level k such that,

    MSS(N) > MSS(k+h)

where N is the planetary boundary layer (PBL) and k+h is the half-level at the lower
boundary of the full level k.

2. When a conditional instability exists, the convection scheme is activated when, either,

    a. the actual moist static energy (MSE) at level N-1 (directly above the PBL) is greater
       than the saturation moist static energy at some half-level k+h,

            MSE(N-1) > MSS(k+h)

    b. the humidity in both the PBL and one layer above exceeds a prescribed threshold,

            Q(N)   > RH_cnv * Qˢᵃᵗ(N)
            Q(N-1) > RH_cnv * Qˢᵃᵗ(N-1)

The top-of-convection (TCN) layer, or cloud-top, is the largest value of k for which
condition 1 is satisfied. The cloud-top layer may be subsequently adjusted upwards by the
large-scale condensation parameterization, which is executed after this one.
"""
function diagnose_convection!(
    column::ColumnVariables{NF},
    model::PrimitiveEquation,
) where {NF<:AbstractFloat}
    (; alhc,pres_ref ) = model.parameters
    (; pres_thresh_cnv, RH_thresh_pbl_cnv ) = model.constants
    (; nlev ) = column
    (; humid, pres, sat_humid, dry_static_energy, moist_static_energy,
    sat_moist_static_energy, sat_moist_static_energy_half) = column

    if pres[end] > pres_thresh_cnv
        # First we pre-compute some values which we will need inside the loop
        # 1. Saturation (or super-saturated) moist static energy in the PBL
        sat_moist_static_energy_pbl =
            max(moist_static_energy[nlev], sat_moist_static_energy[nlev])

        # 2. Minimum of moist static energy in the lowest two levels
        moist_static_energy_lower_trop =
            min(moist_static_energy[nlev], moist_static_energy[nlev-1])

        # 3. Humidity threshold for convection, defined in the PBL and one level above
        humid_threshold_pbl = RH_thresh_pbl_cnv * sat_humid[nlev]
        humid_threshold_above_pbl = RH_thresh_pbl_cnv * sat_humid[nlev-1]

        # The range of this loop requires clarification, but in its current form it means
        # that the top-of-convection level may be any tropospheric level, excluding the two
        # layers directly above the PBL.
        for k = (nlev-3):-1:3
            # Condition 1: Conditional instability (MSS in PBL < MSS at this half-level)
            if sat_moist_static_energy_pbl > sat_moist_static_energy_half[k]
                column.conditional_instability = true
                column.cloud_top = k
            end

            # Condition 2a: Gradient of actual moist static energy between lower and upper troposphere
            if moist_static_energy_lower_trop > sat_moist_static_energy_half[k]
                column.activate_convection = true
                column.excess_humidity = max(
                    humid[nlev] - humid_threshold_pbl,
                    (moist_static_energy[nlev] - sat_moist_static_energy_half[k]) / alhc,
                )
            end
        end

        if column.conditional_instability && column.activate_convection
            return nothing  # Condition for convection already satisfied
        end

        # Condition 2b: Humidity exceeds threshold in both PBL and one layer above
        if column.conditional_instability &&
           (humid[nlev] > humid_threshold_pbl) &&
           (humid[nlev-1] > humid_threshold_above_pbl)
            column.activate_convection = true
            column.excess_humidity = humid[nlev] - humid_threshold_pbl
        end
    end
    return nothing
end

"""
    convection!(
        column::ColumnVariables{NF},
        model::PrimitiveEquation,
    )

Compute fluxes and precipitation due to convection in the given atmospheric column.

The scheme computes fluxes of mass, humidity and dry static energy. A part of the upward
moisture flux at the lower boundary of the cloud-top (TCN) layer is converted into
convective precipitation.

For full details of the scheme see: http://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf
"""
function convection!(
    column::ColumnVariables{NF},
    model::PrimitiveEquation,
) where {NF<:AbstractFloat}
    diagnose_convection!(column, model)  # Diagnose convection

    if !(column.conditional_instability && column.activate_convection)
        return nothing  # No convection
    end

    (; gravity ) = model.constants
    (; alhc, pres_ref ) = model.parameters
    (; σ_levels_full, σ_levels_thick ) = model.geometry
    # Constants for convection
    (;RH_thresh_pbl_cnv, RH_thresh_trop_cnv, pres_thresh_cnv, humid_relax_time_cnv,
    max_entrainment, ratio_secondary_mass_flux) = model.constants
    # Column variables for calculating fluxes due to convection
    (;pres, humid, humid_half, sat_humid, sat_humid_half, dry_static_energy,
    dry_static_energy_half, entrainment_profile, cloud_top, excess_humidity,
    nlev) = column
    # Quantities calculated by this parameterization
    (;cloud_base_mass_flux, net_flux_humid, net_flux_dry_static_energy,
    precip_convection) = column

    # 1. Fluxes in the PBL
    humid_top_of_pbl = min(humid_half[nlev-1], humid[nlev])   # Humidity at the upper boundary of the PBL
    max_humid_pbl = max(NF(1.01) * humid[nlev], sat_humid[nlev])  # Maximum specific humidity in the PBL

    # Cloud-base mass flux
    pₛ = pres[end]                               # surface pressure
    Δp = pres_ref * pₛ * σ_levels_thick[nlev]  # Pressure difference between bottom and top of PBL
    mass_flux =
        Δp / (gravity * 3600humid_relax_time_cnv) *
        min(1, 2 * (pₛ - pres_thresh_cnv) / (1 - pres_thresh_cnv)) *
        min(5, excess_humidity / (max_humid_pbl - humid_top_of_pbl))
    column.cloud_base_mass_flux = mass_flux

    # Upward fluxes at upper boundary
    flux_up_humid = mass_flux * max_humid_pbl
    flux_up_dry_static_energy = mass_flux * dry_static_energy[nlev]

    # Downward fluxes at upper boundary
    flux_down_humid = mass_flux * humid_top_of_pbl
    flux_down_dry_static_energy = mass_flux * dry_static_energy_half[nlev-1]

    # Net flux
    net_flux_dry_static_energy[nlev] = flux_down_dry_static_energy - flux_up_dry_static_energy
    net_flux_humid[nlev] = flux_down_humid - flux_up_humid

    # 2. Fluxes for intermediate layers
    for k = (nlev-1):-1:(cloud_top+1)
        # Fluxes at lower boundary
        net_flux_dry_static_energy[k] = flux_up_dry_static_energy - flux_down_dry_static_energy
        net_flux_humid[k] = flux_up_humid - flux_down_humid

        # Mass entrainment
        mass_entrainment = entrainment_profile[k] * pₛ * cloud_base_mass_flux  # Why multiply by pres here?
        mass_flux += mass_entrainment

        # Upward fluxes at upper boundary
        flux_up_dry_static_energy += mass_entrainment * dry_static_energy[k]
        flux_up_humid += mass_entrainment * humid[k]

        # Downward fluxes at upper boundary
        flux_down_dry_static_energy = mass_flux * dry_static_energy_half[k-1]
        flux_down_humid = mass_flux * humid_half[k-1]

        # Net flux of dry static energy and moisture
        net_flux_dry_static_energy[k] += flux_down_dry_static_energy - flux_up_dry_static_energy
        net_flux_humid[k] = flux_down_humid - flux_up_humid

        # Secondary moisture flux representing shallower, non-precipitating convective systems
        # Occurs when RH in an intermediate layer falls below a threshold
        Δhumid = RH_thresh_trop_cnv * sat_humid[k] - humid[k]
        if Δhumid > 0
            Δflux_humid = ratio_secondary_mass_flux * cloud_base_mass_flux * Δhumid
            net_flux_humid[k] += Δflux_humid
            net_flux_humid[nlev] -= Δflux_humid
        end
    end

    # 3. Fluxes for top-of-convection layer
    # Flux of convective precipitation
    column.precip_convection = max(flux_up_humid - mass_flux * sat_humid_half[cloud_top], 0)

    # Net flux of dry static energy and moisture
    net_flux_dry_static_energy[cloud_top] =
        flux_up_dry_static_energy - flux_down_dry_static_energy + alhc * precip_convection
    net_flux_humid[cloud_top] = flux_up_humid - flux_down_humid - precip_convection

    return nothing
end


    # # Compute the entrainment coefficients for the convection parameterization.
    # (;max_entrainment) = P
    # entrainment_profile = zeros(nlev)
    # for k = 2:nlev-1
    #     entrainment_profile[k] = max(0, (σ_levels_full[k] - 0.5)^2)
    # end

    # # profile as fraction of cloud-base mass flux
    # entrainment_profile /= sum(entrainment_profile)  # Normalise
    # entrainment_profile *= max_entrainment           # fraction of max entrainment

    # # PARAMETRIZATIONS
    # # Large-scale condensation (occurs when relative humidity exceeds a given threshold)
    # RH_thresh_pbl_lsc::NF    # Relative humidity threshold for LSC in PBL
    # RH_thresh_range_lsc::NF  # Vertical range of relative humidity threshold
    # RH_thresh_max_lsc ::NF   # Maximum relative humidity threshold
    # humid_relax_time_lsc::NF # Relaxation time for humidity (hours)

    # # Convection
    # pres_thresh_cnv::NF            # Minimum (normalised) surface pressure for the occurrence of convection
    # RH_thresh_pbl_cnv::NF          # Relative humidity threshold for convection in PBL
    # RH_thresh_trop_cnv::NF         # Relative humidity threshold for convection in the troposphere
    # humid_relax_time_cnv::NF       # Relaxation time for PBL humidity (hours)
    # max_entrainment::NF            # Maximum entrainment as a fraction of cloud-base mass flux
    # ratio_secondary_mass_flux::NF  # Ratio between secondary and primary mass flux at cloud-base


    # "For computing saturation vapour pressure"
    # magnus_coefs::Coefficients = MagnusCoefs{NF}()

    # # Large-Scale Condensation (from table B10)
    # "Index of atmospheric level at which large-scale condensation begins"
    # k_lsc::Int = 2

    # "Relative humidity threshold for boundary layer"
    # RH_thresh_pbl_lsc::Float64 = 0.95

    # "Vertical range of relative humidity threshold"
    # RH_thresh_range_lsc::Float64 = 0.1

    # "Maximum relative humidity threshold"
    # RH_thresh_max_lsc::Float64 = 0.9

    # "Relaxation time for humidity (hours)"
    # humid_relax_time_lsc::Float64 = 4.0

    # # Convection
    # "Minimum (normalised) surface pressure for the occurrence of convection"
    # pres_thresh_cnv::Float64 = 0.8

    # "Relative humidity threshold for convection in PBL"
    # RH_thresh_pbl_cnv::Float64 = 0.9

    # "Relative humidity threshold for convection in the troposphere"
    # RH_thresh_trop_cnv::Float64 = 0.7

    # "Relaxation time for PBL humidity (hours)"
    # humid_relax_time_cnv::Float64 = 6.0

    # "Maximum entrainment as a fraction of cloud-base mass flux"
    # max_entrainment::Float64 = 0.5

    # "Ratio between secondary and primary mass flux at cloud-base"
    # ratio_secondary_mass_flux::Float64 = 0.8
