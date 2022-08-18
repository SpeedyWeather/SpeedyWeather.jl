function conditional_instability!(
    column::ColumnVariables{NF},
    model::PrimitiveEquationModel,
) where {NF<:AbstractFloat}
    @unpack pres_thresh_cnv, RH_thresh_PBL_cnv = model.constants
    @unpack alhc = model.parameters

    @unpack nlev = column
    @unpack humid, pres = column
    @unpack sat_humid, dry_static_energy, moist_static_energy, sat_moist_static_energy,
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
                        (moist_static_energy[nlev] - sat_moist_static_energy_half[k]) / alhc,
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

function convection!(
    column::ColumnVariables{NF},
    model::PrimitiveEquationModel,
) where {NF<:AbstractFloat} end
