abstract type AbstractShortwaveTransmittance <: AbstractShortwave end

# function barrier to dispatch to type of model.transmittance
function transmittance!(column::ColumnVariables, model::AbstractModel)
    transmittance!(column, model.transmittance, model)
end

export TransparentShortwaveTransmittance
struct TransparentShortwaveTransmittance{NF} <: AbstractShortwaveTransmittance end
TransparentShortwaveTransmittance(SG::SpectralGrid) = TransparentShortwaveTransmittance{SG.NF}()
initialize!(::TransparentShortwaveTransmittance, ::AbstractModel) = nothing
function transmittance!(column::ColumnVariables, clouds, ::TransparentShortwaveTransmittance, ::AbstractModel)
    column.transmittance_longwave .= 1
    return view(column.transmittance_shortwave, :, 1)
end

@kwdef struct BackgroundShortwaveTransmittance{NF} <: AbstractShortwaveTransmittance
    "[OPTION] Zenith correction amplitude (SPEEDY azen) [1]"
    zenith_amplitude::NF = 1

    "[OPTION] Zenith correction exponent (SPEEDY nzen)"
    zenith_exponent::NF = 2

    "[OPTION] Absorptivity of dry air [per 10^5 Pa]"
    # Weighted visible + near-IR: 0.95*0.033 + 0.05*0.0 = 0.03135 (SPEEDY absdry, fband weights)
    absorptivity_dry_air::NF = 0.03135

    "[OPTION] Constant aerosol concentration?"
    aerosols::Bool = true

    "[OPTION] Absorptivity of aerosols [per 10^5 Pa]"
    # Weighted visible + near-IR: 0.95*0.033 + 0.05*0.0 = 0.03135 (SPEEDY absaer, fband weights)
    absorptivity_aerosol::NF = 0.03135

    "[OPTION] Absorptivity of water vapor [per kg/kg per 10^5 Pa]"
    # Weighted visible + near-IR: 0.95*0.022 + 0.05*15.0*0.2 = 0.171 per g/kg → 1.7e-4 per kg/kg (SPEEDY abswv1, abswv2)
    absorptivity_water_vapor::NF = 0.00017

    "[OPTION] Base cloud absorptivity [per kg/kg per 10^5 Pa]"
    # Weighted visible band: 0.95*0.015 = 0.014 per g/kg → 1.4e-5 per kg/kg (SPEEDY abscl1)
    absorptivity_cloud_base::NF = 0.000014

    "[OPTION] Maximum cloud absorptivity [per 10^5 Pa]"
    # Weighted one-band scaling: 0.95*0.15 = 0.1425 → rounded to 0.14 (SPEEDY abscl2)
    absorptivity_cloud_limit::NF = 0.14
end

BackgroundShortwaveTransmittance(SG::SpectralGrid; kwargs...) = BackgroundShortwaveTransmittance{SG.NF}(; kwargs...)
initialize!(::BackgroundShortwaveTransmittance, ::AbstractModel) = nothing

function transmittance!(
    column::ColumnVariables,
    clouds,    # NamedTuple from clouds!
    transmittance::BackgroundShortwaveTransmittance,
    model::AbstractModel,
)
    t = view(column.transmittance_shortwave, :, 1)
 
    (; absorptivity_dry_air, absorptivity_aerosol, absorptivity_water_vapor,
    absorptivity_cloud_base, absorptivity_cloud_limit) = transmittance
    (; cloud_top, cloud_cover) = clouds
    (; humid, cos_zenith, nlayers) = column

    sigma_levels = model.geometry.σ_levels_half
    sigma_levels_full = model.geometry.σ_levels_full
    surface_pressure = column.pres[end]  # This is in Pa

    # Zenith angle correction factor 
    azen = transmittance.zenith_amplitude
    nzen = transmittance.zenith_exponent

    # Zenith angle correction to (downward) absorptivity
    zenit_factor = 1 + azen * (1 - cos_zenith)^nzen

    # Cloud absorption term based on cloud base humidity (SPEEDY logic)
    q_base = nlayers > 1 ? humid[nlayers-1] : humid[nlayers]
    cloud_absorptivity_term = min(absorptivity_cloud_base * q_base,
                                absorptivity_cloud_limit)

    for k in 1:nlayers
        q = humid[k]

        # Aerosol factor: use mid-level sigma, squared
        aerosol_factor = transmittance.aerosols ? sigma_levels_full[k]^2 : 0
        
        # Layer absorptivity (all humidity-based parameters are per kg/kg per 10^5 Pa)
        # Aerosol loading increases toward surface (proportional to σ²).
        layer_absorptivity = (absorptivity_dry_air +
                            absorptivity_aerosol * aerosol_factor +
                            absorptivity_water_vapor * q)

        # Add cloud absorption below the FINAL cloud top
        if k >= cloud_top
            layer_absorptivity += cloud_absorptivity_term * cloud_cover
        end

        # Compute differential optical depth with zenith correction
        # CRITICAL: Normalize pressure to 10^5 Pa since absorptivities are per 10^5 Pa
        delta_sigma = sigma_levels[k+1] - sigma_levels[k]
        normalized_pressure = surface_pressure / 100000   # Convert Pa to units of 10^5 Pa
        optical_depth = layer_absorptivity * delta_sigma * normalized_pressure * zenit_factor

        # Transmittance through layer k
        t[k] = exp(-optical_depth)
    end

    return t
end