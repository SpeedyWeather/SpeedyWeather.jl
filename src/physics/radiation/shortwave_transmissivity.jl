abstract type AbstractShortwaveTransmissivity <: AbstractShortwave end

# function barrier to dispatch to type of model.transmissivity
function transmissivity!(column::ColumnVariables, model::AbstractModel)
    return transmissivity!(column, model.transmissivity, model)
end

export TransparentShortwaveTransmissivity
struct TransparentShortwaveTransmissivity{NF} <: AbstractShortwaveTransmissivity end
TransparentShortwaveTransmissivity(SG::SpectralGrid) = TransparentShortwaveTransmissivity{SG.NF}()
initialize!(::TransparentShortwaveTransmissivity, ::AbstractModel) = nothing
function transmissivity!(column::ColumnVariables, clouds, ::TransparentShortwaveTransmissivity, ::AbstractModel, band::Int = 1)
    column.transmissivity_longwave .= 1
    return view(column.transmissivity_shortwave, :, band)
end

export BackgroundShortwaveTransmissivity
@parameterized @kwdef struct BackgroundShortwaveTransmissivity{NF} <: AbstractShortwaveTransmissivity
    "[OPTION] Zenith correction amplitude (SPEEDY azen) [1]"
    @param zenith_amplitude::NF = 1 (bounds=Nonnegative,)

    "[OPTION] Zenith correction exponent (SPEEDY nzen)"
    @param zenith_exponent::NF = 2 (bounds=Nonnegative,)

    "[OPTION] Absorptivity of dry air [per 10^5 Pa]"
    # Weighted visible + near-IR: 0.95*0.033 + 0.05*0.0 = 0.03135 (SPEEDY absdry, fband weights)
    @param absorptivity_dry_air::NF = 0.03135 (bounds=0..1,)

    "[OPTION] Constant aerosol concentration?"
    aerosols::Bool = true

    "[OPTION] Absorptivity of aerosols [per 10^5 Pa]"
    # Weighted visible + near-IR: 0.95*0.033 + 0.05*0.0 = 0.03135 (SPEEDY absaer, fband weights)
    @param absorptivity_aerosol::NF = 0.03135 (bounds=0..1,)

    "[OPTION] Absorptivity of water vapor [per kg/kg per 10^5 Pa]"
    # Weighted visible + near-IR: 0.95*0.022 + 0.05*15.0*0.2 = 0.171 per g/kg → 1.7e-4 per kg/kg (SPEEDY abswv1, abswv2)
    @param absorptivity_water_vapor::NF = 0.00017 (bounds=0..1,)
    "[OPTION] Base cloud absorptivity [per kg/kg per 10^5 Pa]"
    # Weighted visible band: 0.95*0.015 = 0.014 per g/kg → 1.4e-5 per kg/kg (SPEEDY abscl1)
    @param absorptivity_cloud_base::NF = 0.000014 (bounds=0..1,)

    "[OPTION] Maximum cloud absorptivity [per 10^5 Pa]"
    # Weighted one-band scaling: 0.95*0.15 = 0.1425 → rounded to 0.14 (SPEEDY abscl2)
    @param absorptivity_cloud_limit::NF = 0.14 (bounds=0..1,)
end

BackgroundShortwaveTransmissivity(SG::SpectralGrid; kwargs...) = BackgroundShortwaveTransmissivity{SG.NF}(; kwargs...)
initialize!(::BackgroundShortwaveTransmissivity, ::AbstractModel) = nothing

function transmissivity!(
        column::ColumnVariables,
        clouds,    # NamedTuple from clouds!
        transmissivity::BackgroundShortwaveTransmissivity,
        model::AbstractModel,
        band::Int = 1,  # Which spectral band to compute
    )
    t = view(column.transmissivity_shortwave, :, band)

    (;
        absorptivity_dry_air, absorptivity_aerosol, absorptivity_water_vapor,
        absorptivity_cloud_base, absorptivity_cloud_limit,
    ) = transmissivity
    (; cloud_top, cloud_cover) = clouds
    (; humid, cos_zenith, nlayers) = column

    sigma_levels = model.geometry.σ_levels_half
    sigma_levels_full = model.geometry.σ_levels_full
    surface_pressure = column.pres[end]  # This is in Pa

    # Zenith angle correction factor
    azen = transmissivity.zenith_amplitude
    nzen = transmissivity.zenith_exponent

    # Zenith angle correction to (downward) absorptivity
    zenit_factor = 1 + azen * (1 - cos_zenith)^nzen

    # Cloud absorption term based on cloud base humidity (SPEEDY logic)
    q_base = nlayers > 1 ? humid[nlayers - 1] : humid[nlayers]
    cloud_absorptivity_term = min(
        absorptivity_cloud_base * q_base,
        absorptivity_cloud_limit
    )

    for k in 1:nlayers
        q = humid[k]

        # Aerosol factor: use mid-level sigma, squared
        aerosol_factor = transmissivity.aerosols ? sigma_levels_full[k]^2 : 0

        # Layer absorptivity (all humidity-based parameters are per kg/kg per 10^5 Pa)
        # Aerosol loading increases toward surface (proportional to σ²).
        layer_absorptivity = (
            absorptivity_dry_air +
                absorptivity_aerosol * aerosol_factor +
                absorptivity_water_vapor * q
        )

        # Add cloud absorption below the FINAL cloud top
        if k >= cloud_top
            layer_absorptivity += cloud_absorptivity_term * cloud_cover
        end

        # Compute differential optical depth with zenith correction
        # CRITICAL: Normalize pressure to 10^5 Pa since absorptivities are per 10^5 Pa
        delta_sigma = sigma_levels[k + 1] - sigma_levels[k]
        normalized_pressure = surface_pressure / 100000   # Convert Pa to units of 10^5 Pa
        optical_depth = layer_absorptivity * delta_sigma * normalized_pressure * zenit_factor

        # Transmissivity through layer k
        t[k] = exp(-optical_depth)
    end

    return t
end
