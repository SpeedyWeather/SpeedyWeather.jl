abstract type AbstractShortwaveTransmissivity <: AbstractShortwave end

export TransparentShortwaveTransmissivity
struct TransparentShortwaveTransmissivity <: AbstractShortwaveTransmissivity end
Adapt.@adapt_structure TransparentShortwaveTransmissivity
TransparentShortwaveTransmissivity(SG::SpectralGrid) = TransparentShortwaveTransmissivity()
initialize!(::TransparentShortwaveTransmissivity, ::AbstractModel) = nothing

@propagate_inbounds function transmissivity!(
        ij,
        diagn,
        progn,
        clouds,
        ::TransparentShortwaveTransmissivity,
        model,
    )
    t = diagn.dynamics.a_grid
    nlayers = size(t, 2)
    for k in 1:nlayers
        t[ij, k] = one(eltype(t))
    end
    return t
end

export BackgroundShortwaveTransmissivity

"""    BackgroundShortwaveTransmissivity <: AbstractShortwaveTransmissivity
$(TYPEDFIELDS)."""
@parameterized @kwdef struct BackgroundShortwaveTransmissivity{NF} <: AbstractShortwaveTransmissivity
    "[OPTION] Zenith correction amplitude (SPEEDY azen) [1]"
    @param zenith_amplitude::NF = 1 (bounds = Nonnegative,)

    "[OPTION] Zenith correction exponent (SPEEDY nzen)"
    @param zenith_exponent::NF = 2 (bounds = Nonnegative,)

    "[OPTION] Absorptivity of dry air [per 10^5 Pa]"
    # Weighted visible + near-IR: 0.95*0.033 + 0.05*0.0 = 0.03135 (SPEEDY absdry, fband weights)
    @param absorptivity_dry_air::NF = 0.03135 (bounds = Nonnegative,)

    "[OPTION] Constant aerosol concentration?"
    aerosols::Bool = true

    "[OPTION] Absorptivity of aerosols [per 10^5 Pa]"
    # Weighted visible + near-IR: 0.95*0.033 + 0.05*0.0 = 0.03135 (SPEEDY absaer, fband weights)
    @param absorptivity_aerosol::NF = 0.03135 (bounds = Nonnegative,)

    # Weighted visible + near-IR: 0.95*0.022 + 0.05*15.0*0.2 = 0.171 per g/kg → 1.7e-4 per kg/kg (SPEEDY abswv1, abswv2)
    # Value chosen following PR #974 for a ~75W/m^2 shortwave absorption target
    "[OPTION] Absorptivity of water vapor [per kg/kg per 10^5 Pa]"
    @param absorptivity_water_vapor::NF = 75 (bounds = Nonnegative,)
    
    # Weighted visible band: 0.95*0.015 = 0.014 per g/kg → 1.4e-5 per kg/kg (SPEEDY abscl1)
    "[OPTION] Base cloud absorptivity [per kg/kg per 10^5 Pa]"
    @param absorptivity_cloud_base::NF = 15 (bounds = Nonnegative,)

    "[OPTION] Maximum cloud absorptivity [per 10^5 Pa]"
    # Weighted one-band scaling: 0.95*0.15 = 0.1425 → rounded to 0.14 (SPEEDY abscl2)
    @param absorptivity_cloud_limit::NF = 0.14 (bounds = Nonnegative,)
end

Adapt.@adapt_structure BackgroundShortwaveTransmissivity
BackgroundShortwaveTransmissivity(SG::SpectralGrid; kwargs...) = BackgroundShortwaveTransmissivity{SG.NF}(; kwargs...)
initialize!(::BackgroundShortwaveTransmissivity, ::AbstractModel) = nothing

@propagate_inbounds function transmissivity!(
        ij,
        diagn,
        progn,
        clouds,    # NamedTuple from clouds!
        transmissivity::BackgroundShortwaveTransmissivity,
        model,
    )

    # use scratch array for transmissivity t
    t = diagn.dynamics.a_grid
    NF = eltype(t)

    (;
        absorptivity_dry_air, absorptivity_aerosol, absorptivity_water_vapor,
        absorptivity_cloud_base, absorptivity_cloud_limit,
    ) = transmissivity
    (; cloud_top, cloud_cover) = clouds

    humid = diagn.grid.humid_grid_prev
    cos_zenith = diagn.physics.cos_zenith[ij]
    nlayers = size(t, 2)

    sigma_levels = model.geometry.σ_levels_half
    sigma_levels_full = model.geometry.σ_levels_full
    normalized_surface_pressure = diagn.grid.pres_grid_prev[ij] / 100000

    # Zenith angle correction factor
    azen = transmissivity.zenith_amplitude
    nzen = transmissivity.zenith_exponent
    zenith_factor = 1 + azen * (1 - cos_zenith)^nzen

    # Cloud absorption term based on cloud base humidity (SPEEDY logic)
    q_base = nlayers > 1 ? humid[ij, nlayers - 1] : humid[ij, nlayers]
    cloud_absorptivity_term = min(
        absorptivity_cloud_base * q_base,
        absorptivity_cloud_limit
    )

    for k in 1:nlayers
        q = humid[ij, k]

        # Aerosol factor: use mid-level sigma, squared
        aerosol_factor = transmissivity.aerosols ? sigma_levels_full[k]^2 : zero(NF)

        # Layer absorptivity (all humidity-based parameters are per kg/kg per 10^5 Pa)
        # Aerosol loading increases toward surface (proportional to σ²).
        layer_absorptivity = (
            absorptivity_dry_air +
                absorptivity_aerosol * aerosol_factor +
                absorptivity_water_vapor * q
        )

        # Add cloud absorption below the final cloud top
        if k >= cloud_top
            layer_absorptivity += cloud_absorptivity_term * cloud_cover
        end

        # Compute differential optical depth with zenith correction
        # Normalize pressure to 1e5 Pa since absorptivities are per 1e5 Pa
        delta_sigma = sigma_levels[k + 1] - sigma_levels[k]
        optical_depth = layer_absorptivity * delta_sigma * normalized_surface_pressure * zenith_factor

        # Transmissivity through layer k
        t[ij, k] = exp(-optical_depth)
    end

    return t
end
