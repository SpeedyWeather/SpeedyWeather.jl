abstract type AbstractShortwaveTransmissivity <: AbstractShortwave end

export TransparentShortwaveTransmissivity
TransparentShortwaveTransmissivity(SG::SpectralGrid) = ConstantShortwaveTransmissivity(SG, transmissivity = 1)

export ConstantShortwaveTransmissivity

"""Constant transmissivity for the shortwave, weighted by pressure thickness per layer. ($TYPEDFIELDS)"""
@parameterized @kwdef struct ConstantShortwaveTransmissivity{NF} <: AbstractShortwaveTransmissivity
    "[OPTION] Transmissivity of the whole atmosphere (0 .. 1)"
    @param transmissivity::NF = 0.85 (bounds = 0 .. 1,)
end
Adapt.@adapt_structure ConstantShortwaveTransmissivity
ConstantShortwaveTransmissivity(SG::SpectralGrid; kwargs...) = ConstantShortwaveTransmissivity{SG.NF}(; kwargs...)
initialize!(::ConstantShortwaveTransmissivity, ::AbstractModel) = nothing

@propagate_inbounds function transmissivity!(
        ij,
        vars,
        clouds,
        CST::ConstantShortwaveTransmissivity,
        model,
    )
    t = vars.scratch.grid.a
    nlayers = size(t, 2)

    τ = -log(CST.transmissivity)            # total optical depth of the atmosphere
    coord = model.geometry.vertical_coordinates
    pₛ = vars.parameterizations.surface_pressure[ij]

    for k in 1:nlayers
        Δσₖ = pressure_thickness(k, pₛ, coord) / pₛ
        t[ij, k] = exp(-τ * Δσₖ)            # transmissivity through layer k
    end
    return t
end

export BackgroundShortwaveTransmissivity

"""BackgroundShortwaveTransmissivity <: AbstractShortwaveTransmissivity
$(TYPEDFIELDS)."""
@parameterized @kwdef struct BackgroundShortwaveTransmissivity{NF} <: AbstractShortwaveTransmissivity
    "[OPTION] Zenith correction amplitude (SPEEDY azen) [1]"
    @param zenith_amplitude::NF = 1 (bounds = Nonnegative,)

    "[OPTION] Zenith correction exponent (SPEEDY nzen)"
    @param zenith_exponent::NF = 2 (bounds = Nonnegative,)

    # Weighted visible + near-IR: 0.95*0.033 + 0.05*0.0 = 0.03135 (SPEEDY absdry, fband weights)
    "[OPTION] Absorptivity of dry air [per 10^5 Pa]"
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
    @param absorptivity_cloud_base::NF = 10 (bounds = Nonnegative,)

    # Weighted one-band scaling: 0.95*0.15 = 0.1425 → rounded to 0.14 (SPEEDY abscl2)
    "[OPTION] Maximum cloud absorptivity [per 10^5 Pa]"
    @param absorptivity_cloud_limit::NF = 0.14 (bounds = Nonnegative,)
end

Adapt.@adapt_structure BackgroundShortwaveTransmissivity
BackgroundShortwaveTransmissivity(SG::SpectralGrid; kwargs...) = BackgroundShortwaveTransmissivity{SG.NF}(; kwargs...)
initialize!(::BackgroundShortwaveTransmissivity, ::AbstractModel) = nothing

@propagate_inbounds function transmissivity!(
        ij,
        vars,
        clouds,    # NamedTuple from clouds!
        transmissivity::BackgroundShortwaveTransmissivity,
        model,
    )

    # use scratch array for transmissivity t
    t = vars.scratch.grid.a
    NF = eltype(t)

    (;
        absorptivity_dry_air, absorptivity_aerosol, absorptivity_water_vapor,
        absorptivity_cloud_base, absorptivity_cloud_limit,
    ) = transmissivity
    (; cloud_top, cloud_cover) = clouds

    humid = get_prognostic_step(vars.grid.humidity, model.time_stepping, transmissivity)
    cos_zenith = vars.parameterizations.cos_zenith[ij]
    nlayers = size(t, 2)

    coord = model.geometry.vertical_coordinates
    pₛ = vars.parameterizations.surface_pressure[ij]          # surface pressure [Pa]
    normalized_surface_pressure = pₛ / 100000

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

        # Aerosol factor: use mid-level sigma, squared (aerosol loading increases toward surface)
        aerosol_factor = transmissivity.aerosols ? sigma(k, coord)^2 : zero(NF)

        # Layer absorptivity (all humidity-based parameters are per kg/kg per 10^5 Pa)
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
        Δσₖ = pressure_thickness(k, pₛ, coord) / pₛ
        optical_depth = layer_absorptivity * Δσₖ * normalized_surface_pressure * zenith_factor

        # Transmissivity through layer k
        t[ij, k] = exp(-optical_depth)
    end

    return t
end

export TwoBandShortwaveTransmissivity

"""Shortwave transmissivity for two bands following Fortran SPEEDY (Molteni, 2003):
A visible band absorbed by dry air, aerosols, water vapor and clouds, and a near-infrared
band only absorbed by water vapor. Both are corrected for the slant path through the
atmosphere with a zenith angle correction factor `1 + zenith_amplitude*(1 - cos_zenith)^zenith_exponent`.
Absorptivities are per 10^5 Pa. Humidity-dependent absorptivities per kg/kg (SPEEDY uses g/kg).
$(TYPEDFIELDS)"""
@parameterized @kwdef struct TwoBandShortwaveTransmissivity{NF} <: AbstractShortwaveTransmissivity
    "[OPTION] Zenith correction amplitude (SPEEDY azen) [1]"
    @param zenith_amplitude::NF = 1 (bounds = Nonnegative,)

    "[OPTION] Zenith correction exponent (SPEEDY nzen) [1]"
    @param zenith_exponent::NF = 2 (bounds = Nonnegative,)

    "[OPTION] Absorptivity of dry air, visible band (SPEEDY absdry) [per 10^5 Pa]"
    @param absorptivity_dry_air::NF = 0.033 (bounds = Nonnegative,)

    "[OPTION] Constant aerosol concentration?"
    aerosols::Bool = true

    "[OPTION] Absorptivity of aerosols, visible band, scaled with σ² (SPEEDY absaer) [per 10^5 Pa]"
    @param absorptivity_aerosol::NF = 0.033 (bounds = Nonnegative,)

    "[OPTION] Absorptivity of water vapor, visible band (SPEEDY abswv1) [per kg/kg per 10^5 Pa]"
    @param absorptivity_water_vapor::NF = 22 (bounds = Nonnegative,)

    "[OPTION] Absorptivity of water vapor, near-infrared band (SPEEDY abswv2) [per kg/kg per 10^5 Pa]"
    @param absorptivity_water_vapor_near_infrared::NF = 15000 (bounds = Nonnegative,)

    "[OPTION] Cloud absorptivity per cloud-base humidity, visible band (SPEEDY abscl1) [per kg/kg per 10^5 Pa]"
    @param absorptivity_cloud_base::NF = 15 (bounds = Nonnegative,)

    "[OPTION] Maximum cloud absorptivity, visible band (SPEEDY abscl2) [per 10^5 Pa]"
    @param absorptivity_cloud_limit::NF = 0.15 (bounds = Nonnegative,)
end

Adapt.@adapt_structure TwoBandShortwaveTransmissivity
TwoBandShortwaveTransmissivity(SG::SpectralGrid; kwargs...) = TwoBandShortwaveTransmissivity{SG.NF}(; kwargs...)
initialize!(::TwoBandShortwaveTransmissivity, ::AbstractModel) = nothing

"""$(TYPEDSIGNATURES)
Transmissivities of the visible and near-infrared band for every layer in column `ij`, written
into scratch arrays. Clouds absorb in the visible band from the cloud top down to the layer above
the surface layer. Returns a NamedTuple `(; visible, near_infrared, zenith_factor)`."""
@propagate_inbounds function transmissivity!(
        ij,
        vars,
        clouds,    # NamedTuple from clouds!
        transmissivity::TwoBandShortwaveTransmissivity,
        model,
    )
    # use scratch arrays for transmissivities
    t_vis = vars.scratch.grid.a
    t_nir = vars.scratch.grid.b
    NF = eltype(t_vis)

    (; cloud_top, cloud_cover) = clouds
    humid = get_prognostic_step(vars.grid.humidity, model.time_stepping, transmissivity)
    cos_zenith = vars.parameterizations.cos_zenith[ij]
    nlayers = size(t_vis, 2)

    coord = model.geometry.vertical_coordinates
    pₛ = vars.parameterizations.surface_pressure[ij]    # surface pressure [Pa]

    # Zenith angle correction factor for the slant path
    zenith_factor = 1 + transmissivity.zenith_amplitude * (1 - cos_zenith)^transmissivity.zenith_exponent

    # Cloud absorptivity from humidity at cloud base, taken as the layer above the surface layer
    q_base = humid[ij, max(1, nlayers - 1)]
    cloud_absorptivity = cloud_cover * min(
        transmissivity.absorptivity_cloud_base * q_base,
        transmissivity.absorptivity_cloud_limit,
    )

    for k in 1:nlayers
        q = humid[ij, k]
        aerosol_factor = transmissivity.aerosols ? sigma(k, coord)^2 : zero(NF)

        absorptivity_visible = transmissivity.absorptivity_dry_air +
            transmissivity.absorptivity_aerosol * aerosol_factor +
            transmissivity.absorptivity_water_vapor * q

        # clouds absorb from cloud top down to (excluding) the surface layer
        if cloud_top <= k < nlayers
            absorptivity_visible += cloud_absorptivity
        end

        absorptivity_near_infrared = transmissivity.absorptivity_water_vapor_near_infrared * q

        # pressure thickness normalized by 10^5 Pa as absorptivities are per 10^5 Pa
        Δp = pressure_thickness(k, pₛ, coord) / 100000
        t_vis[ij, k] = exp(-absorptivity_visible * Δp * zenith_factor)
        t_nir[ij, k] = exp(-absorptivity_near_infrared * Δp * zenith_factor)
    end

    return (; visible = t_vis, near_infrared = t_nir, zenith_factor)
end
