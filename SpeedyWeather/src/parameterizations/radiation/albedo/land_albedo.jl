abstract type AbstractSnowCover end

export LinearSnowCover, SaturatingSnowCover
"""Linear ramp: snow cover grows with snow depth from 0 to 1 at `snow_depth_scale`."""
struct LinearSnowCover <: AbstractSnowCover end
Adapt.@adapt_structure LinearSnowCover

"""Saturating ramp: snow cover grows with snow depth S as `S/(S+scale)`."""
struct SaturatingSnowCover <: AbstractSnowCover end
Adapt.@adapt_structure SaturatingSnowCover

"""$(TYPEDSIGNATURES) Snow cover fraction for the linear scheme, clamped to 1."""
@inline (::LinearSnowCover)(snow_depth, scale) = min(snow_depth / scale, 1)

"""$(TYPEDSIGNATURES) Snow cover fraction for the saturating scheme."""
@inline (::SaturatingSnowCover)(snow_depth, scale) = snow_depth / (snow_depth + scale)

# LandSnowAlbedo
export LandSnowAlbedo

"""Albedo over land based on bare soil, vegetation (high and low cover) and snow cover.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct LandSnowAlbedo{NF, Scheme <: AbstractSnowCover} <: AbstractAlbedo
    "Albedo of bare land (excluding vegetation) [1]"
    @param albedo_land::NF = 0.4 (bounds = 0 .. 1,)

    "Albedo of high vegetation [1]"
    @param albedo_high_vegetation::NF = 0.15 (bounds = 0 .. 1,)

    "Albedo of low vegetation [1]"
    @param albedo_low_vegetation::NF = 0.2 (bounds = 0 .. 1,)

    "Albedo of snow [1], additive to land"
    @param albedo_snow::NF = 0.4 (bounds = 0 .. 1,)

    "Conversion from snow depth to snow cover [m]"
    @param snow_depth_scale::NF = 0.05 (bounds = Positive,)

    "Snow cover-albedo scheme"
    @param snow_cover::Scheme = SaturatingSnowCover() (group = :snow_cover,)
end

Adapt.@adapt_structure LandSnowAlbedo
LandSnowAlbedo(SG::SpectralGrid; snow_cover = SaturatingSnowCover(), kwargs...) =
    LandSnowAlbedo{SG.NF, typeof(snow_cover)}(; snow_cover, kwargs...)

initialize!(albedo::LandSnowAlbedo, model::PrimitiveEquation) = nothing

@propagate_inbounds function albedo!(ij, albedo, vars, scheme::LandSnowAlbedo, model)

    # 1. Albedo of vegetation + bare soil (no snow)
    (; albedo_land, albedo_high_vegetation, albedo_low_vegetation) = scheme

    if haskey(vars.parameterizations.land, :high_vegetation_cover) && haskey(vars.parameterizations.land, :low_vegetation_cover)
        (; high_vegetation_cover, low_vegetation_cover) = vars.parameterizations.land

        # linear combination of high and low vegetation and bare soil
        albedo[ij] = high_vegetation_cover[ij] * albedo_high_vegetation +
            low_vegetation_cover[ij] * albedo_low_vegetation +
            albedo_land * (1 - high_vegetation_cover[ij] - low_vegetation_cover[ij])
    else
        albedo[ij] = albedo_land
    end

    # 2. Add snow cover
    if haskey(vars.prognostic, :land) && haskey(vars.prognostic.land, :snow_depth)
        ts = time_stepper(model.time_stepping, :land)
        snow_depth = get_prognostic_step(vars.prognostic.land.snow_depth, ts, scheme)
        (; albedo_snow, snow_depth_scale) = scheme

        # how to compute snow cover from snow depth
        snow_cover_scheme = scheme.snow_cover

        # compute snow-cover fraction using the chosen scheme and clamp to [0, 1]
        snow_cover = snow_cover_scheme(snow_depth[ij], snow_depth_scale)

        # set land albedo linearly between bare land and snow depending on snow cover [0, 1]
        albedo[ij] += snow_cover * albedo_snow
    end
    return nothing
end

export LearnedLandAlbedo

"""Learned albedo over land based on surface properties and some ambient conditions.
This breaks the scheme down into visible and near-infrared components and calculates
scattering parameters for each wavelength that depend on the properties of the surface. 
The parameters are combined into representative broadband shortwave parameters, then used 
to construct the white- and black-sky albedo, which combine depending on the optical 
properties of the atmosphere - to finally produce the blue-sky albedo. 
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct LearnedLandAlbedo{NF, B, Scheme <: AbstractSnowCover} <: AbstractAlbedo
    "Visible wavelength weight to the broadband shortwave albedo"
    @param vis_weight::NF = 0.5143 (bounds = 0 .. 1,)

    "Snow depth scaling for calculating snow cover fraction from snow depth"
    @param snow_depth_scale::NF = 0.05 (bounds = Positive,)

    "Snow cover-albedo scheme"
    @param snow_cover::Scheme = SaturatingSnowCover() (group = :snow_cover,)

    "[OPTION] Switch on stochastic output for the learned albedo"
    noise::B = false

    "Learned visible isotropic scattering expression"
    @component visible_isotropic::VisibleIsotropicScattering{NF} = VisibleIsotropicScattering{NF}()

    "Learned visible volumetric scattering expression"
    @component visible_volumetric::VisibleVolumetricScattering{NF} = VisibleVolumetricScattering{NF}()

    "Learned visible geometric scattering expression"
    @component visible_geometric::VisibleGeometricScattering{NF} = VisibleGeometricScattering{NF}()

    "Learned NIR isotropic scattering expression"
    @component near_infrared_isotropic::NIRIsotropicScattering{NF} = NIRIsotropicScattering{NF}()

    "Learned NIR volumetric scattering expression"
    @component near_infrared_volumetric::NIRVolumetricScattering{NF} = NIRVolumetricScattering{NF}()

    "Learned NIR geometric scattering expression"
    @component near_infrared_geometric::NIRGeometricScattering{NF} = NIRGeometricScattering{NF}()
end

Adapt.@adapt_structure LearnedLandAlbedo
function LearnedLandAlbedo(
        SG::SpectralGrid;
        snow_cover = SaturatingSnowCover(),
        kwargs...,
    )
    return LearnedLandAlbedo{SG.NF, Bool, typeof(snow_cover)}(; snow_cover, kwargs...)
end

function variables(::LearnedLandAlbedo)
    return (
        ParameterizationVariable(:black_sky_albedo, GridXY(), desc = "Black-sky land albedo", units = "1", namespace = :land),
        ParameterizationVariable(:white_sky_albedo, GridXY(), desc = "White-sky land albedo", units = "1", namespace = :land),
    )
end

initialize!(albedo::LearnedLandAlbedo, model::PrimitiveEquation) = nothing


@propagate_inbounds function land_albedo_predictors(ij, vars, scheme::LearnedLandAlbedo, model)
    # Unpack predictors for albedo parameter prediction
    (; land) = vars.parameterizations

    NF = eltype(land.albedo)
    zero_value = zero(NF)
    noise = zero_value #TODO: add noise from Gaussian noise field if scheme.noise == true

    t2m = vars.parameterizations.surface_air_temperature[ij]
    U₁₀ = vars.parameterizations.surface_wind_speed[ij]
    vₕ = haskey(land, :high_vegetation_cover) ? land.high_vegetation_cover[ij] : zero_value
    vₗ = haskey(land, :low_vegetation_cover) ? land.low_vegetation_cover[ij] : zero_value
    laiₕ = haskey(land, :high_vegetation_leaf_area) ? land.high_vegetation_leaf_area[ij] : zero_value
    laiₗ = haskey(land, :low_vegetation_leaf_area) ? land.low_vegetation_leaf_area[ij] : zero_value
    geopotential = NF(9.80665) * model.orography.orography[ij]  # trained on Earth's gravity, so use Earth's gravity

    # Get snow depth and calculate snow cover
    ts = time_stepper(model.time_stepping, :land)
    if haskey(vars.prognostic, :land) && haskey(vars.prognostic.land, :snow_depth)
        snow_depth = get_prognostic_step(vars.prognostic.land.snow_depth, ts, scheme)[ij]

        # how to compute snow cover from snow depth
        snow_cover_scheme = scheme.snow_cover

        # compute snow-cover fraction using the chosen scheme
        snow_cover = snow_cover_scheme(snow_depth, scheme.snow_depth_scale)
    else
        snow_depth = zero_value
        snow_cover = zero_value
    end

    # Get soil moisture at upper and lower layers
    if haskey(vars.prognostic.land, :soil_moisture) && haskey(land, :soil_moisture_capacity)
        soil_moisture = get_prognostic_step(
            vars.prognostic.land.soil_moisture, ts, scheme,
        )
        capacity = land.soil_moisture_capacity[ij]

        # Multiply by soil moisture capacity to get back volume of water in the layer
        sm1 = soil_moisture[ij, begin] * capacity
        sm2 = soil_moisture[ij, end] * capacity
    else
        sm1 = sm2 = zero_value
    end

    # Get soil temperature at upper and lower layers
    if haskey(vars.prognostic.land, :soil_temperature)
        soil_temperature = get_prognostic_step(
            vars.prognostic.land.soil_temperature, ts, scheme,
        )
        st1 = soil_temperature[ij, begin]
        st2 = soil_temperature[ij, end]
    else
        st1 = st2 = zero_value
    end

    # Change in soil temperature helps the model to detect seasonal shift
    ΔTₛ = st1 - st2

    # from precipitation schemes [m/s]
    snow_fall_rate = haskey(vars.parameterizations, :snow_rate) ? vars.parameterizations.snow_rate[ij] : zero_value

    return (; vₕ, vₗ, laiₕ, laiₗ, sm1, sm2, st1, st2, ΔTₛ, geopotential, snow_cover, t2m, U₁₀, snow_fall_rate, noise, snow_depth)
end

@propagate_inbounds function broadband_parameter(vis_weight, vis_param, nir_param)
    # Blend parameters based on empirical VIS/NIR weights (Munday et al., in prep)
    nir_weight = 1 - vis_weight
    broadband_param = (vis_param * vis_weight) + (nir_param * nir_weight)

    return broadband_param
end

@propagate_inbounds function white_sky_albedo(f_iso, f_vol, f_geo)
    # Calculates the SZA-independent White-Sky Albedo (WSA).
    # Based on the Lucht et al. (2000) BRDF model approximation.
    NF = typeof(f_iso)

    # Coefficients taken from Lucht et al. (2000)
    wsa_vol_const = NF(0.189184)
    wsa_geo_const = NF(-1.377622)
    return f_iso + (wsa_vol_const * f_vol) + (wsa_geo_const * f_geo)
end

@propagate_inbounds function black_sky_albedo(f_iso, f_vol, f_geo, θ)
    # Calculate black sky albedo using the Lucht et al. (2000) polynomial
    NF = typeof(f_iso)

    # Coefficients taken from Lucht et al. (2000)
    c_iso = NF(1.0)
    c_vol = evalpoly(θ, (NF(-0.007574), NF(0.0), NF(-0.070987), NF(0.307588)))
    c_geo = evalpoly(θ, (NF(-1.284909), NF(0.0), NF(-0.166314), NF(0.04184)))

    return (c_iso * f_iso) + (c_vol * f_vol) + (c_geo * f_geo)
end

@propagate_inbounds function albedo!(ij, albedo, vars, scheme::LearnedLandAlbedo, model)
    (; land) = vars.parameterizations
    albedo_predictors = land_albedo_predictors(ij, vars, scheme, model)

    isotropic_parameter = broadband_parameter(scheme.vis_weight, scheme.visible_isotropic(albedo_predictors), scheme.near_infrared_isotropic(albedo_predictors))
    volumetric_parameter = broadband_parameter(scheme.vis_weight, scheme.visible_volumetric(albedo_predictors), scheme.near_infrared_volumetric(albedo_predictors))
    geometric_parameter = broadband_parameter(scheme.vis_weight, scheme.visible_geometric(albedo_predictors), scheme.near_infrared_geometric(albedo_predictors))

    # Get zenith angle in radians
    μ = vars.parameterizations.cos_zenith[ij]
    μ = clamp(μ, zero(μ), one(μ))  # only need when the sun is up
    θ = atan(sqrt(1 - μ^2), μ)  # TODO: Workaround for acos, to be reverted back when Reactant raised correctly to stablehlo.acos

    # Then from those we'll calculate the black and white sky albedos
    black_sky = black_sky_albedo(isotropic_parameter, volumetric_parameter, geometric_parameter, θ)
    white_sky = white_sky_albedo(isotropic_parameter, volumetric_parameter, geometric_parameter)
    land.black_sky_albedo[ij] = black_sky
    land.white_sky_albedo[ij] = white_sky

    # Fallback for radiation schemes with no direct/diffuse radiation fraction
    albedo[ij] = white_sky

    return nothing
end