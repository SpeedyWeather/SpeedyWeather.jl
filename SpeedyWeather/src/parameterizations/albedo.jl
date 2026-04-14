"""
Albedo structs should be defined as MyAlbedo <: AbstractAlbedo and implement the following interface:

- `initialize!(albedo::MyAlbedo, model)`
- `albedo!(ij, vars::Variables, albedo::MyAlbedo, model)`

`initialize!` is as for every model component but in contrast to other parameterizations albedos
should not implement `parameterization!` but rather `albedo!`, which is called from within `parameterization!`.
This is because every albedo needs to be used either as part of OceanLandAlbedo or the same albedo
is called over both land and ocean surfaces."""
abstract type AbstractAlbedo <: AbstractModelComponent end

# all albedos need to define albedo fields for ocean and land separately and total
function variables(::AbstractAlbedo)
    return (
        ParameterizationVariable(:albedo, Grid2D(), desc = "Albedo", units = "1"),
        ParameterizationVariable(:albedo, Grid2D(), desc = "Albedo over the ocean", units = "1", namespace = :ocean),
        ParameterizationVariable(:albedo, Grid2D(), desc = "Albedo over the land", units = "1", namespace = :land),
    )
end

export OceanLandAlbedo

"""Composite type: Two albedo parameterizations, one for ocean and one for land surfaces.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct OceanLandAlbedo{Ocean, Land} <: AbstractAlbedo
    "[OPTION] Albedo parameterization for ocean surfaces"
    @component ocean::Ocean

    "[OPTION] Albedo parameterization for land surfaces"
    @component land::Land
end

Adapt.@adapt_structure OceanLandAlbedo
function OceanLandAlbedo(
        SG::SpectralGrid;
        ocean = OceanSeaIceAlbedo(SG),      # default ocean albedo
        land = LandSnowAlbedo(SG),          # default land albedo
    )
    return OceanLandAlbedo(ocean, land)
end

function Base.show(io::IO, A::OceanLandAlbedo)
    println(io, "OceanLandAlbedo <: SpeedyWeather.AbstractAlbedo")
    properties = propertynames(A)
    n = length(properties)
    for (i, key) in enumerate(properties)
        val = getfield(A, key)
        s = i == n ? "└" : "├"  # choose ending └ for last property
        p = i == n ? print : println
        p(io, "$s $key: $(typeof(val))")
    end
    return
end

# TODO deprecate?
export DefaultAlbedo

"""$(TYPEDSIGNATURES) Default albedo parameterization."""
DefaultAlbedo(SG::SpectralGrid; kwargs...) = OceanLandAlbedo(SG; kwargs...)

# composite albedos need to define albedo for ocean/land separately and combined
# plus any custom variables of the ocean and land albedo parameterizations (if specific method is implemented)
function variables(albedo::OceanLandAlbedo)
    return (
        ParameterizationVariable(:albedo, Grid2D(), desc = "Albedo", units = "1"),
        ParameterizationVariable(:albedo, Grid2D(), desc = "Albedo over the ocean", units = "1", namespace = :ocean),
        ParameterizationVariable(:albedo, Grid2D(), desc = "Albedo over the land", units = "1", namespace = :land),
        variables(albedo.ocean)...,
        variables(albedo.land)...,
    )
end

function initialize!(albedo::OceanLandAlbedo, model::PrimitiveEquation)
    initialize!(albedo.ocean, model)
    initialize!(albedo.land, model)
    return nothing
end

"""$(TYPEDSIGNATURES) OceanLandAlbedo: Call ocean and land albedo parameterizations separately."""
@propagate_inbounds function parameterization!(ij, vars::Variables, albedo::OceanLandAlbedo, model)
    albedo!(ij, vars, albedo.ocean, model)
    return albedo!(ij, vars, albedo.land, model)
end

"""$(TYPEDSIGNATURES) Single Albedo: Call same albedo over ocean and land."""
@propagate_inbounds function parameterization!(ij, vars::Variables, albedo::AbstractAlbedo, model)
    albedo!(ij, vars, albedo.ocean, model)
    return albedo!(ij, vars, albedo.land, model)
end

export GlobalConstantAlbedo

"""Global constant albedo parameterization. To be used for land and ocean 
or only one of them within a `OceanLandAlbedo`. Fields are $(TYPEDFIELDS)"""
@kwdef struct GlobalConstantAlbedo{NF} <: AbstractAlbedo
    "[OPTION] Albedo value [1]"
    albedo::NF = 0.3
end

Adapt.@adapt_structure GlobalConstantAlbedo
GlobalConstantAlbedo(SG::SpectralGrid; kwargs...) = GlobalConstantAlbedo{SG.NF}(; kwargs...)
initialize!(albedo::GlobalConstantAlbedo, ::PrimitiveEquation) = nothing
@propagate_inbounds albedo!(ij, vars, scheme::GlobalConstantAlbedo, model) = (vars.parameterizations.albedo[ij] = scheme.albedo)

export ManualAlbedo

"""Manual albedo field, to be used with `set!` and is copied into the diagnostic variables on every time step.
Defined so that parameterizations can change the albedo at every time step (e.g. snow cover) without
losing the information of the original surface albedo. Fields are
$(TYPEDFIELDS)"""
struct ManualAlbedo{GridVariable2D} <: AbstractAlbedo
    "Albedo field [1]"
    albedo::GridVariable2D
end

Adapt.@adapt_structure ManualAlbedo
ManualAlbedo(SG::SpectralGrid) = ManualAlbedo{SG.GridVariable2D}(zeros(SG.GridVariable2D, SG.grid))
initialize!(albedo::ManualAlbedo, model::PrimitiveEquation) = nothing
@propagate_inbounds albedo!(ij, vars, scheme::ManualAlbedo, model) = (vars.parameterizations.albedo[ij] = scheme.albedo[ij])

export AlbedoClimatology

"""Albedo climatology loaded from netcdf file. Fields are $(TYPEDFIELDS)"""
@kwdef struct AlbedoClimatology{GridVariable2D} <: AbstractAlbedo
    "[OPTION] filename of albedo"
    file::String = "albedo.nc"

    "[OPTION] path to the folder containing the soil moisture"
    path::String = joinpath("data", "boundary_conditions", file)

    "[OPTION] flag to check for soil moisture in SpeedyWeatherAssets or locally"
    from_assets::Bool = true

    "[OPTION] variable name in netcdf file"
    varname::String = "alb"

    "[OPTION] SpeedyWeatherAssets version number"
    version::VersionNumber = DEFAULT_ASSETS_VERSION

    "[OPTION] Grid the albedo file comes on"
    FieldType::Type{<:AbstractField} = FullGaussianField

    "Albedo climatology"
    albedo::GridVariable2D
end

# For GPU usage just discard the extra information and treat it as a `ManualAlbedo`
Adapt.adapt_structure(to, albedo::AlbedoClimatology) = adapt(to, ManualAlbedo(albedo.albedo))

function AlbedoClimatology(SG::SpectralGrid; kwargs...)
    (; GridVariable2D, grid) = SG
    albedo = zeros(GridVariable2D, grid)
    return AlbedoClimatology{GridVariable2D}(; albedo, kwargs...)
end

# set albedo with grid, scalar, function; just define path `albedo.albedo` to grid here
set!(albedo::AbstractAlbedo, args...; kwargs...) = set!(albedo.albedo, args...; kwargs...)

function initialize!(albedo::AlbedoClimatology, model::PrimitiveEquation)

    # LOAD NETCDF FILE
    a = get_asset(
        albedo.path;
        from_assets = albedo.from_assets,
        name = albedo.varname,
        ArrayType = albedo.FieldType,
        FileFormat = NCDataset,
        version = albedo.version
    )

    return interpolate!(albedo.albedo, a)
end

@propagate_inbounds albedo!(ij, vars, scheme::AlbedoClimatology, model) = (vars.parameterizations.albedo[ij] = scheme.albedo[ij])

# OceanSeaIceAlbedo
export OceanSeaIceAlbedo

"""Albedo that scales linearly between ocean and ice albedo depending on sea ice concentration.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct OceanSeaIceAlbedo{NF} <: AbstractAlbedo
    "[OPTION] Albedo over open ocean [1]"
    @param albedo_ocean::NF = 0.06 (bounds = 0 .. 1,)

    "[OPTION] Albedo over sea ice at concentration=1 [1]"
    @param albedo_ice::NF = 0.6 (bounds = 0 .. 1,)
end

Adapt.@adapt_structure OceanSeaIceAlbedo
OceanSeaIceAlbedo(SG::SpectralGrid; kwargs...) = OceanSeaIceAlbedo{SG.NF}(; kwargs...)
initialize!(::OceanSeaIceAlbedo, ::PrimitiveEquation) = nothing

@propagate_inbounds function albedo!(ij, vars, scheme::OceanSeaIceAlbedo, model)
    (; albedo_ocean, albedo_ice) = scheme
    NF = eltype(vars.parameterizations.ocean.albedo)
    ℵ = haskey(vars.prognostic.ocean, :sea_ice_concentration) ? vars.prognostic.ocean.sea_ice_concentration[ij] : zero(NF)

    # set ocean albedo linearly between ocean and ice depending on sea ice concentration
    return albedo[ij] = albedo_ocean + ℵ * (albedo_ice - albedo_ocean)
end

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

@propagate_inbounds function albedo!(ij, vars, scheme::LandSnowAlbedo, model)

    # 1. Albedo of vegetation + bare soil (no snow)
    (; albedo_land, albedo_high_vegetation, albedo_low_vegetation) = scheme
    (; albedo) = vars.parameterizations.land

    if haskey(vars.parameterizations.land, :vegetation_high) && haskey(vars.parameterizations.land, :vegetation_low)
        (; vegetation_high, vegetation_low) = vars.parameterizations.land

        # linear combination of high and low vegetation and bare soil
        albedo[ij] = vegetation_high[ij] * albedo_high_vegetation +
            vegetation_low[ij] * albedo_low_vegetation +
            albedo_land * (1 - vegetation_high[ij] - vegetation_low[ij])
    else
        albedo[ij] = albedo_land
    end

    # 2. Add snow cover
    if haskey(vars.prognostic.land, :snow_depth)
        (; snow_depth) = vars.prognostic.land
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

## JinOceanAlbedo (Jin et al., 2011)
export JinOceanAlbedo

"""
Ocean surface albedo parameterization based on Jin et al. (2011).
Calculates broadband albedo as a function of wind speed, surface roughness, 
and solar zenith angle, partitioned into direct and diffuse components.
Fields are $(TYPEDFIELDS)
"""
@parameterized @kwdef struct JinOceanAlbedo{NF} <: AbstractAlbedo
    "[OPTION] Refractive index of water for broadband visible spectrum [1]"
    @param refractive_index::NF = 1.34 (bounds = 1 .. 2,)

    "[OPTION] Constant foam reflectance (Koepke 1984) [1]"
    @param foam_reflectance::NF = 0.55 (bounds = 0 .. 1,)

    "[OPTION] Volume scattering contribution for Case 1 waters [1]"
    @param volume_scattering::NF = 0.006 (bounds = 0 .. 0.1,)

    "[OPTION] Albedo of sea ice at concentration=1 [1]"
    @param albedo_ice::NF = 0.6 (bounds = 0 .. 1,)

    "[OPTION] Use overcast skies formula for diffuse albedo?"
    cloudy_sky::Bool = false
end

Adapt.@adapt_structure JinOceanAlbedo
JinOceanAlbedo(SG::SpectralGrid; kwargs...) = JinOceanAlbedo{SG.NF}(; kwargs...)

initialize!(::JinOceanAlbedo, ::PrimitiveEquation) = nothing

"""$(TYPEDSIGNATURES) Compute unpolarized Fresnel reflectance."""
@inline function fresnel_reflectance(n::NF, mu::NF) where {NF}
    # Mu is bounded to prevent domain errors
    μ = clamp(mu, NF(1.0e-5), one(NF))
    θᵢ = acos(μ)

    sin_theta_t = sin(θᵢ) / n
    θₜ = asin(sin_theta_t)

    # Perpendicular and parallel polarizations
    rₛ = ((cos(θᵢ) - n * cos(θₜ)) / (cos(θᵢ) + n * cos(θₜ)))^2
    rₚ = ((n * cos(θᵢ) - cos(θₜ)) / (n * cos(θᵢ) + cos(θₜ)))^2

    return NF(0.5) * (rₛ + rₚ)
end

"""$(TYPEDSIGNATURES) Polynomial regression function f(mu, sigma) from Jin et al. (2011)."""
@inline function f_mu_sigma(mu::NF, sigma::NF) where {NF}
    p1, p2, p3 = NF(0.0152), NF(-1.7873), NF(6.8972)
    p4, p5, p6 = NF(-8.5778), NF(4.071), NF(-7.6446)
    p7, p8, p9 = NF(0.1643), NF(-7.8409), NF(-3.5639)
    p10, p11 = NF(-2.3588), NF(10.0538)

    poly_pre = p1 + p2 * mu + p3 * mu^2 + p4 * mu^3 + p5 * sigma + p6 * mu * sigma
    poly_exp = exp(p7 + p8 * mu + p9 * mu^2 + p10 * sigma + p11 * mu * sigma)

    return poly_pre * poly_exp
end

@propagate_inbounds function albedo!(ij, vars, albedo_scheme::JinOceanAlbedo, model)
    NF = model.spectral_grid.NF
    land_fraction = model.land_sea_mask.mask[ij]
    μ = vars.parameterizations.cos_zenith[ij]

    # Exit if the point is land
    if land_fraction == 1
        vars.parameterizations.ocean.albedo[ij] = zero(land_fraction)
        return
    end

    # Enforce night/horizon boundary
    if μ <= zero(NF)
        return
    end

    (; refractive_index, foam_reflectance, volume_scattering, cloudy_sky, albedo_ice) = albedo_scheme

    wind_speed = vars.parameterizations.surface_wind_speed[ij]
    f_dir = vars.parameterizations.direct_radiation_fraction[ij]
    f_dif = vars.parameterizations.diffuse_radiation_fraction[ij]
    σ = vars.parameterizations.ocean.surface_roughness[ij]

    n₀ = refractive_index

    # 1. Direct Surface Albedo
    r_f = fresnel_reflectance(n₀, μ)
    alpha_dir_s = r_f - f_mu_sigma(μ, σ)

    # 2. Diffuse Surface Albedo
    if cloudy_sky
        # Eq 5b (Cloudy sky)
        alpha_dif_s = NF(-0.1479) + NF(0.1502) * n₀ - NF(0.0176) * n₀ * σ
    else
        # Eq 5a (Isotropic clear sky)
        alpha_dif_s = NF(-0.1482) - NF(0.012) * n₀ + NF(0.1608) * n₀^2 - NF(0.0244) * n₀ * σ
    end

    # 3. Assemble Broadband Albedo (Eq 15)
    αᵦ = f_dir * alpha_dir_s + f_dif * alpha_dif_s + volume_scattering

    # 4. Foam Adjustment (Eqs 16 & 17)
    f_wc = NF(2.95e-6) * (wind_speed^NF(3.52))
    f_wc = clamp(f_wc, zero(NF), one(NF))

    alpha_final = f_wc * foam_reflectance + (one(NF) - f_wc) * αᵦ

    # Also account for sea ice if sea ice concentration is available
    if haskey(vars.prognostic.ocean, :sea_ice_concentration)
        sea_ice_conc = vars.prognostic.ocean.sea_ice_concentration[ij]
        vars.parameterizations.ocean.albedo[ij] = clamp((one(NF) - sea_ice_conc) * alpha_final + sea_ice_conc * albedo_ice, zero(NF), one(NF))
    else
        vars.parameterizations.ocean.albedo[ij] = clamp(alpha_final, zero(NF), one(NF))
    end

    return
end

export LearnedLandAlbedo
@parameterized @kwdef struct LearnedLandAlbedo{NF, LNN, LP, LS, Scheme <: AbstractSnowCover} <: AbstractAlbedo
    "[OPTION] filename of land weights"
    file::String = "brdf.npz"

    "[OPTION] path to the folder containing the neural network weights"
    path::String = joinpath("data", "weights", file)

    "[OPTION] flag to check for weights in SWA or locally"
    from_assets::Bool = true

    "[OPTION] SpeedyWeatherAssets version number"
    version::VersionNumber = DEFAULT_ASSETS_VERSION

    "Conversion from snow depth to snow cover [m]"
    @param snow_depth_scale::NF = 0.05 (bounds = Positive,)

    "Snow cover-albedo scheme"
    @param snow_cover::Scheme = SaturatingSnowCover() (group = :snow_cover,)

    input_buffer::Matrix{NF}

    # Normalisation parameters
    norm_means::Vector{Float32} = zeros(Float32, 10)
    norm_stds::Vector{Float32} = zeros(Float32, 10)
    unnorm_means::Vector{Float32} = zeros(Float32, 6)
    unnorm_stds::Vector{Float32} = zeros(Float32, 6)

    brdf_nn::LNN
    brdf_params::LP
    brdf_states::LS
end

function Base.show(io::IO, scheme::LearnedLandAlbedo)
    print(io, "LearnedLandAlbedo{$(eltype(scheme.input_buffer))}")
    println(io)
    n_layers = length(keys(scheme.brdf_params))
    return println(io, "└ BRDF: Neural network ($n_layers layers)")
end
