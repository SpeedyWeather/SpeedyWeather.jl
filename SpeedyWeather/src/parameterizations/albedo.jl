"""
Albedo structs should be defined as MyAlbedo <: AbstractAlbedo and implement the following interface:

- `initialize!(albedo::MyAlbedo, model)`
- `albedo!(ij, albedofield, vars::Variables, albedo::MyAlbedo, model)`

`initialize!` is as for every model component but in contrast to other parameterizations albedos
should not implement `parameterization!` but rather `albedo!`, which is called from within `parameterization!`.
This is because every albedo needs to be used either as part of OceanLandAlbedo or the same albedo
is called over both land and ocean surfaces. If a surface namespace (`vars.parameterizations.land` or `.ocean`) 
contains both `black_sky_albedo` and `white_sky_albedo`, a shortwave scheme that knows the direct-beam fraction 
blends them into that namespace's albedo. Otherwise it leaves albedo alone."""
abstract type AbstractAlbedo <: AbstractParameterization end

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
        ocean = JinOceanAlbedo(SG),      # default ocean albedo
        land = LearnedLandAlbedo(SG),       # default land albedo
    )
    return OceanLandAlbedo(ocean, land)
end

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

# abbreviate show for OceanLandAlbedo to avoid printing all fields of ocean and land albedo in the main show
Base.show(io::IO, A::OceanLandAlbedo) = show(io, A, values = false)

function initialize!(albedo::OceanLandAlbedo, model::PrimitiveEquation)
    initialize!(albedo.ocean, model)
    initialize!(albedo.land, model)
    return nothing
end

"""$(TYPEDSIGNATURES) Composite OceanLandAlbedo: Call albedo.ocean over ocean and albedo.land over land."""
@propagate_inbounds function parameterization!(ij, vars::Variables, albedo::OceanLandAlbedo, model)
    albedo!(ij, vars.parameterizations.ocean.albedo, vars, albedo.ocean, model)
    return albedo!(ij, vars.parameterizations.land.albedo, vars, albedo.land, model)
end

"""$(TYPEDSIGNATURES) Single Albedo: Call same albedo over ocean and land."""
@propagate_inbounds function parameterization!(ij, vars::Variables, albedo::AbstractAlbedo, model)
    albedo!(ij, vars.parameterizations.ocean.albedo, vars, albedo, model)
    return albedo!(ij, vars.parameterizations.land.albedo, vars, albedo, model)
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
@propagate_inbounds albedo!(ij, albedo, vars, scheme::GlobalConstantAlbedo, model) = (albedo[ij] = scheme.albedo)

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
@propagate_inbounds albedo!(ij, albedo, vars, scheme::ManualAlbedo, model) = (albedo[ij] = scheme.albedo[ij])

export PrescribedAlbedo

"""Prescribed albedo field, to be used with `set!`. In contrast to `ManualAlbedo` the albedo is not 
hold in the struct and not copied into the diagnostic variables automatically. It's expected that the
albedo variable (e.g `vars.land.albedo`) is set externally at every time step.
Defined so that the albedo can be prescribed externally e.g. when coupling the model. """
struct PrescribedAlbedo <: AbstractAlbedo
end

Adapt.@adapt_structure PrescribedAlbedo
PrescribedAlbedo(SG::SpectralGrid) = PrescribedAlbedo()
initialize!(albedo::PrescribedAlbedo, model::PrimitiveEquation) = nothing
@propagate_inbounds albedo!(ij, albedo, vars, scheme::PrescribedAlbedo, model) = nothing

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

    a = on_architecture(model.architecture, a)
    return interpolate!(albedo.albedo, a)
end

@propagate_inbounds albedo!(ij, albedo, vars, scheme::AlbedoClimatology, model) = (albedo[ij] = scheme.albedo[ij])

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

@propagate_inbounds function albedo!(ij, albedo, vars, scheme::OceanSeaIceAlbedo, model)
    (; albedo_ocean, albedo_ice) = scheme
    NF = eltype(albedo)
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

abstract type AbstractLearnedScattering end

"""Visible radiation isotropic scattering parameter calculation. This expression was derived by
distilling a trained neural network using symbolic regression. Isotropic scattering means that
light at the surface is reflected equally in all directions, regardless of where the sun is.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct VisibleIsotropicScattering{NF} <: AbstractLearnedScattering
    "Coefficient fit via symbolic regression"
    @param c1::NF = 0.13926611446264833

    "Coefficient fit via symbolic regression"
    @param c2::NF = 0.049407226600965096

    "Coefficient fit via symbolic regression"
    @param c3::NF = 0.03012000138970739

    "Coefficient fit via symbolic regression"
    @param c4::NF = 0.007003960811839407

    "Coefficient fit via symbolic regression"
    @param c5::NF = 0.0030819533608978183

    "Coefficient fit via symbolic regression"
    @param c6::NF = 98.73431990602951

    "Coefficient fit via symbolic regression"
    @param c7::NF = 0.0007767787214513013

    "Coefficient fit via symbolic regression"
    @param c8::NF = 0.4746273568956652

    "Coefficient fit via symbolic regression"
    @param c9::NF = 0.004805110395235223

    "Coefficient fit via symbolic regression"
    @param c10::NF = 1.1381888960412316

    "Coefficient fit via symbolic regression"
    @param c11::NF = 0.0203040824498551

end

Adapt.@adapt_structure VisibleIsotropicScattering
VisibleIsotropicScattering(SG::SpectralGrid; kwargs...) = VisibleIsotropicScattering{SG.NF}(; kwargs...)

@inline function (c::VisibleIsotropicScattering)(p)
    (; laiₗ, sm1, sm2, snow_cover, st2, t2m, vₕ, vₗ, ΔTₛ, noise) = p
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11) = c

    # Mean calculation
    μ = (
        -(c1 * vₕ) + vₗ * (-(c2 * laiₗ) - c2 * snow_cover - c2 * sm1 - c3) + snow_cover * (c4 - c5 * vₕ)
            * (ΔTₛ + c6) - c1 * sm2 - c7 * t2m + c8
    )

    # Std deviation
    raw_σ = snow_cover * sm2 * (c9 * st2 - c10) + c11
    σ = max(raw_σ, zero(raw_σ))  # shouldn't be negative...
    return μ + σ * noise
end

"""Visible radiation volumetric scattering parameter calculation. This expression was derived by
distilling a trained neural network using symbolic regression. Volumetric scattering represents the
reflection mechanism of light interacting with horizontally homogenous vegetation canopy at the surface - 
light meets the canopy, bounces around a bit beneath the canopy and exits back towards the sky. This 
does depend on where the sun is and the viewing angle.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct VisibleVolumetricScattering{NF} <: AbstractLearnedScattering
    "Coefficient fit via symbolic regression"
    @param c1::NF = 0.005137089001311701

    "Coefficient fit via symbolic regression"
    @param c2::NF = 8.359579016151135e-7

    "Coefficient fit via symbolic regression"
    @param c3::NF = 547.1867710914254

    "Coefficient fit via symbolic regression"
    @param c4::NF = 4.88214507804199

    "Coefficient fit via symbolic regression"
    @param c5::NF = 0.006321233141025383

    "Coefficient fit via symbolic regression"
    @param c6::NF = 547.1858817550591

    "Coefficient fit via symbolic regression"
    @param c7::NF = 7.789539669853218e-8

    "Coefficient fit via symbolic regression"
    @param c8::NF = 0.0012984891007226625

    "Coefficient fit via symbolic regression"
    @param c9::NF = 0.013366006624365259

    "Coefficient fit via symbolic regression"
    @param c10::NF = 0.037107779858204085

    "Coefficient fit via symbolic regression"
    @param c11::NF = 0.001076347986980297

    "Coefficient fit via symbolic regression"
    @param c12::NF = 0.02154070145308177

    "Coefficient fit via symbolic regression"
    @param c13::NF = 0.013735265859374146

end

Adapt.@adapt_structure VisibleVolumetricScattering
VisibleVolumetricScattering(SG::SpectralGrid; kwargs...) = VisibleVolumetricScattering{SG.NF}(; kwargs...)

@inline function (c::VisibleVolumetricScattering)(p)
    (; geopotential, laiₗ, sm1, snow_cover, st1, t2m, vₕ, vₗ, ΔTₛ, noise) = p
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13) = c
    # Vegetation cover term
    q0 = vₕ + vₗ

    # Mean calculation
    μ = (
        q0 * (
            -(c1 * laiₗ) + c2 * t2m * (
                ΔTₛ + snow_cover * (
                    -(c3 * snow_cover) + c4 * t2m + c5 *
                        geopotential - c6
                ) + c7 * geopotential^2
            ) / (q0 + c8) - c9
        ) + c8 * st1 - c10 * sm1 - c11 * t2m
    )

    # Std devation calculation
    raw_σ = c12 * vₗ * snow_cover + c13
    σ = max(raw_σ, zero(raw_σ))  # can't be negative
    return μ + σ * noise
end

"""Visible radiation geometric-optical scattering parameter calculation. This expression was derived by
distilling a trained neural network using symbolic regression. Geometric scattering represents the 
reflection due to shadowcasting of 3D objects at the surface, like tree crowns or big rocks, and 
obscuration of shadows due to viewing angle.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct VisibleGeometricScattering{NF} <: AbstractLearnedScattering
    "Coefficient fit via symbolic regression"
    @param c1::NF = 98.30689271056275

    "Coefficient fit via symbolic regression"
    @param c2::NF = 0.0002465039456962893

    "Coefficient fit via symbolic regression"
    @param c3::NF = 1.74846262394399e-8

    "Coefficient fit via symbolic regression"
    @param c4::NF = 0.2655102008045942

    "Coefficient fit via symbolic regression"
    @param c5::NF = 0.00016722961466941218

    "Coefficient fit via symbolic regression"
    @param c6::NF = 1.8198639960846736e-6

    "Coefficient fit via symbolic regression"
    @param c7::NF = 0.009976526656995078

    "Coefficient fit via symbolic regression"
    @param c8::NF = 0.000145787380943221

    "Coefficient fit via symbolic regression"
    @param c9::NF = 0.012478157859566873

    "Coefficient fit via symbolic regression"
    @param c10::NF = 0.01037571803805248

    "Coefficient fit via symbolic regression"
    @param c11::NF = 0.007018418642708496

    "Coefficient fit via symbolic regression"
    @param c12::NF = 0.000694635309383557

    "Coefficient fit via symbolic regression"
    @param c13::NF = 0.005474843497829615

end

Adapt.@adapt_structure VisibleGeometricScattering
VisibleGeometricScattering(SG::SpectralGrid; kwargs...) = VisibleGeometricScattering{SG.NF}(; kwargs...)

@inline function (c::VisibleGeometricScattering)(p)
    (; U₁₀, geopotential, laiₗ, sm1, snow_cover, snow_depth, vₕ, vₗ, ΔTₛ, noise) = p
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13) = c

    # Mean calculation
    μ = (
        snow_cover * (snow_depth * (ΔTₛ - U₁₀) + c1) * (
            -(c2 * vₗ) + c3 * geopotential * (
                vₕ + laiₗ * (vₗ + snow_depth) - c4
            ) + c5
        ) + (
            c6 * geopotential + (-(c7 * vₕ) - c7 * vₗ) * (
                laiₗ * sm1 + c8 * geopotential
            ) + c9
        ) * sqrt(vₗ) + c10
    )

    # Std deviation
    raw_σ = (
        snow_cover * (c11 - c12 * snow_depth) + c13
    )
    σ = max(raw_σ, zero(raw_σ))  # can't be negative
    return μ + σ * noise
end

"""Near-infrared radiation isotropic scattering parameter calculation. This expression was derived by
distilling a trained neural network using symbolic regression. Isotropic scattering means that
light at the surface is reflected equally in all directions, regardless of where the sun is.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct NIRIsotropicScattering{NF} <: AbstractLearnedScattering
    "Coefficient fit via symbolic regression"
    @param c1::NF = 0.13250859400402135

    "Coefficient fit via symbolic regression"
    @param c2::NF = 0.08409939463271902

    "Coefficient fit via symbolic regression"
    @param c3::NF = 0.03616646161624834

    "Coefficient fit via symbolic regression"
    @param c4::NF = 0.10051756565160144

    "Coefficient fit via symbolic regression"
    @param c5::NF = 0.200655651610176

    "Coefficient fit via symbolic regression"
    @param c6::NF = 0.13125739630648842

    "Coefficient fit via symbolic regression"
    @param c7::NF = 1.2330577987093237e-5

    "Coefficient fit via symbolic regression"
    @param c8::NF = 0.308490327357136

    "Coefficient fit via symbolic regression"
    @param c9::NF = 0.38361527926198996

    "Coefficient fit via symbolic regression"
    @param c10::NF = 0.5104495709936557

    "Coefficient fit via symbolic regression"
    @param c11::NF = 0.07677930381946513

    "Coefficient fit via symbolic regression"
    @param c12::NF = 0.052701989554117486

    "Coefficient fit via symbolic regression"
    @param c13::NF = 0.25930825098307064

    "Coefficient fit via symbolic regression"
    @param c14::NF = 0.047761528867464

    "Coefficient fit via symbolic regression"
    @param c15::NF = 0.29940528501431796

    "Coefficient fit via symbolic regression"
    @param c16::NF = 0.03240815843538468

end

Adapt.@adapt_structure NIRIsotropicScattering
NIRIsotropicScattering(SG::SpectralGrid; kwargs...) = NIRIsotropicScattering{SG.NF}(; kwargs...)

@inline function (c::NIRIsotropicScattering)(p)
    (; geopotential, laiₗ, sm2, snow_cover, vₕ, vₗ, noise) = p
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16) = c

    # Mean calculation
    μ = (
        vₕ * (-(c1 * snow_cover) - c2) + laiₗ * (-(c3 * vₗ) + c4 * sm2) + c5 * snow_cover^3 + (
            - (c6 * laiₗ) + c7 * geopotential * (c8 * snow_cover + sm2 - c9)
        ) / (laiₗ + snow_cover + c10) - (c11 * sm2 - c12) / sqrt(c13 * sm2 + c14) + c15
    )

    # Std deviation, in this case a constant
    σ = c16
    return μ + σ * noise
end

"""Near-infrared radiation volumetric scattering parameter calculation. This expression was derived by
distilling a trained neural network using symbolic regression. Volumetric scattering represents the
reflection mechanism of light interacting with horizontally homogenous vegetation canopy at the surface - 
light meets the canopy, bounces around a bit beneath the canopy and exits back towards the sky. This 
does depend on where the sun is and the viewing angle.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct NIRVolumetricScattering{NF} <: AbstractLearnedScattering
    "Coefficient fit via symbolic regression"
    @param c1::NF = 0.02644856368910439

    "Coefficient fit via symbolic regression"
    @param c2::NF = 2.7037305935790665

    "Coefficient fit via symbolic regression"
    @param c3::NF = 2.2973971005827364e-5

    "Coefficient fit via symbolic regression"
    @param c4::NF = 0.0063074532998523255

    "Coefficient fit via symbolic regression"
    @param c5::NF = 0.015870141431257815

    "Coefficient fit via symbolic regression"
    @param c6::NF = 7.900349778266809e-8

    "Coefficient fit via symbolic regression"
    @param c7::NF = 1.3580861013567708e-11

    "Coefficient fit via symbolic regression"
    @param c8::NF = 0.1966651392522715

    "Coefficient fit via symbolic regression"
    @param c9::NF = 0.026503362065879634

end

Adapt.@adapt_structure NIRVolumetricScattering
NIRVolumetricScattering(SG::SpectralGrid; kwargs...) = NIRVolumetricScattering{SG.NF}(; kwargs...)

@inline function (c::NIRVolumetricScattering)(p)
    (; sm2, snow_cover, st1, st2, t2m, vₗ, noise) = p
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9) = c

    # Mean calculation
    μ = (
        snow_cover * (c1 * vₗ + c2) + sm2 * (c3 * st1^2 - c4 * st2) + t2m * (
            - (c5 * snow_cover) + t2m^2 * (c6 * snow_cover - c7 * st1)
        ) + c8
    )

    # Std deviation, in this case just a constant
    σ = c9
    return μ + σ * noise
end

"""Near-infrared radiation geometric-optical scattering parameter calculation. This expression was derived by
distilling a trained neural network using symbolic regression. Geometric scattering represents the 
reflection due to shadowcasting of 3D objects at the surface, like tree crowns or big rocks, and 
obscuration of shadows due to viewing angle.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct NIRGeometricScattering{NF} <: AbstractLearnedScattering
    "Coefficient fit via symbolic regression"
    @param c1::NF = 0.00016970668634889394

    "Coefficient fit via symbolic regression"
    @param c2::NF = 4.1657057100227406e-10

    "Coefficient fit via symbolic regression"
    @param c3::NF = 1.8792166707362004e-7

    "Coefficient fit via symbolic regression"
    @param c4::NF = 2.0653707823044138e-5

    "Coefficient fit via symbolic regression"
    @param c5::NF = 0.04338449131472477

    "Coefficient fit via symbolic regression"
    @param c6::NF = 0.10172346530005257

    "Coefficient fit via symbolic regression"
    @param c7::NF = 0.07360077673775481

    "Coefficient fit via symbolic regression"
    @param c8::NF = 3.49400847746966e-5

    "Coefficient fit via symbolic regression"
    @param c9::NF = 0.005699938844269787

    "Coefficient fit via symbolic regression"
    @param c10::NF = 0.0008285983510980352

    "Coefficient fit via symbolic regression"
    @param c11::NF = 0.614697095826219

    "Coefficient fit via symbolic regression"
    @param c12::NF = 0.0790756374361683

    "Coefficient fit via symbolic regression"
    @param c13::NF = 0.009887490927614274

end

Adapt.@adapt_structure NIRGeometricScattering
NIRGeometricScattering(SG::SpectralGrid; kwargs...) = NIRGeometricScattering{SG.NF}(; kwargs...)

@inline function (c::NIRGeometricScattering)(p)
    (; geopotential, laiₕ, sm2, snow_cover, st2, t2m, vₕ, vₗ, noise) = p
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13) = c

    # Mean calculation
    μ = (
        -(c1 * t2m) + geopotential * (c2 * st2^2 - c3 * st2 + c4) + (c5 - c6 * snow_cover) * (
            vₕ * (
                c7 * laiₕ - c8 * geopotential
            ) + sm2 - c9 * t2m * (-(c10 * st2) + sm2) * (vₗ + 2 * sm2 - c11)
        ) + c12
    )

    # Std deviation, in this case a constant
    σ = c13
    return μ + σ * noise
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

abstract type AbstractOceanAlbedo end

@parameterized @kwdef struct OceanDirectReflection{NF} <: AbstractOceanAlbedo

    "Regression coefficient from Jin et al."
    @param c0::NF = 0.0152

    "Regression coefficient from Jin et al."
    @param c1::NF = -1.7873

    "Regression coefficient from Jin et al."
    @param c2::NF = 6.8972

    "Regression coefficient from Jin et al."
    @param c3::NF = -8.5778

    "Regression coefficient from Jin et al."
    @param c4::NF = 4.071

    "Regression coefficient from Jin et al."
    @param c5::NF = -7.6446

    "Regression coefficient from Jin et al."
    @param c6::NF = 0.1643

    "Regression coefficient from Jin et al."
    @param c7::NF = -7.8409

    "Regression coefficient from Jin et al."
    @param c8::NF = -3.5639

    "Regression coefficient from Jin et al."
    @param c9::NF = -2.3588

    "Regression coefficient from Jin et al."
    @param c10::NF = 10.0538
end

Adapt.@adapt_structure OceanDirectReflection
OceanDirectReflection(SG::SpectralGrid; kwargs...) = OceanDirectReflection{SG.NF}(; kwargs...)

"""
    fresnel_reflectance(n, μ)

Unpolarized Fresnel reflectance of a flat air–water interface for light incident from air
with cosine of incidence angle `μ` onto a medium of relative refractive index `n > 1`.
Average of the perpendicular (s) and parallel (p) polarized reflectances. Returns
((n-1)/(n+1))² at normal incidence (μ = 1) and 1 at grazing incidence (μ = 0)."""
@inline function fresnel_reflectance(n, μ)
    s = sqrt(n^2 - 1 + μ^2)                 # n cos θₜ via Snell's law, θₜ refraction angle
    r_s = (μ - s) / (μ + s)                 # perpendicular (s) polarization
    r_p = (n^2 * μ - s) / (n^2 * μ + s)     # parallel (p) polarization
    return (r_s^2 + r_p^2) / 2
end

@inline function (c::OceanDirectReflection)(p)
    (; μ, σ, n, n₀) = p
    (; c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10) = c

    # Formulation following Jin et al. 2011
    # Albedo is simply the Fresnel reflectance when the surface roughness σ == 0
    # and updates when windy/wavey with function f(μ, σ)

    # Calculate surface roughness dependent (μ, σ)
    term1 = c0 + c1 * μ + c2 * μ^2 + c3 * μ^3 + c4 * σ + c5 * μ * σ
    term2 = exp(c6 + c7 * μ + c8 * μ^2 + c9 * σ + c10 * μ * σ)
    f = term1 * term2

    # Now fresnel reflectances
    r = fresnel_reflectance(n, μ)
    r₀ = fresnel_reflectance(n₀, μ)

    albedo_direct_reflection = r - (r / r₀) * f

    return albedo_direct_reflection
end

@parameterized @kwdef struct OceanDiffuseReflection{NF} <: AbstractOceanAlbedo
    "Regression coefficient from Jin et al."
    @param c0::NF = -0.1479

    "Regression coefficient from Jin et al."
    @param c1::NF = 0.1502

    "Regression coefficient from Jin et al."
    @param c2::NF = -0.0176

end

Adapt.@adapt_structure OceanDiffuseReflection
OceanDiffuseReflection(SG::SpectralGrid; kwargs...) = OceanDiffuseReflection{SG.NF}(; kwargs...)

@inline function (c::OceanDiffuseReflection)(p)
    (; σ, n) = p
    (; c0, c1, c2) = c

    # Formulation following Jin et al. 2011
    # Using the cloudy diffuse albedo parameterisation instead of isotropic diffuse
    # since we define direct/diffuse from cloud cover fraction in SpeedyWeather
    albedo_diffuse_reflection = c0 + c1 * n + c2 * n * σ

    return albedo_diffuse_reflection
end

@parameterized @kwdef struct ConstantOceanVolumeAlbedo{NF} <: AbstractOceanAlbedo
    "Broadband volume-scattering albedo for case 1 waters, Jin et al. eq. (15)"
    @param albedo::NF = 0.006 (bounds = 0 .. 1,)
end

@inline (c::ConstantOceanVolumeAlbedo)(p, μ) = c.albedo

Adapt.@adapt_structure ConstantOceanVolumeAlbedo
ConstantOceanVolumeAlbedo(SG::SpectralGrid; kwargs...) = ConstantOceanVolumeAlbedo{SG.NF}(; kwargs...)

@parameterized @kwdef struct SpectralOceanVolumeAlbedo{NF, VectorType} <: AbstractOceanAlbedo
    "[OPTION] Chlorophyll concentration [mg/m³], clamped to ≥ 0.02 (M&M eq. 14 lower limit)"
    @param chl::NF = 0.1 (bounds = Positive,)

    "Regression coefficient from Jin et al."
    @param r0::NF = 0.4817

    "Regression coefficient from Jin et al."
    @param r1::NF = 0.0149

    "Regression coefficient from Jin et al."
    @param r2::NF = 0.207

    "Regression coefficient from Jin et al."
    @param β0::NF = 0.6279

    "Regression coefficient from Jin et al."
    @param β1::NF = 0.2227

    "Regression coefficient from Jin et al."
    @param β2::NF = 0.0513

    "Regression coefficient from Jin et al."
    @param β3::NF = 0.2465

    "Regression coefficient from Jin et al."
    @param β4::NF = 0.3119

    "Coefficient from Morel & Maritorena, eqn (13)"
    @param b0::NF = 0.002

    "Coefficient from Morel & Maritorena, eqn (13)"
    @param b1::NF = 0.01

    "Coefficient from Morel & Maritorena, eqn (13)"
    @param b2::NF = 0.5

    "Coefficient from Morel & Maritorena, eqn (13)"
    @param b3::NF = 0.25

    "Coefficient from Morel & Maritorena, eqn (13), reference wavelength [nm]"
    @param b4::NF = 550

    "Coefficient from Morel & Maritorena, eqn (12)"
    @param b5::NF = 0.416

    "Coefficient from Morel & Maritorena, eqn (12)"
    @param b6::NF = 0.766

    "Coefficient from Morel & Maritorena, eqn (14)"
    @param b7::NF = 0.3

    "Coefficient from Morel & Maritorena, eqn (16)"
    @param a0::NF = 0.06

    "Coefficient from Morel & Maritorena, eqn (16)"
    @param a1::NF = 0.65

    "Coefficient from Morel & Maritorena, eqn (17)"
    @param a2::NF = 0.014

    "Coefficient from Morel & Maritorena, eqn (17), reference wavelength [nm]"
    @param a3::NF = 440

    "Coefficient from Morel & Maritorena, eqn (16)"
    @param a4::NF = 0.2

    "Pure water absorption at 440 nm in the CDOM term, Pope & Fry (1997) as in Jin et al.'s code [m⁻¹]"
    @param a_w_440::NF = 0.00635

    "Band centre wavelengths [nm]"
    wavelengths::VectorType

    "Per-band R₀ intercept A in R₀ = A + Bμ, computed in initialize!"
    R₀_intercept::VectorType

    "Per-band R₀ slope B in R₀ = A + Bμ, computed in initialize!"
    R₀_slope::VectorType
end

Adapt.@adapt_structure SpectralOceanVolumeAlbedo
function SpectralOceanVolumeAlbedo(SG::SpectralGrid; wavelengths = [450, 550, 650], kwargs...)
    (; NF, architecture) = SG
    nbands = length(wavelengths)
    on_arch(x) = on_architecture(architecture, x)
    λs = on_arch(NF.(wavelengths))
    return SpectralOceanVolumeAlbedo{NF, typeof(λs)}(;
        wavelengths = λs,
        R₀_intercept = on_arch(zeros(NF, nbands)),
        R₀_slope = on_arch(zeros(NF, nbands)),
        kwargs...,
    )
end

"""$(TYPEDSIGNATURES) Linear interpolation of tabulated `ys` at ascending `xs`,
clamped to the end values outside the table (no extrapolation)."""
@inline function interpolate_table(xs::NTuple{N}, ys::NTuple{N}, x) where {N}
    x <= xs[1] && return oftype(x, ys[1])
    x >= xs[N] && return oftype(x, ys[N])
    lo, hi = 1, N                       # binary search for xs[lo] ≤ x < xs[hi]
    while hi - lo > 1
        mid = (lo + hi) ÷ 2
        xs[mid] <= x ? (lo = mid) : (hi = mid)
    end
    t = (x - xs[lo]) / (xs[hi] - xs[lo])
    return oftype(x, (1 - t) * ys[lo] + t * ys[hi])
end

# Normalized phytoplankton absorption A_chl(λ), A_chl(440) = 1, from Jin et al. (2011) supplementary code
const ACHL_WAVELENGTHS = Tuple(350.0:10.0:700.0)
const ACHL = (
    0.56, 0.57, 0.59, 0.61, 0.66, 0.687, 0.828, 0.913, 0.973, 1.0, 0.944, 0.917,
    0.87, 0.798, 0.75, 0.668, 0.618, 0.528, 0.474, 0.416, 0.357, 0.294, 0.276, 0.291,
    0.282, 0.236, 0.252, 0.276, 0.317, 0.334, 0.356, 0.441, 0.595, 0.502, 0.329, 0.215,
)

# Pure water absorption a_w(λ) [m⁻¹], from Jin et al. (2011) supplementary code:
# Smith & Baker (1981) for λ ≤ 510 nm, Pope & Fry (1997) for 512.5–727.5 nm
const AW_WAVELENGTHS = (Tuple(200.0:10.0:510.0)..., Tuple(512.5:2.5:730.0)..., Tuple(740.0:10.0:800.0)...)
const AW = (
    3.07, 1.99, 1.31, 0.927, 0.72, 0.559, 0.457, 0.373,
    0.288, 0.215, 0.141, 0.105, 0.0844, 0.0678, 0.0561, 0.0463,
    0.0379, 0.03, 0.022, 0.0191, 0.0171, 0.0162, 0.0153, 0.0144,
    0.0145, 0.0145, 0.0156, 0.0156, 0.0176, 0.0196, 0.0257, 0.0357,
    0.0372, 0.0396, 0.0399, 0.0409, 0.0416, 0.0417, 0.0428, 0.0434,
    0.0447, 0.0452, 0.0466, 0.0474, 0.0489, 0.0511, 0.0537, 0.0565,
    0.0593, 0.0596, 0.0606, 0.0619, 0.064, 0.0642, 0.0672, 0.0695,
    0.0733, 0.0772, 0.0836, 0.0896, 0.0989, 0.11, 0.122, 0.1351,
    0.1516, 0.1672, 0.1925, 0.2224, 0.247, 0.2577, 0.2629, 0.2644,
    0.2665, 0.2678, 0.2707, 0.2755, 0.281, 0.2834, 0.2904, 0.2916,
    0.2995, 0.3012, 0.3077, 0.3108, 0.322, 0.325, 0.335, 0.34,
    0.358, 0.371, 0.393, 0.41, 0.424, 0.429, 0.436, 0.439,
    0.448, 0.448, 0.461, 0.465, 0.478, 0.486, 0.502, 0.516,
    0.538, 0.559, 0.592, 0.624, 0.663, 0.704, 0.756, 0.827,
    0.914, 1.007, 1.119, 1.231, 1.356, 1.489, 1.678, 1.799,
    2.38, 2.47, 2.55, 2.51, 2.36, 2.16, 2.07,
)

@inline A_chl(λnm) = interpolate_table(ACHL_WAVELENGTHS, ACHL, λnm)
@inline a_w(λnm) = interpolate_table(AW_WAVELENGTHS, AW, λnm)

# From Morel (1974) and Bricaud (1998)
@inline b_w(λnm) = oftype(λnm, 0.00288) * (λnm / 500)^oftype(λnm, -4.32)

"""$(TYPEDSIGNATURES) Coefficients A, B of the subsurface irradiance reflectance R₀ = A + Bμ
at wavelength `λ` [nm] for chlorophyll `chl` [mg/m³], Jin et al. eqs. (8)–(9) with
Morel & Maritorena (2001) eqs. (11)–(18). Called once per band at initialization."""
function R₀_coefficients(c::SpectralOceanVolumeAlbedo, λ, chl)
    (; β0, β1, β2, β3, β4, b0, b1, b2, b3, b4, b5, b6, b7, a0, a1, a2, a3, a4, a_w_440) = c
    chl = max(chl, oftype(chl, 0.02))           # lower limit of M&M eq. (14)

    # Backscattering, M&M eqs. (11)–(14)
    v = chl > 2 ? zero(chl) : (log10(chl) - b7) / 2
    bbp = (b0 + b1 * (b2 - b3 * log10(chl)) * (λ / b4)^v) * b5 * chl^b6
    bbw = b_w(λ) / 2
    bb = bbw + bbp
    η_b = bbw / bb

    # Absorption, M&M eqs. (16)–(18), with A_chl(440) = 1
    # NOTE: Jin et al.'s supplementary code uses Smith & Baker (1981) a_w below 510 nm but
    # the Pope & Fry (1997) value a_w(440) = 0.00635 in the CDOM term. Kept for parity;
    # Pope & Fry throughout would be consistent with Morel & Maritorena (2001).
    a_y440 = a4 * (a_w_440 + a0 * chl^a1)
    a = a_w(λ) + a0 * A_chl(λ) * chl^a1 + a_y440 * exp(-a2 * (λ - a3))

    # β = β0 - β1η - β2η² + (β3η - β4)μ is linear in μ, so R₀ = A + Bμ
    A = (β0 - β1 * η_b - β2 * η_b^2) * bb / a
    B = (β3 * η_b - β4) * bb / a
    return A, B
end

function initialize!(scheme::SpectralOceanVolumeAlbedo, model::PrimitiveEquation)
    λs = Array(scheme.wavelengths)              # compute on CPU
    A = similar(λs)
    B = similar(λs)
    for k in eachindex(λs)
        A[k], B[k] = R₀_coefficients(scheme, λs[k], scheme.chl)
    end
    copyto!(scheme.R₀_intercept, A)             # works for CPU and GPU arrays
    copyto!(scheme.R₀_slope, B)
    return nothing
end

@inline function (c::SpectralOceanVolumeAlbedo)(p, μ)
    (; σ, αₛ, band) = p
    (; r0, r1, r2) = c

    r_w = r0 - r1 * σ - r2 * σ^2                                # Jin eq. (7)
    R₀ = c.R₀_intercept[band] + c.R₀_slope[band] * μ            # Jin eqs. (8)–(9), precomputed
    return R₀ * (1 - r_w) * (1 - αₛ) / (1 - r_w * R₀)           # Jin eq. (6)
end

export WhitecapFoam, NoFoam

"""Whitecap (foam) correction following Koepke (1984) as used by Jin et al. (2011) eqs. (16)–(17):
constant foam albedo and fractional whitecap coverage f = a * w^b from wind speed w [m/s].
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct WhitecapFoam{NF} <: AbstractOceanAlbedo
    "[OPTION] Albedo of whitecaps/foam [1], Koepke (1984)"
    @param albedo_foam::NF = 0.55 (bounds = 0 .. 1,)

    "[OPTION] Whitecap coverage coefficient [(s/m)^b], Jin et al. eq. (16)"
    @param coverage_coefficient::NF = 2.95e-6 (bounds = Positive,)

    "[OPTION] Whitecap coverage wind speed exponent [1], Jin et al. eq. (16)"
    @param coverage_exponent::NF = 3.52 (bounds = Positive,)
end

Adapt.@adapt_structure WhitecapFoam
WhitecapFoam(SG::SpectralGrid; kwargs...) = WhitecapFoam{SG.NF}(; kwargs...)

@inline function whitecap_coverage(foam::WhitecapFoam, w)
    (; coverage_coefficient, coverage_exponent) = foam
    # Fractional whitecap coverage in [0, 1] at wind speed `w` [m/s], Jin et al. eq. (16).
    return min(coverage_coefficient * w^coverage_exponent, one(w))
end

@inline function (foam::WhitecapFoam)(α, w)
    # Albedo `α` corrected for whitecaps at wind speed `w` [m/s], Jin et al. eq. (17).
    f = whitecap_coverage(foam, w)
    return f * foam.albedo_foam + (1 - f) * α
end

"""No whitecap correction."""
@parameterized struct NoFoam <: AbstractOceanAlbedo end
Adapt.@adapt_structure NoFoam
NoFoam(SG::SpectralGrid) = NoFoam()
@inline (::NoFoam)(α, w) = α

export JinOceanAlbedo, ConstantOceanVolumeAlbedo, SpectralOceanVolumeAlbedo

@parameterized @kwdef struct JinOceanAlbedo{NF, V, F} <: AbstractAlbedo
    "Ratio of air-water refractive indices"
    @param n::NF = 1.34

    "Refractive index of water in visible spectrum"
    @param n₀::NF = 1.34

    "Effective cosine for diffuse incidence, Jin et al. eq. (10)"
    @param μₑ::NF = 0.676

    "[OPTION] Albedo over sea ice at concentration=1 [1]"
    @param albedo_ice::NF = 0.6 (bounds = 0 .. 1,)

    @component ocean_direct_reflection::OceanDirectReflection{NF} = OceanDirectReflection{NF}()

    @component ocean_diffuse_reflection::OceanDiffuseReflection{NF} = OceanDiffuseReflection{NF}()

    @component ocean_volume::V = ConstantOceanVolumeAlbedo{NF}()

    @component foam::F = WhitecapFoam{NF}()
end

Adapt.@adapt_structure JinOceanAlbedo
function JinOceanAlbedo(
        SG::SpectralGrid;
        ocean_volume = ConstantOceanVolumeAlbedo(SG),
        foam = WhitecapFoam(SG),
        kwargs...
    )
    ocean_volume isa SpectralOceanVolumeAlbedo &&
        throw(ArgumentError("SpectralOceanVolumeAlbedo requires a spectral shortwave scheme"))
    return JinOceanAlbedo{SG.NF, typeof(ocean_volume), typeof(foam)}(; ocean_volume, foam, kwargs...)
end


function variables(::JinOceanAlbedo)
    return (
        ParameterizationVariable(:black_sky_albedo, GridXY(), desc = "Black-sky ocean albedo", units = "1", namespace = :ocean),
        ParameterizationVariable(:white_sky_albedo, GridXY(), desc = "White-sky ocean albedo", units = "1", namespace = :ocean),
    )
end

function initialize!(albedo::JinOceanAlbedo, model::PrimitiveEquation)
    initialize!(albedo.ocean_volume, model)
    return nothing
end

initialize!(::ConstantOceanVolumeAlbedo, ::PrimitiveEquation) = nothing

@propagate_inbounds function albedo!(ij, albedo, vars, scheme::JinOceanAlbedo, model)
    (; ocean) = vars.parameterizations
    (; n, n₀, μₑ) = scheme
    NF = eltype(albedo)

    # Cosine of solar zenith angle, only when the sun is up
    μ = clamp(vars.parameterizations.cos_zenith[ij], zero(NF), one(NF))
    w = vars.parameterizations.surface_wind_speed[ij]
    σ = sqrt(NF(0.003) + NF(0.00512) * w)     # Cox & Munk, Jin et al. eq. (2)

    p = (; μ, σ, n, n₀)

    # Surface (Fresnel) reflection, Jin eqs. (1) and (5b)
    αs_dir = scheme.ocean_direct_reflection(p)
    αs_dif = scheme.ocean_diffuse_reflection(p)

    # Ocean volume scattering, Jin eqs. (6) and (10)
    αw_dir = scheme.ocean_volume((; p..., αₛ = αs_dir), μ)
    αw_dif = scheme.ocean_volume((; p..., αₛ = αs_dif), μₑ)

    # Direct (black-sky) and diffuse (white-sky) albedo, Jin eqs. (11)–(12),
    # blended by the shortwave scheme via combine_surface_albedo!
    black_sky = αs_dir + αw_dir
    white_sky = αs_dif + αw_dif

    # Whitecap (foam) correction, Jin et al. eqs. (16)–(17)
    black_sky = scheme.foam(black_sky, w)
    white_sky = scheme.foam(white_sky, w)

    ℵ = haskey(vars.prognostic.ocean, :sea_ice_concentration) ? vars.prognostic.ocean.sea_ice_concentration[ij] : zero(NF)
    black_sky += ℵ * (scheme.albedo_ice - black_sky)
    white_sky += ℵ * (scheme.albedo_ice - white_sky)

    ocean.black_sky_albedo[ij] = black_sky
    ocean.white_sky_albedo[ij] = white_sky

    # Fallback for radiation schemes with no direct/diffuse radiation fraction
    albedo[ij] = white_sky

    return nothing
end


"""$(TYPEDSIGNATURES)
Blend direct (black-sky) and diffuse (white-sky) albedo into the broadband `albedo` of
`surface` (a namespace of `vars.parameterizations`, e.g. `.land`) given the fraction of
direct shortwave radiation `fraction_direct`. Does nothing if `surface` does not define
both, decided at compile time."""
@propagate_inbounds function combine_surface_albedo!(ij, surface::NamedTuple, fraction_direct)
    if haskey(surface, :black_sky_albedo) && haskey(surface, :white_sky_albedo)
        black_sky = surface.black_sky_albedo[ij]
        white_sky = surface.white_sky_albedo[ij]
        surface.albedo[ij] = fraction_direct * black_sky + (1 - fraction_direct) * white_sky
    end
    return nothing
end

"""$(TYPEDSIGNATURES) Blend black/white-sky albedos over ocean and land wherever defined."""
@propagate_inbounds function combine_albedos!(ij, vars, fraction_direct)
    combine_surface_albedo!(ij, vars.parameterizations.ocean, fraction_direct)
    combine_surface_albedo!(ij, vars.parameterizations.land, fraction_direct)
    return nothing
end
