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
        ocean = OceanSeaIceAlbedo(SG),      # default ocean albedo
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
    θ = atan(sqrt(1 - μ^2), μ)  # = acos(μ) but Reactant-compatible (no acos in XLA)

    # Then from those we'll calculate the black and white sky albedos
    black_sky = black_sky_albedo(isotropic_parameter, volumetric_parameter, geometric_parameter, θ)
    white_sky = white_sky_albedo(isotropic_parameter, volumetric_parameter, geometric_parameter)
    land.black_sky_albedo[ij] = black_sky
    land.white_sky_albedo[ij] = white_sky

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
