"""
Albedo structs should be defined as MyAlbedo <: AbstractAlbedo and implement the following interface:

- `initialize!(albedo::MyAlbedo, model)`
- `albedo!(ij, albedofield, vars::Variables, albedo::MyAlbedo, model)`

`initialize!` is as for every model component but in contrast to other parameterizations albedos
should not implement `parameterization!` but rather `albedo!`, which is called from within `parameterization!`.
This is because every albedo needs to be used either as part of OceanLandAlbedo or the same albedo
is called over both land and ocean surfaces. If a surface namespace (vars.parameterizations.land or .ocean) 
contains both black_sky_albedo and white_sky_albedo, a shortwave scheme that knows the direct-beam fraction 
blends them into that namespace's albedo. Otherwise it leaves albedo alone."""
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
        snow_depth = get_prognostic_step(vars.prognostic.land.snow_depth, ts, DummyParameterization())
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

@parameterized @kwdef struct LearnedLandAlbedo{NF, Scheme <: AbstractSnowCover} <: AbstractAlbedo
    # Put symbolic regression coefficients in here
    "Visible wavelength weight to the broadband shortwave albedo"
    @param vis_weight::NF = 0.5143 (bounds = 0 .. 1,)

    "Snow depth scaling for calculation snow cover fraction from snow depth"
    @param snow_depth_scale::NF = 0.05 (bounds = Positive,)

    "Snow cover-albedo scheme"
    @param snow_cover::Scheme = SaturatingSnowCover() (group = :snow_cover,)

    "[OPTION] Switch on stochastic output for the learned albedo"
    noise::Bool = false
end

Adapt.@adapt_structure LearnedLandAlbedo
function LearnedLandAlbedo(
        SG::SpectralGrid;
        snow_cover = SaturatingSnowCover(),
        kwargs...,
    )
    return LearnedLandAlbedo{SG.NF, typeof(snow_cover)}(; snow_cover, kwargs...)
end

function variables(::LearnedLandAlbedo)
    return (
        ParameterizationVariable(:black_sky_albedo, GridXY(), desc = "Black-sky land albedo", units = "1", namespace = :land),
        ParameterizationVariable(:white_sky_albedo, GridXY(), desc = "White-sky land albedo", units = "1", namespace = :land),
    )
end

initialize!(albedo::LearnedLandAlbedo, model::PrimitiveEquation) = nothing

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

@inline _safe_sqrt(x) = x >= zero(x) ? sqrt(x) : oftype(x, NaN)

@inline function visible_isotropic_scattering(vₕ, vₗ, laiₗ, sm2, st2, snow_cover, noise)
    T = promote_type(Float32, typeof(vₕ), typeof(vₗ), typeof(laiₗ), typeof(sm2), typeof(st2), typeof(snow_cover), typeof(noise))
    return (((noise * T(0.03896045440440107) + T(-0.1218158966586572)) * T(0.025728017593462257)) * snow_cover + T(-0.1346102784222333)) * vₕ + (((T(0.006908593126432138) * snow_cover + T(-0.13673067396327065)) + (T(0.683619830705067) + ((T(-0.0009017510759633551) * st2 + noise * T(0.0285256619975839)) + T(-0.1346102784222333) * sm2))) + T(-0.052392875442762124) * (vₗ * ((T(0.683619830705067) + laiₗ) + (T(-0.14487101003784875) + T(0.11589160767036352) * noise) * (T(-0.13673067396327065) * snow_cover)) + T(0.7161992148232291)))
end

@inline function visible_volumetric_scattering(vₕ, vₗ, laiₗ, sm1, st1, st2, geopotential, snow_cover, t2m, noise)
    T = promote_type(Float32, typeof(vₕ), typeof(vₗ), typeof(laiₗ), typeof(sm1), typeof(st1), typeof(st2), typeof(geopotential), typeof(snow_cover), typeof(t2m), typeof(noise))
    return (T(0.016140250524517893) * noise + T(-0.04410325379601232) * sm1) + ((((vₕ + vₗ) * (T(-0.01470440805972792) + (T(9.122606489165302e-5) * (t2m * T(0.008546824863390104) + noise)) * ((((T(0.0025200785769508103) * geopotential) * (geopotential * T(0.0025200785769508103))) * T(0.016140250524517893) + (T(0.005668214606309846) * snow_cover) * (((t2m - snow_cover) - snow_cover) * (st2 * ((t2m * T(0.008546824863390104)) * T(0.005668214606309846))) + T(0.016140250524517893) * geopotential)) / ((vₕ + T(0.0025200785769508103)) + vₗ))) + laiₗ * T(-0.005123590195731909)) + st1 * T(0.0014159332530380526)) + t2m * T(-0.0011845617711790535))
end

@inline function visible_geometric_scattering(vₕ, laiₗ, sm1, geopotential, snow_cover, noise)
    T = promote_type(Float32, typeof(vₕ), typeof(laiₗ), typeof(sm1), typeof(geopotential), typeof(snow_cover), typeof(noise))
    return T(0.013754593690262884) + (_safe_sqrt(_safe_sqrt(laiₗ)) * (T(9.906885493797509e-7) * geopotential + (T(0.006595930444019463) + (sm1 + (geopotential * laiₗ) * T(2.1805730417085534e-5)) * T(-0.018464047185435716))) + (((geopotential * T(4.3663667100446695e-7)) * snow_cover - T(0.04395703327526142)) * ((((laiₗ - vₕ * geopotential) * T(2.1805730417085534e-5) + (laiₗ - (vₕ - T(0.03334265669607637)) * (T(-2.161750513690768) + laiₗ))) - T(0.04395703327526142)) * (T(0.03334265669607637) + noise * T(0.019535615884723058))) + (noise + sm1) * T(0.008867523895391013)))
end

@inline function near_infrared_isotropic_scattering(vₕ, vₗ, laiₕ, sm2, snow_cover, t2m, noise)
    T = promote_type(Float32, typeof(vₕ), typeof(vₗ), typeof(laiₕ), typeof(sm2), typeof(snow_cover), typeof(t2m), typeof(noise))
    return T(-0.08235648263453457) * (_safe_sqrt(vₗ) + (snow_cover / t2m + T(-0.1892657831005156) * vₕ) * (T(-0.07359907041268714) * (noise * laiₕ + snow_cover))) + (((T(-0.09496896272780166) * vₕ + sm2 * (T(-0.741150547535996) + sm2)) + T(0.4934449209685155)) + noise * T(0.045780621419124166))
end

@inline function near_infrared_volumetric_scattering(vₗ, sm1, sm2, st1, st2, geopotential, snow_cover, t2m, noise)
    T = promote_type(Float32, typeof(vₗ), typeof(sm1), typeof(sm2), typeof(st1), typeof(st2), typeof(geopotential), typeof(snow_cover), typeof(t2m), typeof(noise))
    return (((((T(-0.0003633530764346828) * geopotential) * vₗ + t2m) * ((t2m * t2m) * (st1 * T(-1.5726052867019688e-11) + (snow_cover + sm1) * T(6.961606055872098e-10)) + snow_cover * T(-0.00014194486565677876)) + sm2 * ((T(-0.021328337298982213) * st2 + (st1 * T(7.673078268962595e-5)) * st1) * T(0.3118690306291605))) + T(0.02426093826350491) * snow_cover) + T(0.21359294919077146)) + T(0.03500057089164972) * noise
end

@inline function near_infrared_geometric_scattering(vₕ, vₗ, laiₕ, sm2, st1, st2, geopotential, snow_cover, noise)
    T = promote_type(Float32, typeof(vₕ), typeof(vₗ), typeof(laiₕ), typeof(sm2), typeof(st1), typeof(st2), typeof(geopotential), typeof(snow_cover), typeof(noise))
    return (st2 * ((geopotential * T(3.818568357356443e-10)) * st2 + (T(-0.008433043494284033) + T(0.008482159093250554)) * noise) + ((snow_cover * T(-0.0005895508304773011) + T(0.02588756791821848)) * ((((T(-0.0005895508304773011) * (st2 + st2)) * ((sm2 + (st2 + st2) * T(-0.0005895508304773011)) + (snow_cover * (T(-3.7415065307307705e-5) * geopotential) + snow_cover) * T(-0.008433043494284033))) * (((geopotential * T(7.820596733390864e-5) + snow_cover * T(-0.008433043494284033)) + vₗ) + (sm2 + vₗ)) + geopotential * (T(-3.7415065307307705e-5) * vₕ)) + (((T(0.02588756791821848) * laiₕ) * laiₕ) * vₕ + sm2)) + ((T(0.08179264552082091) + (st2 * (geopotential * T(-1.7057584346842938e-7)) + T(0.008482159093250554) * vₗ)) + st1 * T(-0.00020412859377729543)))) + (vₕ * T(0.009331484441448655) + geopotential * T(1.853145610962989e-5))
end

@propagate_inbounds function albedo!(ij, albedo, vars, scheme::LearnedLandAlbedo, model)
    # Unpack predictors for albedo parameter prediction
    (; land) = vars.parameterizations
    field_capacity = model.land.thermodynamics.field_capacity

    NF = eltype(albedo)
    zero_value = zero(NF)

    noise = zero_value
    t2m = vars.parameterizations.surface_air_temperature[ij]
    vₕ = haskey(land, :high_vegetation_cover) ? land.high_vegetation_cover[ij] : zero_value
    vₗ = haskey(land, :low_vegetation_cover) ? land.low_vegetation_cover[ij] : zero_value
    laiₕ = haskey(land, :high_vegetation_leaf_area) ? land.high_vegetation_leaf_area[ij] : zero_value
    laiₗ = haskey(land, :low_vegetation_leaf_area) ? land.low_vegetation_leaf_area[ij] : zero_value
    geopotential = model.planet.gravity * model.orography.orography[ij]

    # Get snow depth and calculate snow cover
    ts = time_stepper(model.time_stepping, :land)
    if haskey(vars.prognostic, :land) && haskey(vars.prognostic.land, :snow_depth)
        ts = time_stepper(model.time_stepping, :land)
        snow_depth = get_prognostic_step(vars.prognostic.land.snow_depth, ts, DummyParameterization())

        # how to compute snow cover from snow depth
        snow_cover_scheme = scheme.snow_cover

        # compute snow-cover fraction using the chosen scheme and convert to percentage
        # as the learned albedos use snow-cover in percentage
        snow_cover = NF(100) * snow_cover_scheme(snow_depth[ij], scheme.snow_depth_scale)
    else
        snow_cover = zero_value
    end

    # Get soil moisture at upper and lower layers
    if haskey(vars.prognostic.land, :soil_moisture)
        soil_moisture = get_prognostic_step(
            vars.prognostic.land.soil_moisture, ts, scheme,
        )
        capacity = land.soil_moisture_capacity[ij]
        
        # Multiply by field capacity to get back volume of water in the layer
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

    # Calculate parameters!
    vis_iso = visible_isotropic_scattering(vₕ, vₗ, laiₗ, sm2, st2, snow_cover, noise)
    nir_iso = near_infrared_isotropic_scattering(vₕ, vₗ, laiₕ, sm2, snow_cover, t2m, noise)

    vis_vol = visible_volumetric_scattering(vₕ, vₗ, laiₗ, sm1, st1, st2, geopotential, snow_cover, t2m, noise)
    nir_vol = near_infrared_volumetric_scattering(vₗ, sm1, sm2, st1, st2, geopotential, snow_cover, t2m, noise)

    vis_geo = visible_geometric_scattering(vₕ, laiₗ, sm1, geopotential, snow_cover, noise)
    nir_geo = near_infrared_geometric_scattering(vₕ, vₗ, laiₕ, sm2, st1, st2, geopotential, snow_cover, noise)

    # Near infrared parameter weight
    nir_weight = 1 - scheme.vis_weight

    # Weight vis/nir parameters to get the shortwave broadband albedo
    sw_iso = (vis_iso * scheme.vis_weight) + (nir_iso * nir_weight)
    sw_vol = (vis_vol * scheme.vis_weight) + (nir_vol * nir_weight)
    sw_geo = (vis_geo * scheme.vis_weight) + (nir_geo * nir_weight)

    # Get zenith angle in radians
    μ = vars.parameterizations.cos_zenith[ij]
    θ = acos(clamp(μ, zero(μ), one(μ))) # only need when the sun is up

    # Then from those we'll calculate the black and white sky albedos
    black_sky = black_sky_albedo(sw_iso, sw_vol, sw_geo, θ)
    white_sky = white_sky_albedo(sw_iso, sw_vol, sw_geo)
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
