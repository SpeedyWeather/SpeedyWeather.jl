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

"""Some albedos rely on scattering computations, they will then have subcomponents of these types."""
abstract type AbstractScattering <: AbstractModelComponent end

# all albedos need to define albedo fields for ocean and land separately and total
function variables(::AbstractAlbedo)
    return (
        ParameterizationVariable(:albedo, Grid2D(), desc = "Albedo", units = "1"),
        ParameterizationVariable(:albedo, Grid2D(), desc = "Albedo over the ocean", units = "1", namespace = :ocean),
        ParameterizationVariable(:albedo, Grid2D(), desc = "Albedo over the land", units = "1", namespace = :land),
    )
end

# Every radiation scheme should call combine_albedos! to blend white and black sky
# depending on the fraction of direct/diffusive radiation (depends on clouds)
"""$(TYPEDSIGNATURES) Blend black/white-sky albedos over ocean and land wherever defined."""
@propagate_inbounds function combine_albedos!(ij, vars, fraction_direct)
    combine_surface_albedo!(ij, vars.parameterizations.ocean, fraction_direct)
    combine_surface_albedo!(ij, vars.parameterizations.land, fraction_direct)
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