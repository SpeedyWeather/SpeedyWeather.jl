abstract type AbstractAlbedo <: AbstractModelComponent end

export OceanLandAlbedo

"""Composite type: Two albedo parameterizations, one for ocean and one for land surfaces.
Fields are $(TYPEDFIELDS)"""
@kwdef struct OceanLandAlbedo{Ocean, Land} <: AbstractAlbedo
    "[OPTION] Albedo parameterization for ocean surfaces"
    ocean::Ocean

    "[OPTION] Albedo parameterization for land surfaces"
    land::Land
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
end

Adapt.@adapt_structure OceanLandAlbedo

export DefaultAlbedo

"""$(TYPEDSIGNATURES)
Default albedo parameterization with `OceanSeaIceAlbedo`
for ocean and `AlbedoClimatology` for land surfaces."""
function DefaultAlbedo(SG::SpectralGrid;
    ocean = OceanSeaIceAlbedo(SG),
    land = LandSnowAlbedo(SG))
    return OceanLandAlbedo(ocean, land)
end

function initialize!(albedo::OceanLandAlbedo, model::PrimitiveEquation)
    initialize!(albedo.ocean, model)
    initialize!(albedo.land, model)
end

# composite OceanLandAlbedo: call separately for ocean and land with .ocean and .land
@propagate_inbounds function parameterization!(ij, diagn::DiagnosticVariables, progn, albedo::OceanLandAlbedo, model)
    parameterization!(ij, diagn.physics.ocean, progn, albedo.ocean, model)
    parameterization!(ij, diagn.physics.land, progn, albedo.land, model)
end

function variables(::OceanLandAlbedo)
    return (
        DiagnosticVariable(name=:albedo, dims=Grid2D(), desc="Albedo", units="1"),
        DiagnosticVariable(name=:albedo, dims=Grid2D(), desc="Albedo over the ocean", units="1", namespace=:ocean),
        DiagnosticVariable(name=:albedo, dims=Grid2D(), desc="Albedo over the land", units="1", namespace=:land),
    )
end

# single albedo: call separately for ocean and land with the same albedo
 @propagate_inbounds function parameterization!(ij, diagn::DiagnosticVariables, progn, albedo::AbstractAlbedo, model)
    parameterization!(ij, diagn.physics.ocean, progn, albedo, model)
    parameterization!(ij, diagn.physics.land, progn, albedo, model)
end

## GLOBAL CONSTANT ALBEDO
export GlobalConstantAlbedo

"""Global constant albedo parameterization. To be used for land and ocean 
or only one of them within a `OceanLandAlbedo`.
Fields are $(TYPEDFIELDS)"""
@kwdef struct GlobalConstantAlbedo{NF} <: AbstractAlbedo
    "[OPTION] Albedo value [1]"
    albedo::NF = 0.3
end

GlobalConstantAlbedo(SG::SpectralGrid; kwargs...) = GlobalConstantAlbedo{SG.NF}(; kwargs...)
initialize!(albedo::GlobalConstantAlbedo, ::PrimitiveEquation) = nothing
@propagate_inbounds parameterization!(ij, diagn, progn, albedo::GlobalConstantAlbedo, model) =
    albedo!(ij, diagn.albedo, albedo.albedo)
@propagate_inbounds function albedo!(ij, diagn_albedo::AbstractArray, albedo::Real)
    diagn_albedo[ij] = albedo
end

Adapt.@adapt_structure GlobalConstantAlbedo

## MANUAL ALBEDO
export ManualAlbedo

"""Manual albedo field, to be used with `set!` and is copied into the diagnostic variables on every time step.
Defined so that parameterizations can change the albedo at every time step (e.g. snow cover) without
losing the information of the original surface albedo. Fields are
$(TYPEDFIELDS)"""
struct ManualAlbedo{GridVariable2D} <: AbstractAlbedo
    "Albedo field [1]"
    albedo::GridVariable2D
end

ManualAlbedo(SG::SpectralGrid) = ManualAlbedo{SG.GridVariable2D}(zeros(SG.GridVariable2D, SG.grid))
initialize!(albedo::ManualAlbedo, model::PrimitiveEquation) = nothing
@propagate_inbounds parameterization!(ij, diagn, progn, albedo::ManualAlbedo, model) = albedo!(ij, diagn.albedo, albedo.albedo)
@propagate_inbounds function albedo!(ij, diagn_albedo::AbstractArray, albedo::AbstractArray)
    diagn_albedo[ij] = albedo[ij]
end

Adapt.@adapt_structure ManualAlbedo

## ALBEDO CLIMATOLOGY
export AlbedoClimatology

"""Albedo climatology loaded from netcdf file.
Fields are $(TYPEDFIELDS)"""
@kwdef struct AlbedoClimatology{GridVariable2D} <: AbstractAlbedo
    "[OPTION] path to the folder containing the albedo file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "[OPTION] filename of albedo"
    file::String = "albedo.nc"

    "[OPTION] variable name in netcdf file"
    varname::String = "alb"

    "[OPTION] Grid the albedo file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "Albedo climatology"
    albedo::GridVariable2D
end

# For GPU usage just discard the extra information and treat it as a `ManualAlbedo`
Adapt.adapt_structure(to, albedo::AlbedoClimatology) = adapt(to, ManualAlbedo(albedo.albedo))

function AlbedoClimatology(SG::SpectralGrid; kwargs...)
    (; GridVariable2D, grid) = SG
    albedo = zeros(GridVariable2D, grid)
    AlbedoClimatology{GridVariable2D}(; albedo, kwargs...)
end

# set albedo with grid, scalar, function; just define path `albedo.albedo` to grid here
set!(albedo::AbstractAlbedo, args...; kwargs...) = set!(albedo.albedo, args...; kwargs...)

function initialize!(albedo::AlbedoClimatology, model::PrimitiveEquation)

    # LOAD NETCDF FILE
    if albedo.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__, "../../input_data", albedo.file)
    else
        path = joinpath(albedo.path, albedo.file)
    end
    ncfile = NCDataset(path)

    a = on_architecture(model.architecture, albedo.file_Grid(ncfile[albedo.varname].var[:, :], input_as=Matrix))
    interpolate!(albedo.albedo, a)
end

@propagate_inbounds parameterization!(ij, diagn, progn, albedo::AlbedoClimatology, model) = albedo!(ij, diagn.albedo, albedo.albedo)

## OceanSeaIceAlbedo
export OceanSeaIceAlbedo

"""Albedo that scales linearly between ocean and ice albedo depending on sea ice concentration.
Fields are $(TYPEDFIELDS)"""
@kwdef struct OceanSeaIceAlbedo{NF} <: AbstractAlbedo
    "[OPTION] Albedo over open ocean [1]"
    albedo_ocean::NF = 0.06

    "[OPTION] Albedo over sea ice at concentration=1 [1]"
    albedo_ice::NF = 0.6
end

Adapt.@adapt_structure OceanSeaIceAlbedo
OceanSeaIceAlbedo(SG::SpectralGrid; kwargs...) = OceanSeaIceAlbedo{SG.NF}(;kwargs...)
initialize!(::OceanSeaIceAlbedo, ::PrimitiveEquation) = nothing
@propagate_inbounds parameterization!(ij, diagn, progn, albedo::OceanSeaIceAlbedo, model) = albedo!(ij, diagn.albedo, progn.ocean, albedo)

@propagate_inbounds function albedo!(ij, diagn_albedo::AbstractArray, ocean, albedo::OceanSeaIceAlbedo)
    (; sea_ice_concentration ) = ocean
    (; albedo_ocean, albedo_ice) = albedo

    # set ocean albedo linearly between ocean and ice depending on sea ice concentration
    diagn_albedo[ij] = albedo_ocean + sea_ice_concentration[ij] * (albedo_ice - albedo_ocean)
end

function variables(::AbstractAlbedo)
    return (
        DiagnosticVariable(name=:albedo, dims=Grid2D(), desc="Albedo", units="1"),
    )
end

abstract type AbstractSnowCover end

export LinearSnowCover, SaturatingSnowCover
"""Linear ramp: snow cover grows with snow depth from 0 to 1 at `snow_depth_scale`."""
struct LinearSnowCover <: AbstractSnowCover end
Adapt.@adapt_structure LinearSnowCover

"""Saturating ramp: snow cover grows with snow depth S as `S/(S+scale)`."""
struct SaturatingSnowCover  <: AbstractSnowCover end
Adapt.@adapt_structure SaturatingSnowCover

"""$(TYPEDSIGNATURES) Snow cover fraction for the linear scheme, clamped to 1."""
@inline (::LinearSnowCover)(snow_depth, scale) = min(snow_depth / scale, 1)

"""$(TYPEDSIGNATURES) Snow cover fraction for the saturating scheme."""
@inline (::SaturatingSnowCover)(snow_depth, scale) = snow_depth / (snow_depth + scale)

## LandSnowAlbedo
export LandSnowAlbedo

@kwdef struct LandSnowAlbedo{NF, Scheme <: AbstractSnowCover} <: AbstractAlbedo
    "Albedo of bare land (excluding vegetation) [1]"
    albedo_land::NF = 0.4

    "Albedo of high vegetation [1]"
    albedo_high_vegetation::NF = 0.15

    "Albedo of low vegetation [1]"
    albedo_low_vegetation::NF = 0.20

    "Albedo of snow [1], additive to land"
    albedo_snow::NF = 0.4

    "Conversion from snow depth to snow cover [m]"
    snow_depth_scale::NF = 0.05

    "Snow cover-albedo scheme"
    snow_cover::Scheme = SaturatingSnowCover()
end

Adapt.@adapt_structure LandSnowAlbedo
LandSnowAlbedo(SG::SpectralGrid; snow_cover = SaturatingSnowCover(), kwargs...) =
    LandSnowAlbedo{SG.NF, typeof(snow_cover)}(; snow_cover, kwargs...)

initialize!(albedo::LandSnowAlbedo, model::PrimitiveEquation) = nothing

@propagate_inbounds parameterization!(ij, diagn, progn, albedo::LandSnowAlbedo, model) = 
    albedo!(ij, diagn, progn.land, albedo)

@propagate_inbounds function albedo!(
    ij,                             # grid point index
    diagn,                          # diagnostic physics land variables
    land,                           # land prognostic variables
    albedo_scheme::LandSnowAlbedo,
)
    # 1. Albedo of vegetation + bare soil (no snow)
    (; albedo_land, albedo_high_vegetation, albedo_low_vegetation) = albedo_scheme

    if haskey(diagn, :vegetation_high) && haskey(diagn, :vegetation_low)
        (; vegetation_high, vegetation_low) = diagn

        # linear combination of high and low vegetation and bare soil
        diagn.albedo[ij] = vegetation_high[ij] * albedo_high_vegetation +
                        vegetation_low[ij] * albedo_low_vegetation +
                        albedo_land * (1 - vegetation_high[ij] - vegetation_low[ij])
    else
        diagn.albedo[ij] = albedo_land
    end

    # 2. Add snow cover
    if haskey(land, :snow_depth)
        (; snow_depth) = land
        (; albedo_snow, snow_depth_scale) = albedo_scheme

        # how to compute snow cover from snow depth
        snow_cover_scheme = albedo_scheme.snow_cover

        # compute snow-cover fraction using the chosen scheme and clamp to [0, 1]
        snow_cover = snow_cover_scheme(snow_depth[ij], snow_depth_scale)

        # set land albedo linearly between bare land and snow depending on snow cover [0, 1]
        diagn.albedo[ij] += snow_cover * albedo_snow
    end
end
