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
    land = AlbedoClimatology(SG))
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

# For GPU usage just discard the extra information and treat it as a `ManualAlbedo`
Adapt.adapt_structure(to, albedo::AlbedoClimatology) = adapt(to, ManualAlbedo(albedo.albedo))

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