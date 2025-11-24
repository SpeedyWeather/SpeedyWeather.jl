abstract type AbstractAlbedo <: AbstractModelComponent end

export Albedo
@kwdef struct Albedo{Ocean, Land} <: AbstractAlbedo
    ocean::Ocean
    land::Land
end

function Base.show(io::IO, A::Albedo)
    println(io, "Albedo <: SpeedyWeather.AbstractAlbedo")
    properties = propertynames(A)
    n = length(properties)
    for (i, key) in enumerate(properties)
        val = getfield(A, key)
        s = i == n ? "└" : "├"  # choose ending └ for last property
        p = i == n ? print : println
        p(io, "$s $key: $(typeof(val))")
    end
end

export DefaultAlbedo
function DefaultAlbedo(SG::SpectralGrid;
    ocean = OceanSeaIceAlbedo(SG),
    land = AlbedoClimatology(SG))
    return Albedo(ocean, land)
end

function initialize!(albedo::Albedo, model::PrimitiveEquation)
    initialize!(albedo.ocean, model)
    initialize!(albedo.land, model)
end

# dispatch over model.albedo
albedo!(diagn::DiagnosticVariables, progn::PrognosticVariables, model::PrimitiveEquation) =
    albedo!(diagn, progn, model.albedo, model)

# composite albedos: call separately for ocean and land with .ocean and .land
function albedo!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    albedo::Albedo,
    model::PrimitiveEquation,
)
    albedo!(diagn.physics.ocean, diagn, progn, albedo.ocean, model)
    albedo!(diagn.physics.land, diagn, progn, albedo.land, model)
end

# single albedo: call separately for ocean and land with the same albedo
function albedo!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    albedo::AbstractAlbedo,
    model::PrimitiveEquation,
)
    albedo!(diagn.physics.ocean, diagn, progn, albedo, model)
    albedo!(diagn.physics.land, diagn, progn, albedo, model)
end

# single albedo, copy over the constant albedo
function albedo!(
    diagn::AbstractDiagnosticVariables,     # ocean or land!
    diagn_all::DiagnosticVariables,         # all diagnostics
    progn::PrognosticVariables,
    albedo::AbstractAlbedo,
    model::PrimitiveEquation,
)
    diagn.albedo .= albedo.albedo           # copy over the constant albedo
end

## GLOBAL CONSTANT ALBEDO
export GlobalConstantAlbedo
@kwdef mutable struct GlobalConstantAlbedo{NF} <: AbstractAlbedo
    albedo::NF = 0.3
end

GlobalConstantAlbedo(SG::SpectralGrid; kwargs...) = GlobalConstantAlbedo{SG.NF}(; kwargs...)
initialize!(albedo::GlobalConstantAlbedo, ::PrimitiveEquation) = nothing

## MANUAL ALBEDO
export ManualAlbedo

"""Manual albedo field, to be used with `set!` and is copied into the diagnostic variables on every time step.
Defined so that parameterizations can change the albedo at every time step (e.g. snow cover) without
losing the information of the original surface albedo. Fields are
$(TYPEDFIELDS)"""
struct ManualAlbedo{NF, GridVariable2D} <: AbstractAlbedo
    albedo::GridVariable2D
end

ManualAlbedo(SG::SpectralGrid) = ManualAlbedo{SG.NF, SG.GridVariable2D}(zeros(SG.GridVariable2D, SG.grid))
initialize!(albedo::ManualAlbedo, model::PrimitiveEquation) = nothing

## ALEBDO CLIMATOLOGY

export AlbedoClimatology
@kwdef struct AlbedoClimatology{NF, GridVariable2D} <: AbstractAlbedo
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
    (; NF, GridVariable2D, grid) = SG
    albedo = zeros(GridVariable2D, grid)
    AlbedoClimatology{NF, GridVariable2D}(; albedo, kwargs...)
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

    a = albedo.file_Grid(ncfile[albedo.varname].var[:, :], input_as=Matrix)
    interpolate!(albedo.albedo, a)
end

## OceanSeaIceAlbedo
export OceanSeaIceAlbedo

@kwdef struct OceanSeaIceAlbedo{NF} <: AbstractAlbedo
    albedo_ocean::NF = 0.06
    albedo_ice::NF = 0.6
end

OceanSeaIceAlbedo(SG::SpectralGrid; kwargs...) = OceanSeaIceAlbedo{SG.NF}(;kwargs...)
initialize!(albedo::OceanSeaIceAlbedo, model::PrimitiveEquation) = nothing

function albedo!(
    diagn::AbstractDiagnosticVariables,     # ocean only
    diagn_all::DiagnosticVariables,         # all diagnostics
    progn::PrognosticVariables,
    albedo::OceanSeaIceAlbedo,
    model::PrimitiveEquation,
)
    (; sea_ice_concentration ) = progn.ocean
    (; albedo_ocean, albedo_ice) = albedo

    # set ocean albedo linearly between ocean and ice depending on sea ice concentration
    diagn.albedo .= albedo_ocean .+ sea_ice_concentration .* (albedo_ice .- albedo_ocean)
end


export AbstractSnowAlbedoScheme, LinearSnowAlbedo, SaturatingSnowAlbedo
abstract type AbstractSnowAlbedoScheme end

"""Linear ramp: snow cover grows from 0 to 1 at `snow_depth_scale`."""
struct LinearSnowAlbedo <: AbstractSnowAlbedoScheme end

"""Saturating ramp: snow cover grows as `S/(S+scale)`."""
struct SaturatingSnowAlbedo  <: AbstractSnowAlbedoScheme end

"""$(TYPEDSIGNATURES) Snow cover fraction for the linear scheme, clamped to 1."""
@inline snow_cover_fraction(::LinearSnowAlbedo, S, scale) =
    min(S / scale, 1)

"""$(TYPEDSIGNATURES) Snow cover fraction for the saturating scheme."""
@inline snow_cover_fraction(::SaturatingSnowAlbedo, S, scale) =
    S / (S + scale)

## LandSnowAlbedo
export LandSnowAlbedo

@kwdef struct LandSnowAlbedo{NF, Scheme<:AbstractSnowAlbedoScheme} <: AbstractAlbedo
    "Albedo of bare land (excluding vegetation) [1]"
    albedo_land::NF = 0.3

    "Albedo of vegetation [1]"
    albedo_vegetation::NF = 0.15

    "Albedo of snow [1]"
    albedo_snow::NF = 0.8

    "Conversion from snow depth to snow cover [m]"
    snow_depth_scale::NF = 0.1

    "Snow cover-albedo scheme"
    scheme::Scheme = SaturatingSnowAlbedo()

end

function LandSnowAlbedo(SG::SpectralGrid; scheme = SaturatingSnowAlbedo(), kwargs...)
    return LandSnowAlbedo{SG.NF, typeof(scheme)}(; scheme, kwargs...)
end

initialize!(albedo::LandSnowAlbedo, model::PrimitiveEquation) = nothing

function albedo!(
    diagn::AbstractDiagnosticVariables,     # land only
    diagn_all::DiagnosticVariables,         # all diagnostics
    progn::PrognosticVariables,
    albedo::LandSnowAlbedo,
    model::PrimitiveEquation,
)
    (; snow_depth) = progn.land
    (; albedo_land, albedo_vegetation, albedo_snow, snow_depth_scale, scheme) = albedo

    snow_cover = diagn_all.dynamics.a_2D_grid       # scratch memory
    
    # compute snow-cover fraction using the chosen scheme and clamp to [0, 1]
    snow_cover .= snow_cover_fraction.(Ref(scheme), snow_depth, snow_depth_scale)
    @. snow_cover = min(max(snow_cover, 0), 1)
	# ## Now compute the snow area cover fraction based on snow depth
	# #σₛ = Sₐ / (10. * drag_snow + Sₐ) #JULES, from Betts et al.
	# σₛ[ij] = min(1.0, Sₐ[ij] / snow_depth_scale)   # e.g. snow_depth_scale ≈ 0.1 m
	# ## Now compute the change in albedo
	# albedo_snow[ij] = snow_albedo_old + (snow_albedo_fresh - snow_albedo_old) * exp(-Δt_snow / τˢ) # time needs to be expresssed in days. 
	# column.albedo_land[ij] = (1-σₛ[ij]) * albedo_land + σₛ[ij] * albedo_snow[ij]

    # set land albedo linearly between bare land and snow depending on snow cover [0, 1]
    diagn.albedo .= albedo_land .+ snow_cover .* (albedo_snow - albedo_land)
end
