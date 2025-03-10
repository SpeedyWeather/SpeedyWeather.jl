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
    ocean = GlobalConstantAlbedo(SG, albedo=0.06),
    land = GlobalConstantAlbedo(SG, albedo=0.4))
    return Albedo(ocean, land)
end

function initialize!(albedo::Albedo, model::PrimitiveEquation)
    initialize!(albedo.ocean, model)
    initialize!(albedo.land, model)
end

# dispatch over model.albedo
albedo!(diagn::DiagnosticVariables, progn::PrognosticVariables, model::PrimitiveEquation) =
    albedo!(diagn, progn, model.albedo, model)

function albedo!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    albedo::Albedo,
    model::PrimitiveEquation,
)
    albedo!(diagn.physics.ocean, progn, albedo.ocean, model)
    albedo!(diagn.physics.land, progn, albedo.land, model)
end

## GLOBAL CONSTANT ALBEDO
export GlobalConstantAlbedo
@kwdef mutable struct GlobalConstantAlbedo{NF} <: AbstractAlbedo
    albedo::NF = 0.3
end

GlobalConstantAlbedo(SG::SpectralGrid; kwargs...) = GlobalConstantAlbedo{SG.NF}(; kwargs...)
initialize!(albedo::GlobalConstantAlbedo, ::PrimitiveEquation) = nothing
function albedo!(
    diagn::AbstractDiagnosticVariables,     # ocean or land!
    progn::PrognosticVariables,
    albedo::GlobalConstantAlbedo,
    model::PrimitiveEquation,
)
    diagn.albedo .= albedo.albedo           # copy over the constant albedo
end

## MANUAL ALBEDO
export ManualAlbedo

"""Manual albedo field, to be used with set! and is copied into the diagnostic variables on every time step.
Defined so that parameterizations can change the albedo at every time step (e.g. snow cover) without
losing the information of the original surface albedo. Fields are
$(TYPEDFIELDS)"""
struct ManualAlbedo{NF, Grid} <: AbstractAlbedo
    albedo::Grid
end

ManualAlbedo(SG::SpectralGrid) = ManualAlbedo{SG.NF, SG.GridVariable2D}(zeros(SG.GridVariable2D, SG.nlat_half))
initialize!(albedo::ManualAlbedo, model::PrimitiveEquation) = nothing
function albedo!(
    diagn::AbstractDiagnosticVariables,     # ocean or land!
    progn::PrognosticVariables,
    albedo::ManualAlbedo,
    model::PrimitiveEquation,
)
    diagn.albedo .= albedo.albedo           # copy over the manual albedo
end

## ALEBDO CLIMATOLOGY

export AlbedoClimatology
@kwdef struct AlbedoClimatology{NF, Grid} <: AbstractAlbedo
    "number of latitudes on one hemisphere, Equator included"
    nlat_half::Int

    "[OPTION] path to the folder containing the albedo file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "[OPTION] filename of albedo"
    file::String = "albedo.nc"

    "[OPTION] variable name in netcdf file"
    varname::String = "alb"

    "[OPTION] Grid the albedo file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "Albedo climatology"
    albedo::Grid = zeros(Grid, nlat_half)
end

function AlbedoClimatology(SG::SpectralGrid; kwargs...)
    (; NF, Grid, nlat_half) = SG
    AlbedoClimatology{NF,Grid{NF}}(;nlat_half, kwargs...)
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

function albedo!(
    diagn::AbstractDiagnosticVariables,     # ocean or land!
    progn::PrognosticVariables,
    albedo::AlbedoClimatology,
    model::PrimitiveEquation,
)
    diagn.albedo .= albedo.albedo           # copy over the albedo field
end