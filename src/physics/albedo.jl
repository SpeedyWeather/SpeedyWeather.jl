abstract type AbstractAlbedo <: AbstractModelComponent end

export Albedo
@kwdef struct Albedo{Ocean, Land} <: AbstractAlbedo
    ocean::Ocean
    land::Land
end

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
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    albedo::AlbedoClimatology,
    model::PrimitiveEquation,
)
    albedo!(diagn.physics.ocean, progn, albedo, model)
    albedo!(diagn.physics.land, progn, albedo, model)
end

function albedo!(
    diagn::AbstractDiagnosticVariables,     # ocean or land!
    progn::PrognosticVariables,
    albedo::AlbedoClimatology,
    model::PrimitiveEquation,
)
    diagn.albedo .= albedo.albedo           # copy over the albedo field
end