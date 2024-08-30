abstract type AbstractAlbedo <: AbstractModelComponent end

## GLOBAL CONSTANT ALBEDO

export GlobalConstantAlbedo
Base.@kwdef struct GlobalConstantAlbedo{NF,Grid} <: AbstractAlbedo
    nlat_half::Int
    α::NF = 0.8
    albedo::Grid = zeros(Grid, nlat_half)
end

function GlobalConstantAlbedo(SG::SpectralGrid; kwargs...)
    (;NF, Grid, nlat_half) = SG
    GlobalConstantAlbedo{NF, Grid{NF}}(; nlat_half, kwargs...)
end

function initialize!(albedo::GlobalConstantAlbedo,::PrimitiveEquation)
    albedo.albedo .= albedo.α
end

## ALEBDO CLIMATOLOGY

export AlbedoClimatology
Base.@kwdef struct AlbedoClimatology{NF,Grid} <: AbstractAlbedo
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
