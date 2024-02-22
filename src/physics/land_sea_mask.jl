abstract type AbstractLandSeaMask{NF,Grid} end

# make available when using SpeedyWeather
export LandSeaMask, AquaPlanetMask

"""Land-sea mask, fractional, read from file.
$(TYPEDFIELDS)"""
Base.@kwdef struct LandSeaMask{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractLandSeaMask{NF,Grid}

    # OPTIONS
    "path to the folder containing the land-sea mask file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "filename of land sea mask"
    file::String = "land-sea_mask.nc"

    "Grid the land-sea mask file comes on"
    file_Grid::Type{<:AbstractGrid} = FullClenshawGrid

    # FIELDS (to be initialized in initialize!)
    "Land-sea mask [1] on grid-point space. Land=1, sea=0, land-area fraction in between."
    land_sea_mask::Grid
end

"""
$(TYPEDSIGNATURES)
Generator function pulling the resolution information from `spectral_grid`."""
function LandSeaMask(spectral_grid::SpectralGrid;kwargs...)
    (;NF, Grid, nlat_half) = spectral_grid
    land_sea_mask   = zeros(Grid{NF},nlat_half)
    return LandSeaMask{NF,Grid{NF}}(;land_sea_mask,kwargs...)
end

function initialize!(land_sea_mask::LandSeaMask)

    (;file_Grid) = land_sea_mask

    # LOAD NETCDF FILE
    if land_sea_mask.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__,"../../input_data",land_sea_mask.file)
    else
        path = joinpath(land_sea_mask.path,land_sea_mask.file)
    end
    ncfile = NCDataset(path)
    
    # high resolution land-sea mask
    lsm_highres = file_Grid(ncfile["lsm"][:,:])

    # average onto grid cells of the model
    RingGrids.grid_cell_average!(land_sea_mask.land_sea_mask,lsm_highres)

    # TODO this shoudln't be necessary, but at the moment grid_cell_average! can return values > 1
    clamp!(land_sea_mask.land_sea_mask,0,1)
end

"""Land-sea mask with zero = sea everywhere.
$(TYPEDFIELDS)"""
Base.@kwdef struct AquaPlanetMask{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractLandSeaMask{NF,Grid}
    "Land-sea mask [1] on grid-point space. Land=1, sea=0, land-area fraction in between."
    land_sea_mask::Grid
end

"""
$(TYPEDSIGNATURES)
Generator function pulling the resolution information from `spectral_grid`."""
function AquaPlanetMask(spectral_grid::SpectralGrid;kwargs...)
    (;NF, Grid, nlat_half) = spectral_grid
    land_sea_mask = zeros(Grid{NF},nlat_half)
    return AquaPlanetMask{NF,Grid{NF}}(;land_sea_mask,kwargs...)
end

function initialize!(land_sea_mask::AquaPlanetMask)
    land_sea_mask.land_sea_mask .= 0    # set all to sea
    return nothing
end