"""
Abstract type for land-sea masks. Custom land-sea masks have to be defined as

    CustomMask{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractLandSeaMask{NF,Grid}

and need to have at least a field called `land_sea_mask::Grid` that uses a `Grid` as defined
by the spectral grid object, so of correct size and with the number format NF.
It is therefore recommended (but not required) to write a generator function as follows

    function CustomMask(spectral_grid::SpectralGrid;kwargs...)
        (;NF, Grid, nlat_half) = spectral_grid
        land_sea_mask   = zeros(Grid{NF},nlat_half)
        return LandSeaMask{NF,Grid{NF}}(;land_sea_mask,kwargs...)
    end

to allow a convenient construction like `mask = CustomMask(spectral_grid,option=argument)`.
Then the initialize function has to be extended for that new mask

    initialize!(mask::CustomMask,model::PrimitiveEquation)

which generally is used to tweak the mask.land_sea_mask grid as you like, using
any other options you have included in `CustomMask` as fields
or anything else (preferrably read-only, because this is only to initialize the land-sea mask, nothing else)
from `model`.

The land-sea mask grid is expected to have values between [0,1] as we use a fractional mask,
allowing for grid points being, e.g. quarter land and three quarters sea for 0.25
with 0 (=sea) and 1 (=land). The surface fluxes will weight proportionally the fluxes e.g.
from sea and land surface temperatures. Note however, that the land-sea mask can declare
grid points being (at least partially) ocean even though the sea surface temperatures
aren't defined (=NaN) in that grid point. In that case, not flux is applied."""
abstract type AbstractLandSeaMask{NF,Grid} end

function Base.show(io::IO,L::AbstractLandSeaMask)
    println(io,"$(typeof(L)) <: AbstractLandSeaMask")
    keys = propertynames(L)
    print_fields(io,L,keys)
end

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

"""
$(TYPEDSIGNATURES)
Reads a high-resolution land-sea mask from file and interpolates (grid-call average)
onto the model grid for a fractional sea mask."""
function initialize!(land_sea_mask::LandSeaMask, model::PrimitiveEquation)

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

"""
$(TYPEDSIGNATURES)
Sets all grid points to sea."""
function initialize!(land_sea_mask::AquaPlanetMask, model::PrimitiveEquation)
    land_sea_mask.land_sea_mask .= 0    # set all to sea
    return nothing
end