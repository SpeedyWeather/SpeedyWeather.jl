"""
Abstract super type for land-sea masks. Custom land-sea masks have to be defined as

    CustomMask{NF<:AbstractFloat, Grid<:AbstractGrid{NF}} <: AbstractLandSeaMask

and need to have at least a field called `mask::Grid` that uses a `Grid` as defined
by the spectral grid object, so of correct size and with the number format `NF`.
All `AbstractLandSeaMask` have a convenient generator function to be used like
`mask = CustomMask(spectral_grid, option=argument)`, but you may add your own or customize by
defining `CustomMask(args...)` which should return an instance of type `CustomMask{NF, Grid}` with
parameters matching the spectral grid. Then the initialize function has to be extended for
that new mask

    initialize!(mask::CustomMask, model::PrimitiveEquation)

which generally is used to tweak the mask.mask grid as you like, using
any other options you have included in `CustomMask` as fields or anything else (preferrably read-only,
because this is only to initialize the land-sea mask, nothing else) from `model`. You can
for example read something from file, set some values manually, or use coordinates from `model.geometry`.

The land-sea mask grid is expected to have values between [0, 1] as we use a fractional mask,
allowing for grid points being, e.g. quarter land and three quarters sea for 0.25
with 0 (=sea) and 1 (=land). The surface fluxes will weight proportionally the fluxes e.g.
from sea and land surface temperatures. Note however, that the land-sea mask can declare
grid points being (at least partially) ocean even though the sea surface temperatures
aren't defined (=NaN) in that grid point. In that case, not flux is applied."""
abstract type AbstractLandSeaMask end

function Base.show(io::IO, L::AbstractLandSeaMask)
    println(io, "$(typeof(L)) <: AbstractLandSeaMask")
    keys = propertynames(L)
    print_fields(io, L, keys)
end

function mask!(
    field::AbstractField,
    mask::AbstractField,
    land_or_sea::Symbol;
    masked_value = NaN,
)

    val = land_or_sea == :land ? 1 : 0
    masked_val = convert(eltype(field), masked_value)

    @boundscheck fields_match(field, mask, horizontal_only=true) || throw(DimensionMismatch(field, mask))
    @boundscheck ndims(mask) == 1 || throw(DimensionMismatch(field, mask))

    for k in eachlayer(field)
        for ij in eachgridpoint(mask)
            if mask[ij] == val
                field[ij, k] = masked_val
            end
        end
    end

    return field
end

# also allow for land_sea_mask struct to be passed on, use .mask in that case
mask!(field::AbstractField, mask::AbstractLandSeaMask, args...; kwargs...) =
    mask!(field, mask.mask, args...; kwargs...)

# make available when using SpeedyWeather
export EarthLandSeaMask

"""Land-sea mask, fractional, read from file.
$(TYPEDFIELDS)"""
@kwdef struct EarthLandSeaMask{NF, GridVariable2D} <: AbstractLandSeaMask

    # OPTIONS
    "path to the folder containing the land-sea mask file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "filename of land sea mask"
    file::String = "land-sea_mask.nc"

    "Grid the land-sea mask file comes on"
    file_Grid::Type{<:AbstractGrid} = FullClenshawGrid

    # FIELDS (to be initialized in initialize!)
    "Land-sea mask [1] on grid-point space. Land=1, sea=0, land-area fraction in between."
    mask::GridVariable2D
end

export LandSeaMask
const LandSeaMask = EarthLandSeaMask

"""
$(TYPEDSIGNATURES)
Generator function pulling the resolution information from `spectral_grid`."""
function (L::Type{<:AbstractLandSeaMask})(spectral_grid::SpectralGrid; kwargs...)
    (; NF, GridVariable2D, nlat_half) = spectral_grid
    mask = zeros(GridVariable2D, nlat_half)
    return L{NF, GridVariable2D}(; mask, kwargs...)
end

# set mask with grid, scalar, function; just define path `mask.mask` to grid here
function set!(mask::AbstractLandSeaMask, args...; kwargs...)
    set!(mask.mask, args...; kwargs...)
    lo, hi = extrema(mask.mask)
    (lo < 0 || hi > 1) && @warn "Land-sea mask was not set to values in [0, 1] but in [$lo, $hi]. Clamping."
    clamp!(mask.mask, 0, 1)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Reads a high-resolution land-sea mask from file and interpolates (grid-cell average)
onto the model grid for a fractional sea mask."""
function initialize!(land_sea_mask::EarthLandSeaMask, model::PrimitiveEquation)

    (; file_Grid) = land_sea_mask

    # LOAD NETCDF FILE
    if land_sea_mask.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__, "../../input_data", land_sea_mask.file)
    else
        path = joinpath(land_sea_mask.path, land_sea_mask.file)
    end
    ncfile = NCDataset(path)
    
    # high resolution land-sea mask
    # F = RingGrids.field_type(file_Grid)   # TODO this isn't working, hardcode instead
    lsm_highres = FullClenshawField(ncfile["lsm"].var[:, :], input_as=Matrix)

    # average onto grid cells of the model
    RingGrids.grid_cell_average!(land_sea_mask.mask, lsm_highres)

    # TODO this shoudln't be necessary, but at the moment grid_cell_average! can return values > 1
    # lo, hi = extrema(land_sea_mask.mask)
    # (lo < 0 || hi > 1) && @warn "Land-sea mask has values in [$lo, $hi], clamping to [0, 1]."
    clamp!(land_sea_mask.mask, 0, 1)
end

export AquaPlanetMask

"""Land-sea mask with zero = sea everywhere.
$(TYPEDFIELDS)"""
@kwdef struct AquaPlanetMask{NF, GridVariable2D} <: AbstractLandSeaMask
    "Land-sea mask [1] on grid-point space. Land=1, sea=0, land-area fraction in between."
    mask::GridVariable2D
end

"""
$(TYPEDSIGNATURES)
Sets all grid points to 0 = sea."""
function initialize!(land_sea_mask::AquaPlanetMask, model::PrimitiveEquation)
    land_sea_mask.mask .= 0    # set all to sea
    return nothing
end

export RockyPlanetMask

"""Land-sea mask with one = land everywhere.
$(TYPEDFIELDS)"""
@kwdef struct RockyPlanetMask{NF, GridVariable2D} <: AbstractLandSeaMask
    "Land-sea mask [1] on grid-point space. Land=1, sea=0, land-area fraction in between."
    mask::GridVariable2D
end

"""
$(TYPEDSIGNATURES)
Sets all grid points to 1 = land."""
function initialize!(land_sea_mask::RockyPlanetMask, model::PrimitiveEquation)
    land_sea_mask.mask .= 1    # set all to land
    return nothing
end