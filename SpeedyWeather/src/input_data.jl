"""Default path for SpeedyWeather input data files."""
const INPUT_DATA_PATH = joinpath(@__DIR__, "../input_data")

"""
    input_data_path(path, file) -> String

Resolve the full path to an input data file. If `path` is the default
`"SpeedyWeather.jl/input_data"`, resolve relative to the package's input_data folder.
Otherwise, join `path` and `file` directly."""
function input_data_path(path::String, file::String)
    if path == "SpeedyWeather.jl/input_data"
        return joinpath(INPUT_DATA_PATH, file)
    else
        return joinpath(path, file)
    end
end

"""
$(TYPEDSIGNATURES)

Load a variable from a NetCDF file, wrap it into a grid, optionally replace fill values,
transfer to the given architecture, and interpolate onto the output field `A`.
Works for both 2D (lon, lat) and 3D (lon, lat, time) NetCDF variables.

# Arguments
- `A::AbstractField`: Output field to interpolate into.
- `path::String`: Folder path (use `"SpeedyWeather.jl/input_data"` for package default).
- `file::String`: NetCDF filename.
- `varname::String`: Variable name in the NetCDF file.

# Keyword arguments
- `file_Grid`: Grid type the NetCDF data is on (default `FullGaussianGrid`).
- `missing_value`: Value to replace `_FillValue` entries with (default `NaN`).
- `scale`: Multiplicative scaling applied after interpolation (default `1`).
- `architecture`: Target architecture (default: architecture of `A`).
"""
function load_from_netcdf!(
        A::AbstractField,
        path::String,
        file::String,
        varname::String;
        file_Grid::Type{<:AbstractGrid} = FullGaussianGrid,
        missing_value = NaN,
        scale = 1,
        architecture::AbstractArchitecture = RingGrids.architecture(A),
    )
    # resolve path and open file
    fullpath = input_data_path(path, file)
    ncfile = NCDataset(fullpath)

    # read variable: 2D (lon, lat) or 3D (lon, lat, time) and wrap into grid
    ndims_nc = length(size(ncfile[varname]))
    data_raw = ndims_nc == 2 ? ncfile[varname].var[:, :] : ncfile[varname].var[:, :, :]
    data_grid = file_Grid(data_raw, input_as = Matrix)

    # replace fill values with missing_value (=== to also match NaN fill values)
    if haskey(ncfile[varname].attrib, "_FillValue")
        fill_value = ncfile[varname].attrib["_FillValue"]
        data_grid[data_grid .=== fill_value] .= missing_value
    end

    # transfer to target architecture
    data_grid = on_architecture(architecture, data_grid)

    # interpolate onto output field
    interp = RingGrids.interpolator(A, data_grid, NF = Float32)
    interpolate!(A, data_grid, interp)

    # optional scaling
    scale != 1 && (A .*= scale)

    return A
end
