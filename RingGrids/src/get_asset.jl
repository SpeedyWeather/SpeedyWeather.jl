#= These 2 consts are conceptually the Project.toml for assets.

1. When a new version of SpeedyWeatherAssets is supposed to be the new DEFAULT_ASSETS_VERSION
change it here accordingly. SpeedyWeatherAssets should use semantic versioning, so v"1.1" should
only add new assets. v"2" would be a breaking change.

2. Use `version = "branch"` (a String) to retrieve data from a branch,
e.g. `version = "main"` for the latest on the main branch.
=#
const ASSETS_URL = "https://github.com/SpeedyWeather/SpeedyWeatherAssets/raw/refs"
const DEFAULT_ASSETS_VERSION = v"1"
const AVAILABLE_ASSETS_VERSIONS = (v"1",)

# end SpeedyWeatherAssets "Project.toml"

"""$(TYPEDSIGNATURES)
Downloads a file under `path` from an assets repository (if `from_assets=true`,
otherwise loads locally), adds it to Artifacts.toml in the package directory, and returns
the loaded data. Use `version = "branch"` (a String) to retrieve data from a branch,
e.g. `version = "main"` for the latest on the main branch.
`name` is used to find a variable inside the file, `ArrayType` to determine how to wrap
the loaded data (in most cases into a `Field` subtype, determining the grid).
`FileFormat` determines how to read the file; load NCDatasets.jl to enable reading NetCDF files.
`assets_url` can be set to a custom URL for alternative asset repositories (e.g. for Terrarium).
`version` is the version tag (a `VersionNumber`) or branch name (a `String`) to download from.
`fill_value` is the value used to replace missing data (identified by the file's `_FillValue`
attribute) in the loaded array; defaults to `NaN`.
`output_grid` is the grid to interpolate the loaded data to; defaults to `nothing` (no interpolation). `ArrayType` must be a `AbstractFullField` type for interpolation to work.
"""
function get_asset(
        path::String;
        from_assets::Bool = true,
        name::String = "",
        ArrayType = Array,
        FileFormat = nothing,
        assets_url::String = ASSETS_URL,
        version = DEFAULT_ASSETS_VERSION,
        fill_value = NaN,
        output_grid::Union{Nothing, AbstractGrid} = nothing,
    )

    if !from_assets     # try to load locally
        if isfile(path) # check if path is local (custom input)
            try
                return load_asset(path, name, ArrayType, FileFormat, fill_value)
            catch e
                throw("Local asset loading failed with: $e")
            end
        else
            throw("Local $name file doesn't exist on this path")
        end
    end

    if version isa String
        # interpret as branch name
        target_url = joinpath(assets_url, "heads/$version")
    else
        # assume VersionNumber: construct tag URL and try to download
        target_url = joinpath(assets_url, "tags/v$version")
    end

    filename = basename(path)
    url = joinpath(target_url, path)
    project_root = pkgdir(RingGrids)
    # fallback to @__DIR__ if pkgdir fails
    project_root = isnothing(project_root) ? (@__DIR__) : project_root
    artifact_toml = joinpath(project_root, "Artifacts.toml")

    # Check if the artifact already exists in the TOML
    hash = Artifacts.artifact_hash(filename, artifact_toml)

    # If not found or file not installed create and bind it
    if hash === nothing || !Artifacts.artifact_exists(hash)
        hash = Pkg.Artifacts.create_artifact() do artifact_dir
            dest_path = joinpath(artifact_dir, filename)
            try
                Artifacts.download(url, dest_path)
            catch e
                @error "Download failed for URL: $url" exception = (e)
                throw("Could not download asset '$filename'. Check your internet connection or if the URL asset path exists.")
            end
        end
        Pkg.Artifacts.bind_artifact!(artifact_toml, filename, hash, force = true)
    end

    asset_path = joinpath(Artifacts.artifact_path(hash), filename)

    data = load_asset(asset_path, name, ArrayType, FileFormat, fill_value)

    if !isnothing(output_grid)
        @assert ArrayType <: AbstractFullField "output_grid can only be used with AbstractFullField types as ArrayType, got $ArrayType"
        # interpolate to output_grid
        out_field = zeros(eltype(data), output_grid, size(data)[2:end]...)
        interpolate!(out_field, data)
        return out_field
    end

    return data
end

# load from array into a RingGrid Field
function load_asset(path::String, name::String, ArrayType::Type{<:AbstractFullField}, FileFormat, fill_value)
    data = load_asset(path, name, Array, FileFormat, fill_value)        # first load as Array
    return ArrayType(data, input_as = Matrix)
end

# TODO generalise via some get method?
# TODO the .var drops the Union{Missing, T} into simply T, change that?
load_shape_preserving(ncvar, N::Val{1}) = ncvar.var[:]
load_shape_preserving(ncvar, N::Val{2}) = ncvar.var[:, :]
load_shape_preserving(ncvar, N::Val{3}) = ncvar.var[:, :, :]
load_shape_preserving(ncvar, N::Val{4}) = ncvar.var[:, :, :, :]
load_shape_preserving(ncvar, N::Val{5}) = ncvar.var[:, :, :, :, :]
# add more if needed
