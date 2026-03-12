#= These 3 const are conceptually the Project.toml for assets.

1. When releasing a new version of SpeedyWeatherAssets extend the AVAILABLE_ASSETS_VERSIONS tuple, e.g. (v"1", v"1.1").
Do this manually as we can't automatically detect new versions from the repo.
This allows a user to use a newer version even if not the default yet.

2. When a new version of SpeedyWeatherAssets is supposed to be the new DEFAULT_ASSETS_VERSION change it here accordingly.
SpeedyWeatherAssets should use semantic versioning, so v"1.1" should only add new assets. "v2" would be a breaking change.
=#
const ASSETS_URL = "https://github.com/SpeedyWeather/SpeedyWeatherAssets/raw/refs"
const DEFAULT_ASSETS_VERSION = v"1"
const AVAILABLE_ASSETS_VERSIONS = (v"1", )

# end SpeedyWeatherAssets "Project.toml"

"""$(TYPEDSIGNATURES)
Downloads a file under `path` from the SpeedyWeatherAssets repo (if `from_assets=true`,
otherwise locally), adds it to Artifacts.toml in the project root, and returns the file path.
Use `version = v"1+branch"` to retrieve data from a branch,
e.g. `version = v"1+main"` for the latest on the main branch.
`name` is used to find a variable inside the file, `ArrayType` to determine how to wrap
the loaded data into (in most cases) a `Field`, e.g. `FullGaussianField`, therefore determininig
the grid the data comes on. `FileFormat` determines how to read the file, e.g. `NCDataset` for NetCDF files.
"""
function get_asset(
        path::String;
        from_assets::Bool = true,
        name::String = "",
        ArrayType = Array,
        FileFormat = NCDataset,
        version = DEFAULT_ASSETS_VERSION
    )

    if !from_assets     # try to load locally
        if isfile(path) # check if path is local (custom input)
            try
                return _get_asset(path, name, ArrayType, FileFormat)
            catch e
                throw("Local asset loading failed with: $e")
            end
        else
            throw("Local $name file doesn't exist on this path")
        end
    end

    if version isa String
        # interpret as branch name
        target_url = joinpath(ASSETS_URL, "heads/$version")
    else
        if version in AVAILABLE_ASSETS_VERSIONS
            target_url = joinpath(ASSETS_URL, "tags/v$version")
        else
            available_str = join(AVAILABLE_ASSETS_VERSIONS, ", ")
            msg = "SpeedyWeatherAssets version $version not available. " *
                "Please select from: $available_str, or use `version = \"branch\"` to specify a branch."
            throw(ArgumentError(msg))
        end
    end

    filename = basename(path)
    url = joinpath(target_url, path)
    project_root = pkgdir(SpeedyWeather)
    # fallback to @__DIR__ if pkgdir fails (e.g. in dev mode)
    project_root = isnothing(project_root) ? (@__DIR__) : project_root
    artifact_toml = joinpath(project_root, "Artifacts.toml")

    # Check if the artifact already exists in the TOML
    hash = Artifacts.artifact_hash(filename, artifact_toml)

    # If not found or file not installed create and bind it
    # TODO remove Pkg.Artifacts dependency and use Artifacts.jl directly (just a question of interface presumably) ?
    # but we currently prefer bypassing the artifact"abc" logic
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

    return _get_asset(asset_path, name, ArrayType, FileFormat)
end


"""$(TYPEDSIGNATURES)
Search for a variable name `name` in `ncfile`. If `name` is empty, try to find a single variable that is not a dimension.
If multiple candidates are found, throw an error."""
function get_nc_variable_name(ncfile::NCDataset, name::String)

    if !haskey(ncfile, name) && name != ""
        # Helper to show user available var names
        available = join(keys(ncfile), ", ")
        error("Variable '$name' not found in file. Available variables: [$available]")

    elseif name == ""
        candidates = filter(k -> k ∉ keys(ncfile.dim), keys(ncfile))

        if isempty(candidates)
            error("No suitable variable found in asset (all variables matched dimension names)")
        elseif length(candidates) > 1
            msg = join(candidates, ", ")
            error("Ambiguous asset: found $(length(candidates)) variables ($msg). Please specify a `varname` kwarg")
        end

        target_name = candidates[1]
        @warn "No asset variable name provided; using the only available candidate: $target_name"
        return target_name
    else
        return name
    end
end

# lazy oad from NCDataset into netCDF variable (actually CommonDataModel.CFVariable)
function _get_asset(path::String, name::String, ArrayType::Type{<:NCDataset}, FileFormat::Type{<:NCDataset})
    ds = NCDataset(path)
    # TODO also read lat, lon from file and flip array in case it's not as expected
    target_name = get_nc_variable_name(ds, name)
    return ds[target_name], ds
end

# load from NetCDF into Array
function _get_asset(path::String, name::String, ArrayType::Type{<:Array}, FileFormat::Type{<:NCDataset})
    v, ds = _get_asset(path, name, NCDataset, FileFormat)
    data = load_shape_preserving(v, Val(ndims(v)))
    if eltype(data) <: AbstractFloat    # exclude case of loading integer data, e.g. land-sea mask
        fill_value = get(v.attrib, "_FillValue", NaN)
        # use === to include NaNs but also excpect fill value to have same type as data
        data[data .=== fill_value] .= NaN
    end
    close(ds)
    return data
end

# TODO generalise via some get method?
# TODO the .var drops the Union{Missing, T} into simply T, change that?
load_shape_preserving(ncvar, N::Val{1}) = ncvar.var[:]
load_shape_preserving(ncvar, N::Val{2}) = ncvar.var[:, :]
load_shape_preserving(ncvar, N::Val{3}) = ncvar.var[:, :, :]
load_shape_preserving(ncvar, N::Val{4}) = ncvar.var[:, :, :, :]
load_shape_preserving(ncvar, N::Val{5}) = ncvar.var[:, :, :, :, :]
# add more if needed

function _get_asset(path::String, name::String, ArrayType::Type{<:RingGrids.AbstractFullField}, FileFormat::Type{<:NCDataset})
    data = _get_asset(path, name, Array, FileFormat)        # first load as Array, replace all fill values with NaN
    return ArrayType(data, input_as = Matrix)
end