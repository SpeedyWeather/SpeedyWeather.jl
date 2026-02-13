const ASSETS_URL = "https://github.com/SpeedyWeather/SpeedyWeatherAssets/raw/refs"
const DEFAULT_ASSETS_VERSION = v"1"
const ASSETS_VERSIONS_AVAILABLE = (v"1",)

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

"""
Downloads a file from the SpeedyWeatherAssets repo, adds it to 
Artifacts.toml in the project root, and returns the file path.
"""
function get_asset(
        path::String;
        from_assets::Bool = true,
        name::String = "",
        type::Any = NCDataset,
        format::Any = NCDataset,
        version::VersionNumber = DEFAULT_ASSETS_VERSION
    )

    if !from_assets
        if isfile(path) # check if path is local (custom input)
            try
                return _get_asset(path, name, type, format)
            catch e
                throw("Local asset loading failed with: $e")
            end
        else
            throw("Local $name file doesn't exist on this path")
        end
    end

    if !isempty(version.build)
        branch = version.build[1]
        target_url = joinpath(ASSETS_URL, "heads/$branch")
    else
        if version in ASSETS_VERSIONS_AVAILABLE
            target_url = joinpath(ASSETS_URL, "tags/$version")
        else
            available_str = join(ASSETS_VERSIONS_AVAILABLE, ", ")
            msg = "SpeedyWeatherAssets version $version not available. " *
                "Please select from: $available_str, or the 'main' build"
            throw(ArgumentError(msg))
        end
    end

    filename = basename(path)
    url = joinpath(target_url, path)
    project_root = pkgdir(SpeedyWeather)
    artifact_toml = joinpath(project_root, "Artifacts.toml")

    # Check if the artifact already exists in the TOML
    hash = Artifacts.artifact_hash(filename, artifact_toml)

    # If not found or file not installed create and bind it
    if hash === nothing || !Artifacts.artifact_exists(hash)
        hash = PkgA.create_artifact() do artifact_dir
            dest_path = joinpath(artifact_dir, filename)
            try
                Artifacts.download(url, dest_path)
            catch e
                @error "Download failed for URL: $url" exception = (e)
                throw("Could not download asset '$filename'. Check your internet connection or if the URL asset path exists.")
            end
        end
        PkgA.bind_artifact!(artifact_toml, filename, hash, force = true)
    end

    asset_path = joinpath(Artifacts.artifact_path(hash), filename)

    return _get_asset(asset_path, name, type, format)
end

function _get_asset(path::String, name::String, type::Type{NCDataset}, format::Type{NCDataset})
    return NCDataset(path)
end

function _get_asset(path::String, name::String, type::Type{FullGaussianField}, format::Type{NCDataset})
    ncfile = NCDataset(path)
    target_name = get_nc_variable_name(ncfile, name)
    data = ncfile[target_name].var[:, :]
    close(ncfile)
    return type(data, input_as = Matrix)
end

function _get_asset(path::String, name::String, type::Type{<:AbstractGrid}, format::Type{NCDataset})
    ncfile = NCDataset(path)
    target_name = get_nc_variable_name(ncfile, name)
    v = ncfile[target_name]

    # Handle files with different dims
    if length(size(v)) == 2
        data = v.var[:, :]
    elseif length(size(v)) == 3
        data = v.var[:, :, :]
    end
    fill_value = get(v.attrib, "_FillValue", NaN)
    close(ncfile)
    return type(data, input_as = Matrix), fill_value
end
