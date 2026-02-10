const assets_url = "https://github.com/SpeedyWeather/SpeedyWeatherAssets/raw/refs/heads/main"

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
function get_asset(path::String; name::String = "", type = NCDataset, format = NCDataset)
    if isfile(path) # check if path is local (custom input)
        try
            asset = _get_asset(path, name, type, format)
            return asset
        catch e
            @warn "Custom asset loading failed with: $e"
            @warn "Attempting to load default asset"
        end
    end


    filename = path[end]
    url = joinpath(assets_url, path...)
    project_root = pkgdir(SpeedyWeather)
    artifact_toml = joinpath(project_root, "Artifacts.toml")

    # Check if the artifact already exists in the TOML
    hash = Artifacts.artifact_hash(filename, artifact_toml)

    # If not found or file not installed create and bind it
    if hash === nothing || !Artifacts.artifact_exists(hash)
        hash = PkgA.create_artifact() do artifact_dir
            Artifacts.download(url, joinpath(artifact_dir, filename))
        end

        # Also creates the Artifact.toml if it's missing
        PkgA.bind_artifact!(artifact_toml, filename, hash)
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
