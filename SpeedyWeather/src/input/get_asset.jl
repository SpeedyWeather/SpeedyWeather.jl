const assets_url = "https://github.com/SpeedyWeather/SpeedyWeatherAssets/raw/refs/heads/main"

"""
    get_asset(subfolder, filename)

Downloads a file from the SpeedyWeatherAssets repo, adds it to 
Artifacts.toml in the project root, and returns the file path.
"""
function get_asset(path::String...)
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

    return joinpath(Artifacts.artifact_path(hash), filename)
end
