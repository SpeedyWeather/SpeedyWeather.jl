using Documenter
using DocumenterVitepress
using SpeedyWeatherInternals, LowerTriangularArrays, RingGrids, SpeedyTransforms, SpeedyWeather

makedocs(
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/SpeedyWeather/SpeedyWeather.jl",
        devbranch = "main",
        devurl = "dev",
        deploy_url = "./SpeedyWeatherDocumentation/" # adjust if cross-repo deployment to ./repo_name/
    ),
    sitename = "SpeedyWeather.jl",
    authors = "M Klöwer and SpeedyWeather contributors",
    modules = [SpeedyWeather, SpeedyWeatherInternals, LowerTriangularArrays, RingGrids, SpeedyTransforms],
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Usage" => [
            "Installation" => "installation.md",
            "How to run SpeedyWeather" => "how_to_run_speedy.md",
            "Examples 2D" => "examples_2D.md",
            "Examples 3D" => "examples_3D.md",
            "Initial conditions" => "initial_conditions.md",
            "Time integration" => "time_integration.md",
            "Tracer advection" => "tracers.md",
            "Particle advection" => "particles.md",
            "Stochastic physics" => "stochastic_physics.md",
            "Ocean models" => "ocean.md",
            "Sea ice models" => "sea_ice.md",
            "Land surface models" => "land.md",
            "Analysis" => "analysis.md",
            "Variables" => "variables.md",
            "Models" => "models.md",
            "Differentiability and Adjoint Model" => "differentiability.md",
            "NetCDF output" => "output.md",
            "Other output" => "other_output.md",
            "GPU and Architectures" => "architectures_gpu.md",
        ],
        "Advanced" => [
            "Extensions" => "extensions.md",
            "Variable system" => "variable_system.md",
            "Forcing and drag" => "forcing_drag.md",
            "Parameterizations" => "parameterizations.md",
            "Input data" => "input_data.md",
            "Orography" => "orography.md",
            "Land-Sea Mask" => "land_sea_mask.md",
            "Ocean" => "custom_ocean.md",
            "NetCDF output variables" => "custom_netcdf_output.md",
            "Callbacks" => "callbacks.md",
        ],
        "Physics" => [
            "Barotropic model" => "barotropic.md",
            "Shallow water model" => "shallowwater.md",
            "Primitive equation model" => "primitiveequation.md",
            "Large-scale condensation" => "large_scale_condensation.md",
            "Convection" => "convection.md",
            "Radiation" => "radiation.md",
            "Vertical diffusion" => "vertical_diffusion.md",
            "Surface fluxes" => "surface_fluxes.md",
        ],
        "Numerics" => [
            "Grids" => "grids.md",
            "RingGrids" => "ringgrids.md",
            "LowerTriangularArrays" => "lowertriangularmatrices.md",
            "Spherical Harmonic Transform" => "spectral_transform.md",
            "SpeedyTransforms" => "speedytransforms.md",
            "Gradient operators" => "gradients.md",
        ],
        "API" => "functions.md",
    ]
)

"""
    recursive_find(directory, pattern)

Return list of filepaths within `directory` that contains the `pattern::Regex`.
"""
function recursive_find(directory, pattern)
    return mapreduce(vcat, walkdir(directory)) do (root, dirs, filenames)
        matched_filenames = filter(contains(pattern), filenames)
        map(filename -> joinpath(root, filename), matched_filenames)
    end
end

# remove all .jld2 and .nc files in the docs folder from simulations
for pattern in [r"\.jld2", r"\.nc"]
    filenames = recursive_find(@__DIR__, pattern)

    for filename in filenames
        rm(filename)
    end
end

DocumenterVitepress.deploydocs(
    repo = "github.com/SpeedyWeather/SpeedyWeather.jl",
    deploy_repo = "github.com/SpeedyWeather/SpeedyWeatherDocumentation",
    devbranch = "main",
    push_preview = true,
)
