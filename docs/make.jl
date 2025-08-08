using Documenter
using SpeedyWeather

makedocs(
    format = Documenter.HTML(prettyurls=get(ENV, "CI", nothing)=="true",
                             ansicolor=true,
                             collapselevel=1,
                             canonical = "https://speedyweather.github.io/SpeedyWeatherDocumentation/stable/",
                             size_threshold = 600_000),      # in bytes
    sitename = "SpeedyWeather.jl",
    authors = "M Klöwer and SpeedyWeather contributors",
    modules = [SpeedyWeather],
    pages = ["Home"=>"index.md",
            "Installation"=>"installation.md",
            "Running SpeedyWeather" => [
                "How to run SpeedyWeather"=>"how_to_run_speedy.md",
                "Examples 2D"=>"examples_2D.md",
                "Examples 3D"=>"examples_3D.md",
                "Initial conditions" => "initial_conditions.md",
                "Tracer advection"=>"tracers.md",
                "Particle advection"=>"particles.md",
                "Stochastic physics" => "stochastic_physics.md",
                "Ocean and sea ice" => "ocean_sea_ice.md",
                "Land surface model" => "land.md",
                "Analysis"=>"analysis.md",
                "Tree structure"=>"structure.md",
                "Differentiability and Adjoint Model"=>"differentiability.md",
                "NetCDF output"=>"output.md",
                "GPU & Architectures" => "architectures_gpu.md",
            ],
            "Extending SpeedyWeather" => [
                "Extensions"=>"extensions.md",
                "Forcing and drag"=>"forcing_drag.md",
                "Parameterizations"=>"parameterizations.md",
                "Orography"=>"orography.md",
                "Land-Sea Mask"=>"land_sea_mask.md",
                "Ocean"=>"custom_ocean.md",
                "NetCDF output variables"=>"custom_netcdf_output.md",
                "Callbacks"=>"callbacks.md",
            ],
            "Dynamics" => [
                "Barotropic model"=>"barotropic.md",
                "Shallow water model"=>"shallowwater.md",
                "Primitive equation model"=>"primitiveequation.md",
            ],
            "Physics" => [
                "Large-scale condensation"=>"large_scale_condensation.md",
                "Convection"=>"convection.md",
                "Radiation"=>"radiation.md",
                "Vertical diffusion"=>"vertical_diffusion.md",
                "Surface fluxes"=>"surface_fluxes.md",
            ],
            "Discretization" => [
                "Spherical Harmonic Transform"=>"spectral_transform.md",
                "Grids"=>"grids.md",
            ],
            "RingGrids"=>"ringgrids.md",
            "LowerTriangularArrays"=>"lowertriangularmatrices.md",
            "SpeedyTransforms" => [
                "Spectral transforms" => "speedytransforms.md",
                "Gradient operators" => "gradients.md",
            ],
            "Function and type index"=>"functions.md",
            ]
)

"""
    recursive_find(directory, pattern)

Return list of filepaths within `directory` that contains the `pattern::Regex`.
"""
function recursive_find(directory, pattern)
    mapreduce(vcat, walkdir(directory)) do (root, dirs, filenames)
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

# Replace with below once https://github.com/JuliaDocs/Documenter.jl/pull/2692 is merged and available.
#  deploydocs(repo = "github.com/SpeedyWeather/SpeedyWeather.jl",
#    deploy_repo = "github.com/SpeedyWeather/SpeedyWeatherDocumentation",
#    devbranch = "main",
#    push_preview = true)
if get(ENV, "GITHUB_EVENT_NAME", "") == "pull_request"
    deploydocs(repo = "github.com/SpeedyWeather/SpeedyWeather.jl",
               repo_previews = "github.com/SpeedyWeather/SpeedyWeatherDocumentation",
               devbranch = "main",
               forcepush = true,
               push_preview = true,
               versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"])
else
    repo = "github.com/SpeedyWeather/SpeedyWeatherDocumentation"
    withenv("GITHUB_REPOSITORY" => repo) do
        deploydocs(repo = repo,
                   devbranch = "main",
                   forcepush = true,
                   versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"])
    end
end
