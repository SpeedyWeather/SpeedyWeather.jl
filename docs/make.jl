using Documenter
# using SpeedyWeather

makedocs(
    format = Documenter.HTML(prettyurls=get(ENV, "CI", nothing)=="true",
                             ansicolor=true, collapselevel=1,
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
                "Analysis"=>"analysis.md",
                "Tree structure"=>"structure.md",
                "Particle advection"=>"particles.md",
                "NetCDF output"=>"output.md",
            ],
            "Extending SpeedyWeather" => [
                "Extensions"=>"extensions.md",
                "Forcing and drag"=>"forcing_drag.md",
                "Parameterizations"=>"parameterizations.md",
                "Orography"=>"orography.md",
                "Land-Sea Mask"=>"land_sea_mask.md",
                "Ocean"=>"ocean.md",
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
            "LowerTriangularMatrices"=>"lowertriangularmatrices.md",
            "SpeedyTransforms" => [
                "Spectral transforms" => "speedytransforms.md",
                "Gradient operators" => "gradients.md",
            ],
            "Function and type index"=>"functions.md",
            ]
)

deploydocs(
    repo = "github.com/SpeedyWeather/SpeedyWeather.jl.git",
    devbranch = "main",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"]
)
