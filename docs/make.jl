using Documenter, SpeedyWeather

makedocs(
    format = Documenter.HTML(prettyurls=get(ENV, "CI", nothing)=="true",
                             ansicolor=true, collapselevel=1,
                             size_threshold = 600_000),      # in bytes
    sitename = "SpeedyWeather.jl",
    authors = "M Klöwer and SpeedyWeather contributors",
    modules = [SpeedyWeather],
    pages = ["Home"=>"index.md",
            "Installation"=>"installation.md",
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
            "Running SpeedyWeather" => [
                "How to run SpeedyWeather"=>"how_to_run_speedy.md",
                "Model setups"=>"setups.md",
                "Analysis"=>"analysis.md",
                "Tree structure"=>"structure.md",
                "Particle advection"=>"particles.md",
                "NetCDF output"=>"output.md",
            ],
            "Extending SpeedyWeather" => [
                "Forcing and drag"=>"forcing_drag.md",
                "Parameterizations"=>"parameterizations.md",
                "Orography"=>"orography.md",
                "Land-Sea Mask"=>"land_sea_mask.md",
                "Ocean"=>"ocean.md",
                "Callbacks"=>"callbacks.md",
            ],
            "RingGrids"=>"ringgrids.md",
            "LowerTriangularMatrices"=>"lowertriangularmatrices.md",
            "SpeedyTransforms"=>"speedytransforms.md",
            "Function and type index"=>"functions.md",
            ]
)

deploydocs(
    repo = "github.com/SpeedyWeather/SpeedyWeather.jl.git",
    devbranch = "main",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"]
)
