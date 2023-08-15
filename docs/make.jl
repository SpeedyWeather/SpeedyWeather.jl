using Documenter, SpeedyWeather

makedocs(
    format=Documenter.HTML(prettyurls=get(ENV, "CI", nothing)=="true",
    ansicolor=true),
    sitename="SpeedyWeather.jl",
    authors="M Klöwer and SpeedyWeather contributors",
    modules=[SpeedyWeather],
    pages=["Home"=>"index.md",
            "Installation"=>"installation.md",
            "How to run SpeedyWeather.jl"=>"how_to_run_speedy.md",
            "Spherical Harmonic Transform"=>"spectral_transform.md",
            "Grids"=>"grids.md",
            "Barotropic model"=>"barotropic.md",
            "Shallow water model"=>"shallowwater.md",
            "Primitive equation model"=>"primitiveequation.md",
            "Parameterizations"=>"parameterizations.md",
            "Extending SpeedyWeather"=>"extending.md",
            "NetCDF output"=>"output.md",
            "Submodule: RingGrids"=>"ringgrids.md",
            "Submodule: LowerTriangularMatrices"=>"lowertriangularmatrices.md",
            "Submodule: SpeedyTransforms"=>"speedytransforms.md",
            "Style and convention guide"=>"conventions.md",
            "Development notes"=>"development.md",
            "Function and type index"=>"functions.md"]
)

deploydocs(
    repo = "github.com/SpeedyWeather/SpeedyWeather.jl.git",
    devbranch = "main",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"]
)
