using Documenter, SpeedyWeather

makedocs(
    format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename="SpeedyWeather.jl",
    authors="M Klöwer, T Kimpson, A White",
    modules=[SpeedyWeather],
    pages=["Home"=>"index.md",
            "How to run SpeedyWeather.jl"=>"how_to_run_speedy.md",
            "Spherical harmonic transform"=>"spectral_transform.md",
            "Dynamical core"=>"dynamical_core.md",
            "Parameterizations"=>"parametrizations.md",
            "Boundary conditions"=>"boundary_conditions.md",
            "New model setups"=>"new_model_setups.md",
            "Function and type index"=>"functions.md",
            "Style and convention guide"=>"conventions.md",
            "Development notes"=>"develop.md"]
)

deploydocs(
    repo = "github.com/milankl/SpeedyWeather.jl.git",
    devbranch = "main",
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"]
)