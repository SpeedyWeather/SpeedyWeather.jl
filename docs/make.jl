using Documenter, SpeedyWeather

makedocs(
    format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename="SpeedyWeather.jl",
    authors="M Klöwer, T Kimpson",
    modules=[SpeedyWeather],
    pages=["Home"=>"index.md",
            "How to run SpeedyWeather.jl"=>"how_to_run_speedy.md",
            "Spherical harmonic transform"=>"spectral_transform.md",
            "Dynamical core"=>"dynamical_core.md",
            "Parameterizations"=>"parameterizations.md",
            "Boundary conditions"=>"boundary_conditions.md",
            "New model setups"=>"new_model_setups.md",
            "Function and type index"=>"functions.md",
            "Style and convention guide"=>"conventions.md"]
)

deploydocs(
    repo = "github.com/milankl/SpeedyWeather.jl.git",
    devbranch = "main",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#", devurl => devurl]
)