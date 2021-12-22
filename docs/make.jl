using Documenter, SpeedyWeather

makedocs(
    format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename="SpeedyWeather.jl",
    authors="M Klöwer",
    modules=[SpeedyWeather],
    pages=["Home"=>"index.md",
            "How to run SpeedyWeather.jl"=>"how_to_run_speedy.md",
            "Dynamical core"=>"dynamical_core.md",
            "Time integration"=>"time_integration.md",
            "Parameterizations"=>"parameterizations.md",
            "Boundary conditions"=>"boundary_conditions.md",
            "New model setups"=>"new_model_setups.md",
            "Function and type index"=>"functions.md"]
)

deploydocs(
    repo = "github.com/milankl/SpeedyWeather.jl.git",
    devbranch = "main"
)