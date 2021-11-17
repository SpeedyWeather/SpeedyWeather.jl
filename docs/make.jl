using Documenter, SpeedyWeather

makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename="SpeedyWeather.jl",
    authors="M Klöwer",
    modules=[BitInformation],
    pages=["Home"=>"index.md"]
)

deploydocs(
    repo = "github.com/milankl/SpeedyWeatehr.jl.git",
    devbranch = "main"
)