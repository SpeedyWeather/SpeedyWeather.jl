using Documenter, SpeedyWeather

makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename="SpeedyWeather.jl",
    authors="M Klöwer",
    modules=[SpeedyWeather],
    pages=["Home"=>"index.md"]
)

deploydocs(
    repo = "github.com/milankl/SpeedyWeather.jl.git",
    devbranch = "main"
)