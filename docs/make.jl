using Documenter,
      Literate,
      CairoMakie,   # to avoid capturing precompilation output by Literate
      SpeedyWeather

CairoMakie.activate!(type = "svg")

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")

examples = [
    "basic_example.jl"
]

for example in examples
    example_filepath = joinpath(EXAMPLES_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor())
end

example_pages = [
    "Basic example" => "literated/basic_example.md"
]

makedocs(
    format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename="SpeedyWeather.jl",
    authors="M Klöwer and SpeedyWeather contributors",
    modules=[SpeedyWeather],
    pages=["Home"=>"index.md",
            "How to run SpeedyWeather.jl"=>"how_to_run_speedy.md",
            "Examples" => example_pages,
            "Spherical harmonic transform"=>"spectral_transform.md",
            "Grids"=>"grids.md",
            "Dynamical core"=>"dynamical_core.md",
            "Parameterizations"=>"parametrizations.md",
            "Boundary conditions"=>"boundary_conditions.md",
            "New model setups"=>"new_model_setups.md",
            "Function and type index"=>"functions.md",
            "Style and convention guide"=>"conventions.md",
            "Development notes"=>"development.md"]
)

deploydocs(
    repo = "github.com/SpeedyWeather/SpeedyWeather.jl.git",
    devbranch = "main",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"]
)
