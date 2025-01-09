# define all variables for output
include("variables/dynamics.jl")
include("variables/precipitation.jl")   # collected as PrecipitationOutput()
include("variables/boundaries.jl")      # BoundaryOutput()
include("variables/radiation.jl")       # RadiationOutput()
include("variables/stochastic.jl")      # RandomPatternOutput()
include("variables/surface_fluxes.jl")  # SurfaceFluxesOutput()
include("variables/ocean.jl")           # SeaSurfaceTemperatureOutput()

# collect all together for conveneince
AllOutputVariables() = (
    PrecipitationOutput()...,
    BoundaryOutput()...,
    RadiationOutput()...,
    RandomPatternOutput(),
    SurfaceFluxesOutput()...,
    SeaSurfaceTemperatureOutput(),
)
