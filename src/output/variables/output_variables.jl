# define all variables for output
include("dynamics.jl")
include("precipitation.jl")   # collected as PrecipitationOutput()
include("boundaries.jl")      # BoundaryOutput()
include("radiation.jl")       # RadiationOutput()
include("stochastic.jl")      # RandomPatternOutput()
include("surface_fluxes.jl")  # SurfaceFluxesOutput()
include("land.jl")            # LandOutput()
include("ocean.jl")           # SeaSurfaceTemperatureOutput()
include("tracers.jl")         # TracerOutput()

# collect all together for conveneince
AllOutputVariables() = (
    PrecipitationOutput()...,
    BoundaryOutput()...,
    RadiationOutput()...,
    RandomPatternOutput(),
    SurfaceFluxesOutput()...,
    LandOutput()...,
    OceanOutput(),
)
