# define all variables for output
include("dynamics.jl")        # collected as DynamicsOutput()
include("precipitation.jl")   # PrecipitationOutput()
include("boundaries.jl")      # BoundaryOutput()
include("radiation.jl")       # RadiationOutput()
include("stochastic.jl")      # RandomPatternOutput()
include("surface_fluxes.jl")  # SurfaceFluxesOutput()
include("land.jl")            # LandOutput()
include("ocean.jl")           # SeaSurfaceTemperatureOutput()
include("tracers.jl")         # TracerOutput()
include("boundary_layer.jl")  # BoundaryLayerOutput()

# collect all together for conveneince
AllOutputVariables() = (
    DynamicsOutput()...,
    PrecipitationOutput()...,
    BoundaryOutput()...,
    RadiationOutput()...,
    RandomPatternOutput(),
    SurfaceFluxesOutput()...,
    LandOutput()...,
    OceanOutput()...,
    BoundaryLayerOutput()...,
)
