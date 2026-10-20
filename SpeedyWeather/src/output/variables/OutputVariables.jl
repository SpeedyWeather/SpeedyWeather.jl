"""
Output variables to be added to an output writer like `NetCDFOutput` or `ZarrOutput`,
e.g. `add!(model, OutputVariables.Temperature(), OutputVariables.Humidity())`.
Every output variable is a subtype of `SpeedyWeather.AbstractOutputVariable` and
the module is exported but none of its names are, so always access them
qualified, e.g. `OutputVariables.Vorticity(keepbits = 7)`. Collections of output
variables are returned as tuples, e.g. `OutputVariables.Dynamics()`, `OutputVariables.Land()`
or all of them with `OutputVariables.All()`."""
module OutputVariables

using DocStringExtensions
import GPUArrays
import Dates: Hour

# abstract types and generic functions from SpeedyWeather extended or used here
import ..SpeedyWeather: AbstractOutput, AbstractOutputVariable, AbstractSimulation,
    AbstractVerticalExtrapolation, DryAdiabaticExtrapolation, NetCDFOutput,
    DEFAULT_COMPRESSION_LEVEL, DEFAULT_KEEPBITS, DEFAULT_SHUFFLE,
    define_variable!, field_view, finalize!, get_prognostic_step, hastime,
    interpolate_output!, is_land, output!, path, path_or_nothing, pressure,
    virtual_temperature, write_output!

include("dynamics.jl")        # collected as Dynamics()
include("precipitation.jl")   # Precipitation()
include("boundaries.jl")      # Boundary()
include("radiation.jl")       # Radiation()
include("stochastic.jl")      # RandomPattern()
include("surface_fluxes.jl")  # SurfaceFluxes()
include("land.jl")            # Land()
include("ocean.jl")           # Ocean()
include("tracers.jl")         # Tracer()
include("boundary_layer.jl")  # BoundaryLayer()

"""$(TYPEDSIGNATURES)
All output variables (except tracers and Terrarium) as a tuple,
e.g. `add!(model, OutputVariables.All()...)`."""
All() = (
    Dynamics()...,
    Precipitation()...,
    Boundary()...,
    Radiation()...,
    RandomPattern(),
    SurfaceFluxes()...,
    Land()...,
    Ocean()...,
    BoundaryLayer()...,
)

end
