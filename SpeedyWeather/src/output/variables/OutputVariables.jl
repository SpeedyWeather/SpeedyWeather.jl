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
import KernelAbstractions: @kernel, @index
import SpeedyWeatherInternals.Architectures: architecture
import SpeedyWeatherInternals.KernelLaunching: launch!, LinearWorkOrder

# abstract types and generic functions from SpeedyWeather extended or used here
import ..SpeedyWeather: AbstractOutput, AbstractOutputVariable, AbstractSimulation,
    AbstractVerticalExtrapolation, DryAdiabaticExtrapolation, NetCDFOutput,
    define_variable!, field_view, finalize!, get_prognostic_step,
    interpolate_output!, on_output_architecture, output!, pressure, ConstantExtrapolation,
    virtual_temperature, write_output!

# defaults for output variables that don't define these fields themselves
const DEFAULT_MISSING_VALUE = NaN
const DEFAULT_COMPRESSION_LEVEL = 1
const DEFAULT_SHUFFLE = false
const DEFAULT_KEEPBITS = 15

"""$(TYPEDSIGNATURES) To be extended for every output variable to define
the path where in `simulation` to find that output variable `::AbstractField`,
e.g. `OutputVariables.path(::MyOutputVariable, simulation) = simulation.variables.grid.vorticity`."""
function path end

# properties of output variables, used by the output writers
is2D(variable::AbstractOutputVariable) = ~variable.dims_xyzt[3]
is3D(variable::AbstractOutputVariable) = variable.dims_xyzt[3]
is_land(variable::AbstractOutputVariable) = hasproperty(variable, :is_land) ? variable.is_land : false
hastime(variable::AbstractOutputVariable) = variable.dims_xyzt[4]

"""$(TYPEDSIGNATURES)
`true` if `variable` exists in `simulation` and can therefore be written out. This is the
same `path_or_nothing` check the generic [`output!`](@ref) applies at write time, exposed so
that output writers can apply it at *schema-definition* time too: a variable that is skipped
by `output!` but still defined in the file/store would otherwise be left behind as an
all-`missing_value` phantom array. Variables can be absent when the `Model` type an output
writer was constructed with declares defaults the actual model doesn't have, e.g. a
`PrimitiveWet` output attached to a `BarotropicModel` has no temp/humid/mslp."""
exists_in_simulation(variable::AbstractOutputVariable, simulation) =
    !isnothing(path_or_nothing(variable, simulation))

"""$(TYPEDSIGNATURES) Like `path` but returns `nothing` instead of throwing an error if the variable is not defined in the simulation."""
path_or_nothing(variable::AbstractOutputVariable, simulation) = try
    path(variable, simulation)
catch FieldError
    nothing
end

"""$(TYPEDSIGNATURES)
Extrapolation of `variable` beyond the model's vertical layers: its `extrapolation` field
if it has one (e.g. [`OutputVariables.Temperature`](@ref) descends dry-adiabatically below the
lowest model layer), otherwise constant. Same `hasproperty` pattern as `transform` and
`unscale` in the generic [`output!`](@ref)."""
output_extrapolation(variable::AbstractOutputVariable) =
    hasproperty(variable, :extrapolation) ? variable.extrapolation : ConstantExtrapolation()

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
