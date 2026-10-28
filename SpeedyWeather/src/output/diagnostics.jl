abstract type AbstractDiagnostic <: AbstractModelComponent end

export AbstractDiagnostic, Accumulate, PeriodicMaximum, PrecipitationDiagnostics

# fallbacks: no variables to initialize, nothing to diagnose. Arguments other than the diagnostic
# are untyped so that methods for custom diagnostics, typed or untyped, are never ambiguous
initialize!(vars, ::AbstractDiagnostic, model) = nothing
diagnostic!(vars, ::AbstractDiagnostic, model) = nothing

"""$(TYPEDSIGNATURES)
Calculate all diagnostics in `model.diagnostics`, called once per time step (including spin-up
steps) after the clock has been stepped forward, before the output is written. Diagnostics only
write into `variables.output` which nothing in the model reads, so they cannot change the simulation."""
diagnostics!(vars::AbstractVariables, model::AbstractModel) =
    hasproperty(model, :diagnostics) ? diagnostics!(vars, model.diagnostics, model) : nothing

function diagnostics!(vars::AbstractVariables, diagnostics::NamedTuple, model::AbstractModel)
    foreach(diagnostic -> diagnostic!(vars, diagnostic, model), values(diagnostics))
    return nothing
end

# variable `name` in `group`, within `namespace` unless empty
get_variable(group::NamedTuple, namespace::Symbol, name::Symbol) =
    namespace == Symbol() ? getfield(group, name) : getfield(getfield(group, namespace), name)

"""Diagnostic that accumulates a rate [unit/s] in `variables.parameterizations` over time,
i.e. adds `Δt * rate` every time step into `variables.output`, e.g. accumulated precipitation [m]
from the precipitation rate [m/s]. The accumulation is never reset, also not across restarts.
Fields are $(TYPEDFIELDS)"""
@kwdef struct Accumulate <: AbstractDiagnostic
    "[OPTION] Name of the rate [unit/s] in `variables.parameterizations` to accumulate"
    rate::Symbol

    "[OPTION] Name of the accumulated variable in `variables.output`"
    name::Symbol

    "[OPTION] Namespace of both rate and accumulated variable, e.g. `:land`"
    namespace::Symbol = Symbol()

    "[OPTION] Units of the accumulated variable"
    units::String = ""

    "[OPTION] Description of the accumulated variable"
    desc::String = ""
end

# also declare the rate so that the accumulation works (as zero) without the component writing it
variables(accumulate::Accumulate) = (
    ParameterizationVariable(accumulate.rate, Grid2D(), namespace = accumulate.namespace),
    OutputVariable(
        accumulate.name, Grid2D(), units = accumulate.units, desc = accumulate.desc,
        namespace = accumulate.namespace
    ),
)

function diagnostic!(vars::AbstractVariables, accumulate::Accumulate, model::AbstractModel)
    rate = get_variable(vars.parameterizations, accumulate.namespace, accumulate.rate)
    accumulated = get_variable(vars.output, accumulate.namespace, accumulate.name)
    (; Δt) = model.time_stepping    # the time step [s] the parameterizations use for their rates
    accumulated .+= Δt .* rate
    return nothing
end

"""Diagnostic for the maximum of a 2D grid variable over a period, e.g. the daily maximum
surface air temperature. Keeps a running maximum `<name>_running` in `variables.output` which is
copied into `variables.output.<name>` and reset at the end of every period (as determined by
`schedule`, daily by default). Before the first period has passed `<name>` is NaN.
Fields are $(TYPEDFIELDS)"""
@kwdef struct PeriodicMaximum{S} <: AbstractDiagnostic
    "[OPTION] Name of the 2D grid variable to take the maximum of"
    variable::Symbol = :surface_air_temperature

    "[OPTION] Variable group of `variable`, e.g. `:parameterizations` or `:grid`"
    group::Symbol = :parameterizations

    "[OPTION] Namespace of `variable` and of the maximum, e.g. `:land`"
    namespace::Symbol = Symbol()

    "[OPTION] Name of the maximum in `variables.output`"
    name::Symbol = Symbol(variable, :_maximum)

    "[OPTION] Units of the maximum"
    units::String = "K"

    "[OPTION] Schedule when the maximum is taken and reset, every day by default"
    schedule::S = PeriodicSchedule(every = Day(1))
end

running_name(maximum::PeriodicMaximum) = Symbol(maximum.name, :_running)

variables(maximum::PeriodicMaximum) = (
    OutputVariable(
        maximum.name, Grid2D(), units = maximum.units, namespace = maximum.namespace,
        desc = "Maximum of $(maximum.variable) over the last period"
    ),
    OutputVariable(
        running_name(maximum), Grid2D(), units = maximum.units, namespace = maximum.namespace,
        desc = "Running maximum of $(maximum.variable) since the start of the period"
    ),
)

initialize!(maximum::PeriodicMaximum, model::AbstractModel) = initialize!(maximum.schedule, model.time_stepping)

function initialize!(vars::AbstractVariables, maximum::PeriodicMaximum, ::AbstractModel)
    fill!(get_variable(vars.output, maximum.namespace, maximum.name), NaN)                 # no maximum yet
    fill!(get_variable(vars.output, maximum.namespace, running_name(maximum)), -Inf)       # max(-Inf, x) = x
    return nothing
end

function diagnostic!(vars::AbstractVariables, maximum::PeriodicMaximum, ::AbstractModel)
    variable = get_variable(getfield(vars, maximum.group), maximum.namespace, maximum.variable)
    period_maximum = get_variable(vars.output, maximum.namespace, maximum.name)
    running = get_variable(vars.output, maximum.namespace, running_name(maximum))

    running .= max.(running, variable)

    # at the end of the period: pass on the running maximum and reset it, branchless for Reactant
    end_of_period = isscheduled(maximum.schedule, vars.prognostic.clock)
    period_maximum .= ifelse.(end_of_period, running, period_maximum)
    running .= ifelse.(end_of_period, -Inf, running)
    return nothing
end

"""$(TYPEDSIGNATURES)
Default diagnostics of the `PrimitiveWetModel`: accumulated large-scale and convective rain and
snow [m] and river runoff [m] from their rates [m/s]."""
PrecipitationDiagnostics() = (
    rain_large_scale = Accumulate(
        rate = :rain_rate_large_scale, name = :rain_large_scale, units = "m",
        desc = "Large-scale precipitation (rain, accumulated)"
    ),
    snow_large_scale = Accumulate(
        rate = :snow_rate_large_scale, name = :snow_large_scale, units = "m",
        desc = "Large-scale precipitation (snow, accumulated)"
    ),
    rain_convection = Accumulate(
        rate = :rain_rate_convection, name = :rain_convection, units = "m",
        desc = "Convective precipitation (rain, accumulated)"
    ),
    snow_convection = Accumulate(
        rate = :snow_rate_convection, name = :snow_convection, units = "m",
        desc = "Convective precipitation (snow, accumulated)"
    ),
    river_runoff = Accumulate(
        rate = :river_runoff_rate, name = :river_runoff, namespace = :land, units = "m",
        desc = "River runoff from soil moisture (accumulated)"
    ),
)
