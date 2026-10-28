abstract type AbstractDiagnostic <: AbstractModelComponent end

export AbstractDiagnostic, Accumulate, Maximum, Mean, PrecipitationDiagnostics

# fallbacks: no variables to initialize, nothing to diagnose or reset. Arguments other than the
# diagnostic are untyped so that methods for custom diagnostics, typed or untyped, are never ambiguous.
# Two initialize! methods: `initialize!(diagnostic, model)` with `initialize!(model)` before the variables
# are allocated (fallback for AbstractModelComponent), and `initialize!(vars, diagnostic, model)` after
# allocation with the initial conditions (but before `StartFromFile`) to initialize `variables.output`
initialize!(vars, ::AbstractDiagnostic, model) = nothing
diagnostic!(vars, ::AbstractDiagnostic, model) = nothing
reset_diagnostic!(vars, ::AbstractDiagnostic, model) = nothing

"""$(TYPEDSIGNATURES)
Calculate all diagnostics in `model.diagnostics`, called once per time step (including spin-up
steps) after the time step but before the clock is stepped forward and the output is written.
Diagnostics only write into `variables.output` which nothing in the model reads, so they cannot
change the simulation."""
diagnostics!(vars::AbstractVariables, model::AbstractModel) =
    hasproperty(model, :diagnostics) ? diagnostics!(vars, model.diagnostics, model) : nothing

"""$(TYPEDSIGNATURES)
Reset all diagnostics in `model.diagnostics`, called by the output writer after writing, so
that the output interval determines the period over which the diagnostics are taken, e.g.
precipitation accumulated since the last output, or the maximum temperature over the output interval."""
reset_diagnostics!(vars::AbstractVariables, model::AbstractModel) =
    hasproperty(model, :diagnostics) ? reset_diagnostics!(vars, model.diagnostics, model) : nothing
reset_diagnostics!(simulation::AbstractSimulation) = reset_diagnostics!(simulation.variables, simulation.model)

# unrolled at compile time over the diagnostics in the NamedTuple
for (all_func!, func!) in ((:diagnostics!, :diagnostic!), (:reset_diagnostics!, :reset_diagnostic!))
    @eval @generated function $all_func!(vars::AbstractVariables, diagnostics::NamedTuple{K}, model::AbstractModel) where {K}
        calls = [:($$(QuoteNode(func!))(vars, getfield(diagnostics, $(QuoteNode(k))), model)) for k in K]
        return Expr(:block, calls..., :(return nothing))
    end
end

# variable `name` in `group`, within `namespace` unless empty
get_variable(group::NamedTuple, namespace::Symbol, name::Symbol) =
    namespace == Symbol() ? getfield(group, name) : getfield(getfield(group, namespace), name)

# whether variable `name` exists in `group`, within `namespace` unless empty
has_variable(group::NamedTuple, namespace::Symbol, name::Symbol) =
    namespace == Symbol() ? haskey(group, name) :
    haskey(group, namespace) && haskey(getfield(group, namespace), name)

# warn on initialization if the variable a diagnostic is computed from does not exist
function warn_missing_variable(vars::AbstractVariables, group::Symbol, namespace::Symbol, name::Symbol, diagnostic)
    path = namespace == Symbol() ? "variables.$group.$name" : "variables.$group.$namespace.$name"
    has_variable(getfield(vars, group), namespace, name) ||
        @warn "$path not defined, $(nameof(typeof(diagnostic))) diagnostic is skipped."
    return nothing
end

"""Diagnostic that accumulates a rate [unit/s] in `variables.parameterizations` over time,
i.e. adds `Δt * rate` every time step into `variables.output`, e.g. precipitation [m] from the
precipitation rate [m/s]. `Δt` is the model time that passes per time step. Reset to zero
by the output writer, so it is the accumulation since the last output.
Silently skipped if the rate does not exist (warns on initialization). Fields are $(TYPEDFIELDS)"""
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

variables(accumulate::Accumulate) = (
    OutputVariable(
        accumulate.name, Grid2D(), units = accumulate.units, desc = accumulate.desc,
        namespace = accumulate.namespace
    ),
)

initialize!(vars::AbstractVariables, accumulate::Accumulate, ::AbstractModel) =
    warn_missing_variable(vars, :parameterizations, accumulate.namespace, accumulate.rate, accumulate)

function diagnostic!(vars::AbstractVariables, accumulate::Accumulate, model::AbstractModel)
    (; rate, name, namespace) = accumulate
    has_variable(vars.parameterizations, namespace, rate) || return nothing
    rate_field = get_variable(vars.parameterizations, namespace, rate)
    accumulated = get_variable(vars.output, namespace, name)
    Δt = elapsed_time(model.time_stepping, vars.prognostic.clock)   # model time [s] of this time step
    accumulated .+= Δt .* rate_field
    return nothing
end

reset_diagnostic!(vars::AbstractVariables, accumulate::Accumulate, ::AbstractModel) =
    fill!(get_variable(vars.output, accumulate.namespace, accumulate.name), 0)

"""Diagnostic for the maximum of a 2D grid variable, e.g. the surface air temperature, over the
output interval: The output writer resets the maximum (to -Inf) after writing, so with daily output
this is the daily maximum. NaNs in the variable propagate into the maximum.
Silently skipped if the variable does not exist (warns on initialization). Fields are $(TYPEDFIELDS)"""
@kwdef struct Maximum <: AbstractDiagnostic
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
end

variables(maximum::Maximum) = (
    OutputVariable(
        maximum.name, Grid2D(), units = maximum.units, namespace = maximum.namespace,
        desc = "Maximum of $(maximum.variable) since the last output"
    ),
)

function initialize!(vars::AbstractVariables, maximum::Maximum, model::AbstractModel)
    warn_missing_variable(vars, maximum.group, maximum.namespace, maximum.variable, maximum)
    reset_diagnostic!(vars, maximum, model)
    return nothing
end

function diagnostic!(vars::AbstractVariables, maximum::Maximum, ::AbstractModel)
    (; variable, group, namespace, name) = maximum
    has_variable(getfield(vars, group), namespace, variable) || return nothing
    field = get_variable(getfield(vars, group), namespace, variable)
    running_maximum = get_variable(vars.output, namespace, name)
    running_maximum .= max.(running_maximum, field)     # max(-Inf, x) = x after reset
    return nothing
end

reset_diagnostic!(vars::AbstractVariables, maximum::Maximum, ::AbstractModel) =
    fill!(get_variable(vars.output, maximum.namespace, maximum.name), -Inf)

"""Diagnostic for the time mean of a 2D grid variable, e.g. the surface air temperature, over the
output interval: The output writer resets the mean after writing, so with weekly output this is
the weekly mean. Updated incrementally so that it can be read (e.g. by output) at any time step.
NaN before the first time step after a reset. Every time step has the same weight.
Silently skipped if the variable does not exist (warns on initialization). Fields are $(TYPEDFIELDS)"""
@kwdef struct Mean <: AbstractDiagnostic
    "[OPTION] Name of the 2D grid variable to take the mean of"
    variable::Symbol = :surface_air_temperature

    "[OPTION] Variable group of `variable`, e.g. `:parameterizations` or `:grid`"
    group::Symbol = :parameterizations

    "[OPTION] Namespace of `variable` and of the mean, e.g. `:land`"
    namespace::Symbol = Symbol()

    "[OPTION] Name of the mean in `variables.output`"
    name::Symbol = Symbol(variable, :_mean)

    "[OPTION] Units of the mean"
    units::String = "K"
end

count_name(mean::Mean) = Symbol(mean.name, :_count)

variables(mean::Mean) = (
    OutputVariable(
        mean.name, Grid2D(), units = mean.units, namespace = mean.namespace,
        desc = "Mean of $(mean.variable) since the last output"
    ),
    OutputVariable(
        count_name(mean), ScalarDim(), units = "1", namespace = mean.namespace,
        desc = "Number of time steps in the mean of $(mean.variable)"
    ),
)

function initialize!(vars::AbstractVariables, mean::Mean, model::AbstractModel)
    warn_missing_variable(vars, mean.group, mean.namespace, mean.variable, mean)
    reset_diagnostic!(vars, mean, model)
    return nothing
end

function diagnostic!(vars::AbstractVariables, mean::Mean, ::AbstractModel)
    (; variable, group, namespace, name) = mean
    has_variable(getfield(vars, group), namespace, variable) || return nothing
    field = get_variable(getfield(vars, group), namespace, variable)
    running_mean = get_variable(vars.output, namespace, name)
    count = get_variable(vars.output, namespace, count_name(mean))
    count[] += 1
    n = count[]
    running_mean .= ifelse.(n == 1, field, running_mean .+ (field .- running_mean) ./ n)
    return nothing
end

function reset_diagnostic!(vars::AbstractVariables, mean::Mean, ::AbstractModel)
    fill!(get_variable(vars.output, mean.namespace, mean.name), NaN)
    get_variable(vars.output, mean.namespace, count_name(mean))[] = 0
    return nothing
end

"""$(TYPEDSIGNATURES)
Default diagnostics of the `PrimitiveWetModel`: large-scale and convective rain and snow [m]
and river runoff [m] accumulated from their rates [m/s] since the last output."""
PrecipitationDiagnostics() = (
    rain_large_scale = Accumulate(
        rate = :rain_rate_large_scale, name = :rain_large_scale, units = "m",
        desc = "Large-scale precipitation (rain, accumulated since last output)"
    ),
    snow_large_scale = Accumulate(
        rate = :snow_rate_large_scale, name = :snow_large_scale, units = "m",
        desc = "Large-scale precipitation (snow, accumulated since last output)"
    ),
    rain_convection = Accumulate(
        rate = :rain_rate_convection, name = :rain_convection, units = "m",
        desc = "Convective precipitation (rain, accumulated since last output)"
    ),
    snow_convection = Accumulate(
        rate = :snow_rate_convection, name = :snow_convection, units = "m",
        desc = "Convective precipitation (snow, accumulated since last output)"
    ),
    river_runoff = Accumulate(
        rate = :river_runoff_rate, name = :river_runoff, namespace = :land, units = "m",
        desc = "River runoff from soil moisture (accumulated since last output)"
    ),
)
