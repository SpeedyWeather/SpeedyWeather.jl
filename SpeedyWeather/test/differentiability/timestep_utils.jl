# Out-of-place wrappers around the time stepping for AD / finite-difference tests.
#
# NOTE (2026-06): ported to the post-`time_stepping/`-refactor API. The public single
# step is now `time_step!(vars, model.time_stepping, model)` — there is no externally
# supplied `dt` (it is derived from `vars.prognostic.clock` and `model.time_stepping`)
# and no leapfrog index `lf` (derived from the clock as well). The old
# `SpeedyWeather.timestep!(vars, dt, model, lf1, lf2)` no longer exists.

function timestep_oop!(vars_new::Variables, vars_old::Variables, model)
    copy!(vars_new, vars_old)
    SpeedyWeather.time_step!(vars_new, model.time_stepping, model)
    return nothing
end

function timestep_oop!(vars_new::Variables, vars_old::Variables, model, p::ComponentVector)
    copy!(vars_new, vars_old)
    new_model = SpeedyWeather.reconstruct(model, p)
    SpeedyWeather.time_step!(vars_new, new_model.time_stepping, new_model)
    return nothing
end

# for FiniteDifferences.jl, we need to copy all inputs that are mutated
# because this function is called many times by FiniteDifferences
function timestep_oop(vars, model)
    vars_copy = deepcopy(vars)
    model_copy = deepcopy(model) # just to be safe, as we have some temporary memory in there as well

    SpeedyWeather.time_step!(vars_copy, model_copy.time_stepping, model_copy)
    return vars_copy
end

function timestep_oop(vars, model, p::ComponentVector)
    vars_copy = deepcopy(vars)
    model_copy = deepcopy(model) # just to be safe, as we have some temporary memory in there as well

    new_model = SpeedyWeather.reconstruct(model_copy, p)

    SpeedyWeather.time_step!(vars_copy, new_model.time_stepping, new_model)
    return vars_copy
end
