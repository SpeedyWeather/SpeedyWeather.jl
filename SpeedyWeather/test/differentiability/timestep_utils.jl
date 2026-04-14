function timestep_oop!(vars_new::Variables, vars_old::Variables, dt, model, lf1 = 2, lf2 = 2)
    copy!(vars_new, vars_old)
    SpeedyWeather.timestep!(vars_new, dt, model, lf1, lf2)
    return nothing
end

function timestep_oop!(vars_new::Variables, vars_old::Variables, dt, model, p::ComponentVector, lf1 = 2, lf2 = 2)
    copy!(vars_new, vars_old)
    new_model = SpeedyWeather.reconstruct(model, p)
    SpeedyWeather.timestep!(vars_new, dt, new_model, lf1, lf2)
    return nothing
end

# for FiniteDifferences.jl, we need to copy all inputs that are mutated
# because this function is called many times by FiniteDifferences
function timestep_oop(vars, dt, model, lf1 = 2, lf2 = 2)

    vars_copy = deepcopy(vars)
    model_copy = deepcopy(model) # just to be safe, as we have some temporary memory in there as well

    SpeedyWeather.timestep!(vars_copy, dt, model_copy, lf1, lf2)
    return vars_copy
end

function timestep_oop(vars, dt, model, p::ComponentVector, lf1 = 2, lf2 = 2)

    vars_copy = deepcopy(vars)
    model_copy = deepcopy(model) # just to be safe, as we have some temporary memory in there as well

    new_model = SpeedyWeather.reconstruct(model_copy, p)

    SpeedyWeather.timestep!(vars_copy, dt, new_model, lf1, lf2)
    return vars_copy
end
