function timestep_oop!(progn_new::PrognosticVariables, progn_old::PrognosticVariables, diagn, dt, model, lf1=2, lf2=2)
    copy!(progn_new, progn_old)
    SpeedyWeather.timestep!(progn_new, diagn, dt, model, lf1, lf2)
    return nothing
end 

# for FiniteDifferences.jl, we need to copy all inputs that are mutated 
# because this function is called many times by FiniteDifferences
function timestep_oop(progn, diagn, dt, model, lf1=2, lf2=2)

    progn_copy = deepcopy(progn)
    diagn_copy = deepcopy(diagn)

    model_copy = deepcopy(model) # just to be save, as we have some temporary memory in there as well

    SpeedyWeather.timestep!(progn_copy, diagn_copy, dt, model_copy, lf1, lf2)
    return progn_copy
end 
