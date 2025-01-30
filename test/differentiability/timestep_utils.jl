function timestep_oop!(progn_new::PrognosticVariables, progn_old::PrognosticVariables, diagn, dt, model, lf1=2, lf2=2)
    copy!(progn_new, progn_old)
    SpeedyWeather.timestep!(progn_new, diagn, dt, model, lf1, lf2)
    return nothing
end 

function timestep_oop(progn, diagn, dt, model, lf1=2, lf2=2)
    progn_new = zero(progn)
    timestep_oop!(progn_new, progn, diagn, dt, model, lf1, lf2)
    return progn_new
end 
