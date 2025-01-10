function timestep_oop!(progn_new, progn_old, diagn, dt, model)
    copy!(progn_new, progn_old)
    SpeedyWeather.timestep!(progn_new, diagn, dt, model)
    return nothing
end 

function timestep_oop(progn, diagn, dt, model)
    progn_new = zero(progn)
    timestep_oop!(progn_new, progn, diagn, dt, model)
    return progn_new
end 