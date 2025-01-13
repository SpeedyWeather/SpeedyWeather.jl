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

# for comparisions with FiniteDifferences.jl we need inputs and outputs that are easy to convert
# to Vectors. PrognosticVariables, and even more so a Model struct are not really easy to convert
# Vectors. So, for the comparision we look at more limited problems. E.g. a single step and just 
# looking at the vorticiity
function step_vorticity!(vorticity_new::LowerTriangularArray, vorticity::LowerTriangularArray, progn, diagn, dt, model)
    #set!(progn, model.geometry, vor=vorticity, lf=1)
    progn.vor[1] .= vorticity
    SpeedyWeather.timestep!(progn, diagn, dt, model)
    copy!(vorticity_new, progn.vor[1])
    #vorticity_new.data .= progn.vor[1].data
    return nothing
end 

function step_vorticity(vorticity::LowerTriangularArray, progn, diagn, dt, model)
    vorticity_new = zero(vorticity)
    step_vorticity!(vorticity_new, vorticity, progn, diagn, dt, model)
    return vorticity_new
end