"""
$(TYPEDSIGNATURES)
Run a SpeedyWeather.jl `simulation`. The `simulation.model` is assumed to be initialized,
otherwise use `initialize=true` as keyword argument."""
function run!(  simulation::Simulation;
                n_days::Real = 10,
                startdate::Union{Nothing,DateTime} = nothing,
                output::Bool = false)
    
    (;prognostic_variables, diagnostic_variables, model) = simulation
    (;clock) = prognostic_variables

    # set the clock
    if typeof(startdate) == DateTime
        clock.time = startdate
    end
    clock.n_days = n_days
    initialize!(clock,model.time_stepping)

    model.output.output = output            # enable/disable output

    # run it, yeah!
    time_stepping!(prognostic_variables,diagnostic_variables,model)
end