"""
$(TYPEDSIGNATURES)
Run a SpeedyWeather.jl `simulation`. The `simulation.model` is assumed to be initialized,
otherwise use `initialize=true` as keyword argument."""
function run!(  simulation::Simulation;
                initialize::Bool = false,
                n_days::Real = 10,
                startdate::Union{Nothing,DateTime} = nothing,
                output::Bool = false)
    
    (;prognostic_variables, diagnostic_variables, model) = simulation

    # set the clock
    if typeof(startdate) == DateTime prognostic_variables.clock.time = startdate end
    prognostic_variables.clock.n_days = n_days
    initialize!(prognostic_variables.clock,model.time_stepping)

    model.output.output = output            # enable/disable output
    initialize && initialize!(model)        # initialize again?

    # run it, yeah!
    time_stepping!(prognostic_variables,diagnostic_variables,model)
end