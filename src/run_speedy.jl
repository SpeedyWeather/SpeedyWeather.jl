"""
$(TYPEDSIGNATURES)
Run a SpeedyWeather.jl `simulation`. The `simulation.model` is assumed to be initialized."""
function run!(  simulation::Simulation;
                period = Day(10),
                output::Bool = false, 
                n_days::Union{Nothing, Real}=nothing)
    
    if !isnothing(n_days)
        @warn "run!: n_days keyword is deprecated, use period = Day(n_days) instead."
        period = Day(n_days) 
    end 

    (;prognostic_variables, diagnostic_variables, model) = simulation
    (;clock) = prognostic_variables

    # set the clock's enddate
    set_period!(clock,period)
    initialize!(clock,model.time_stepping)

    model.output.output = output            # enable/disable output

    # run it, yeah!
    time_stepping!(prognostic_variables,diagnostic_variables,model)
end
