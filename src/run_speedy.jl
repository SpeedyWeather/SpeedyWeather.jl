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
    if typeof(startdate) == DateTime model.clock.time = startdate end
    model.clock.n_days = n_days
    initialize!(model.clock,model.time_stepping)

    model.output.output = output            # enable/disable output
    initialize && initialize!(model)        # initialize again?

    # run it, yeah!
    time_stepping!(prognostic_variables,diagnostic_variables,model)
end

"""
    progn_vars = run_speedy(NF,Model;kwargs...)     or
    progn_vars = run_speedy(NF;kwargs...)           or
    progn_vars = run_speedy(Model;kwargs...)

Runs SpeedyWeather.jl with number format `NF` and the model `Model` and any additional parameters
in the keyword arguments `kwargs...`. Unspecified parameters use the default values."""
function run_speedy(::Type{NF} = DEFAULT_NF,                    # default number format
                    ::Type{Model} = DEFAULT_MODEL;              # default model
                    spectral_grid::NamedTuple = NamedTuple(),   # some keyword arguments to be
                    planet::NamedTuple = NamedTuple(),          # passed on
                    atmosphere::NamedTuple = NamedTuple(),
                    time_stepping::NamedTuple = NamedTuple(),
                    feedback::NamedTuple = NamedTuple(),
                    output::NamedTuple = NamedTuple(),
                    clock::NamedTuple = NamedTuple(),
                    kwargs...
                    ) where {NF<:AbstractFloat,Model<:ModelSetup}

    # pass on some keyword arguments to the default structs for convenience
    spectral_grid = SpectralGrid{Model}(;NF,spectral_grid...)
    planet = Earth(;planet...)
    atmosphere = EarthAtmosphere(;atmosphere...)
    time_stepping = Leapfrog(spectral_grid;time_stepping...)
    clock = Clock(time_stepping;clock...)
    output = OutputWriter(spectral_grid;output...)
    feedback = Feedback(output;feedback...)

    # create model with mostly defaults and initalize
    ConcreteModel = default_concrete_model(Model)
    model = ConcreteModel(;spectral_grid,planet,atmosphere,time_stepping,clock,
                            output,feedback,kwargs...)
    simulation = initialize!(model)

    # run it, yeah!
    (;prognostic_variables, diagnostic_variables, model) = simulation
    time_stepping!(prognostic_variables,diagnostic_variables,model)
    return simulation.prognostic_variables  # return prognostic variables when finished
end

# if only Model M provided, use default number format NF
run_speedy(::Type{Model};kwargs...) where {Model<:ModelSetup} = run_speedy(DEFAULT_NF,Model;kwargs...)
