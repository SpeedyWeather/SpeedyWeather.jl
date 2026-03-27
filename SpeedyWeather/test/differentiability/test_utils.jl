# additional utility functions and types

using SpeedyWeather: AbstractSimulation

# TODO: think about how this might be adapted into an extension simulation type;
# this type could potentially also implement the `AbstractSimulation` interface?
struct ADSimulation{
        ModelType <: AbstractModel,
        VarsType,
    }
    model::ModelType
    vars::VarsType
    dvars::VarsType
end

function ADSimulation(sim::AbstractSimulation)
    (; variables, model) = sim
    vars = deepcopy(variables)
    return ADSimulation(
        deepcopy(model),
        vars,
        make_zero(vars),
    )
end

# expose compat accessors so existing callers still work
vars(adsim::ADSimulation) = adsim.vars
dvars(adsim::ADSimulation) = adsim.dvars

function ADseed(adsim::ADSimulation, name::Symbol)
    seed = make_zero(adsim.vars)

    # seed dvars_new with ones (output seed)
    for k in keys(getfield(seed, name))
        field = getfield(getfield(seed, name), k)
        if field isa AbstractArray
            field .= one(eltype(field))
        end
    end

    return deepcopy(adsim.vars), seed
end 

function initialize_with_spinup!(model::AbstractModel, spinup_period = Day(5), init_period = Day(1))
    simulation = initialize!(model)
    initialize!(simulation)
    run!(simulation, period = spinup_period) # spin-up to get nonzero values for all fields
    initialize!(simulation; period = init_period)
    return simulation
end