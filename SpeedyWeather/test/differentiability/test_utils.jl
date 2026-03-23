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
progvars(adsim::ADSimulation) = adsim.vars
diagvars(adsim::ADSimulation) = adsim.vars
dprogvars(adsim::ADSimulation) = adsim.dvars
ddiagvars(adsim::ADSimulation) = adsim.dvars

prognosticseed(adsim::ADSimulation) = deepcopy(adsim.vars), make_zero(adsim.vars)

diagnosticseed(adsim::ADSimulation) = deepcopy(adsim.vars), make_zero(adsim.vars)

function initialize_with_spinup!(model::AbstractModel, spinup_period = Day(5), init_period = Day(1))
    simulation = initialize!(model)
    initialize!(simulation)
    run!(simulation, period = spinup_period) # spin-up to get nonzero values for all fields
    initialize!(simulation; period = init_period)
    return simulation
end
