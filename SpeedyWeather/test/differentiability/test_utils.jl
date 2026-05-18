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

"""
    ADseed(adsim::ADSimulation, name::Symbol)

Return a copy of the `adsim.vars` and a seed for the `name` field of `adsim.dvars`. 
The seed should be the output of the function it is fed into.
"""
function ADseed(adsim::ADSimulation, name::Symbol)
    seed = make_zero(adsim.dvars)

    for k in keys(getfield(seed, name))
        field = getfield(getfield(seed, name), k)
        if field isa AbstractArray
            if eltype(field) <: Complex
                field .= one(eltype(real(field))) + im * one(eltype(real(field)))
            else
                field .= one(eltype(field))
            end
        elseif field isa NamedTuple
            for k2 in keys(field)
                field2 = getfield(field, k2)
                if eltype(field2) <: Complex
                    field2 .= one(eltype(real(field2))) + im * one(eltype(real(field2)))
                else
                    field2 .= one(eltype(field2))
                end
            end
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
