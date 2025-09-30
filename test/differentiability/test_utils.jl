# additional utility functions and types

using SpeedyWeather: AbstractSimulation

# TODO: think about how this might be adapted into an extension simulation type;
# this type could potentially also implement the `AbstractSimulation` interface?
struct ADSimulation{
    ModelType<:AbstractModel,
    PrognosticType,
    DiagnosticType,
}
    model::ModelType
    progvars::PrognosticType
    diagvars::DiagnosticType
    dprogvars::PrognosticType
    ddiagvars::DiagnosticType
end

function ADSimulation(sim::AbstractSimulation)
    (; prognostic_variables, diagnostic_variables, model) = simulation
    progn = deepcopy(prognostic_variables)
    diagn = deepcopy(diagnostic_variables)
    return ADSimulation(
        deepcopy(model),
        progn,
        diagn,
        zero(progn),
        make_zero(diagn),
    )
end

prognosticseed(adsim::ADSimulation) = deepcopy(adsim.progvars), one(adsim.progvars)

diagnosticseed(adsim::ADSimulation) = deepcopy(adsim.diagvars), one(adsim.diagvars)

function Base.one(diag::DiagnosticVariables{NF}) where NF
    vec, re = to_vec(diag)
    vec .= NF(1)
    return re(vec)
end 

function initialize_with_spinup!(model::AbstractModel, spinup_period=Day(5), init_period=Day(1))
    simulation = initialize!(model)  
    initialize!(simulation)
    run!(simulation, period=spinup_period) # spin-up to get nonzero values for all fields
    initialize!(simulation; period=init_period)
    return simulation
end
