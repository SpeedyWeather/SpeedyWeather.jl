# Setup a model and compute sensitivyt of a single grid point

const N_steps = 10
const i_grid = 100

model = create_reactant_model(BarotropicModel, output=nothing)
simulation = initialize!(model)
run!(simulation; steps= 10)

function run_sim(simulation)
    run!(simulation; steps = N_steps)
    return simulation.prognostic_variables.vor[i_grid, 1, 2]
end

dsimulation = make_zero(simulation)

try
    results = @jit autodiff(Enzyme.set_runtime_activity(Reverse), run_sim, Active, Duplicated(simulation, dsimulation))
catch err
    @show err
end
