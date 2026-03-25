# Setup a model and compute sensitivyt of a single grid point

const N_steps = 10
const i_grid = 100

model = create_reactant_model(BarotropicModel, output = nothing)
simulation = initialize!(model)
run!(simulation; steps = 10)

# very high level entry point,
#= 


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

=#

# lower entry point for differentation
#initialize!(simulation, steps = 5)

function run_sim_lowlevel(simulation)
    SpeedyWeather.later_timestep!(simulation)
    return simulation.prognostic_variables.vor[i_grid, 1, 2]
end

dsimulation = make_zero(simulation)


#results = @jit autodiff(Enzyme.set_runtime_activity(Reverse), run_sim_lowlevel, Active, Duplicated(simulation, dsimulation))

@jit raise = true raise_first = true autodiff(Enzyme.set_runtime_activity(Reverse), SpeedyWeather.later_timestep!, Const, Duplicated(simulation, dsimulation))
