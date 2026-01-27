import Pkg 
Pkg.activate("test/differentiability/sensitivity_examples")

using SpeedyWeather, Enzyme

spectral_grid = SpectralGrid(trunc=32, nlayers=8)          # define resolution
model = PrimitiveWetModel(; spectral_grid)                 # construct model
simulation = initialize!(model)  
initialize!(simulation)
run!(simulation, period=Day(20))
        
(; prognostic_variables, diagnostic_variables, model) = simulation
(; Δt, Δt_millisec) = model.time_stepping
dt = 2Δt

const i_point = 443 # pick this kind of random point (it's ≈ Copenhagen)

progn = prognostic_variables
diagn = diagnostic_variables

# do the scaling again because we need it for the timestepping when calling it manually 
SpeedyWeather.scale!(progn, diagn, model.planet.radius)

function temperature_at_gridpoint!(progn, diagn, model, dt)
    # just do a single timestep here for this example, later this will be replaced with run!(... , period=Day(N_days))
    SpeedyWeather.timestep!(progn, diagn, dt, model)
    return diagn.grid.temp_grid[i_point, 8]
end 

progn = simulation.prognostic_variables
dprogn = zero(progn)

diagn = simulation.diagnostic_variables
ddiag = make_zero(diagn)

# just to trigger some compilation
temperature_at_gridpoint!(progn, diagn, model, dt)

@time autodiff(Reverse, temperature_at_gridpoint!, Active, Duplicated(progn, dprogn), Duplicated(diagn, ddiag), Const(model), Const(dt))