import Pkg 
Pkg.activate("test/differentiability/sensitivity_examples")

using SpeedyWeather, Enzyme, JLD2, Checkpointing

# Parse command line argument for N (number of timesteps)
const N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100

println("Running Sensitivity Analyis with N = $N")
savename_base = "sensitivity-$N"

spectral_grid = SpectralGrid(trunc=32, nlayers=8)          # define resolution
model = PrimitiveWetModel(; spectral_grid)                 # construct model
simulation = initialize!(model)  
initialize!(simulation)
run!(simulation, period=Day(20))
        
(; prognostic_variables, diagnostic_variables, model) = simulation
(; Δt, Δt_millisec) = model.time_stepping
dt = 2Δt

progn = prognostic_variables
diagn = diagnostic_variables

# do the scaling again because we need it for the timestepping when calling it manually 
SpeedyWeather.scale!(progn, diagn, model.planet.radius)

function checkpointed_timesteps!(progn::PrognosticVariables, diagn, model, N_steps, checkpoint_scheme::Scheme, lf1=2, lf2=2)
    
    @ad_checkpoint checkpoint_scheme for _ in 1:N_steps 
        SpeedyWeather.timestep!(progn, diagn, 2*model.time_stepping.Δt, model, lf1, lf2)
    end 

    return nothing  
end 

checkpoint_scheme = Revolve(N)

# Temperature One-Hot 
d_progn = zero(progn)
d_model = make_zero(model)
d_diag = make_zero(diagn)
d_diag.grid.temp_grid[443, 8] = 1

jldsave(string(savename_base,"temp-ic.jld2"); progn, diagn)

println("Starting sensitivity computation...")

# occasioanlly this gives a SegmentationFault (espacially on x86 and for large N), but not always
@time autodiff(Enzyme.Reverse, checkpointed_timesteps!, Const, Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Duplicated(model, d_model), Const(N), Const(checkpoint_scheme))

jldsave(string(savename_base,"temp.jld2"); d_progn)
jldsave(string(savename_base,"temp-fc.jld2"); progn, diagn)

# Precip One-Hot 
d_progn = zero(progn)
d_model = make_zero(model)
d_diag = make_zero(diagn)
d_diag.physics.total_precipitation_rate[443] = 1

jldsave(string(savename_base,"precip-ic.jld2"); progn, diagn)

println("Starting sensitivity computation...")

@time autodiff(Enzyme.Reverse, checkpointed_timesteps!, Const, Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Duplicated(model, d_model), Const(N), Const(checkpoint_scheme))

jldsave(string(savename_base,"precip.jld2"); d_progn)
jldsave(string(savename_base,"precip-fc.jld2"); progn, diagn)
