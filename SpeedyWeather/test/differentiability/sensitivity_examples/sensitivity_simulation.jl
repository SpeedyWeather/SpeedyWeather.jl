import Pkg 
Pkg.activate(".")

using SpeedyWeather, Enzyme, JLD2, Checkpointing

# Parse command line argument for N (number of timesteps)
const N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 9

println("Running Sensitivity Analyis with N = $N")
savename_base = "sensitivity-$N"

spectral_grid = SpectralGrid(trunc=32, nlayers=8)          # define resolution
model = PrimitiveWetModel(; spectral_grid)                 # construct model
simulation = initialize!(model)  
initialize!(simulation)
run!(simulation, period=Day(20))

d_sim = make_zero(simulation)

# whereas calling reverse on the more low-level function shown in the other script mostly works
# and just has occational problems when rolling out longer trajectories
# this one just times out and never seems to complete even on a single time step
@time autodiff(Enzyme.Reverse, SpeedyWeather.timestep!, Const, Duplicated(simulation, d_sim))
