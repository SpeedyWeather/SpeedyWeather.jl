import Pkg 
Pkg.activate("SpeedyWeather/test/differentiability/sensitivity_examples")

using SpeedyWeather, Enzyme, JLD2, Checkpointing

# Parse command line argument for N (number of timesteps)
const N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 5

println("Running Sensitivity Analyis with N = $N")
savename_base = "new-sensitivity-$N"

spectral_grid = SpectralGrid(trunc=32, nlayers=8)          # define resolution
model = PrimitiveWetModel(; spectral_grid)                 # construct model
simulation = initialize!(model)  
initialize!(simulation)
run!(simulation, period=Day(20))

# do the scaling again because we need it for the timestepping when calling it manually 
initialize!(simulation, steps=N)
(; prognostic_variables, diagnostic_variables, model) = simulation
(; Δt, Δt_millisec) = model.time_stepping
dt = 2Δt

const i_point = 443 # pick this kind of random point (it's ≈ Copenhagen)

progn = prognostic_variables
diagn = diagnostic_variables

function checkpointed_timesteps!(progn::PrognosticVariables, diagn, model, N_steps, checkpoint_scheme::Scheme, lf1=2, lf2=2)
    
    @ad_checkpoint checkpoint_scheme for _ in 1:N_steps 
        SpeedyWeather.timestep!(progn, diagn, 2*model.time_stepping.Δt, model, lf1, lf2)
    end 

    return diagn.grid.temp_grid[i_point, 8] 
end 

checkpoint_scheme = Revolve(N)

progn = simulation.prognostic_variables
dprogn = zero(progn)

diagn = simulation.diagnostic_variables
ddiag = make_zero(diagn)
dmodel = make_zero(model)

jldsave(string(savename_base,"temp-ic.jld2"); progn, diagn)

println("Starting sensitivity computation...")

# occasioanlly this gives a SegmentationFault (espacially on x86 and for large N), but not always
@time res = autodiff(Enzyme.Reverse, checkpointed_timesteps!, Active, Duplicated(progn, dprogn), Duplicated(diagn, ddiag), Duplicated(model, dmodel), Const(N), Const(checkpoint_scheme))

jldsave(string(savename_base,"temp.jld2"); dprogn)
