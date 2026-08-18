# use Julia 1.10 for this 

import Pkg
Pkg.activate("SpeedyWeather/test/differentiability/sensitivity_examples")

using SpeedyWeather, Enzyme, JLD2, Checkpointing, StatsBase

# Enzyme runs the LLVM Attributor pass on Julia < 1.12
# (`Enzyme.Compiler.RunAttributor = Ref(VERSION < v"1.12")`). Stepping the clock inside the
# `@ad_checkpoint` loop sends the Attributor's AAPotentialValues analysis into unbounded
# recursion (`AAPotentialValuesFloating::updateImpl` -> `getAssumedSimplified` -> ... -> itself),
# which overflows the C++ stack and shows up as a segfault. Disabling the pass avoids it.
Enzyme.Compiler.RunAttributor[] = false

# Parse command line argument for N (number of timesteps)
const N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 5

println("Running Sensitivity Analyis with N = $N")
savename_base = "new-sensitivity-$N"

spectral_grid = SpectralGrid(truncation = 33, nlayers = 8)          # define resolution
model = PrimitiveWetModel(; spectral_grid)                 # construct model
simulation = initialize!(model)
initialize!(simulation)
run!(simulation, period = Day(20))

# do the scaling again because we need it for the timestepping when calling it manually
initialize!(simulation, steps = N)
(; variables, model) = simulation
vars = variables

function checkpointed_timesteps!(vars::Variables, model, N_steps, checkpoint_scheme::Scheme)

     @ad_checkpoint checkpoint_scheme for _ in 1:N_steps
        SpeedyWeather.time_step!(vars, model.time_stepping, model)     # calculate tendencies and step forward
        SpeedyWeather.time_step!(vars.prognostic.clock, model.time_stepping)                # then step the clock forward
    end

    return mean(vars.grid.temperature[:, 8, 2])
end

checkpoint_scheme = Revolve(N)

dvars = make_zero(vars)
dmodel = make_zero(model)

println("Starting sensitivity computation...")

@time res = autodiff(Enzyme.Reverse, checkpointed_timesteps!, Active, Duplicated(vars, dvars), Duplicated(model, dmodel), Const(N), Const(checkpoint_scheme))
