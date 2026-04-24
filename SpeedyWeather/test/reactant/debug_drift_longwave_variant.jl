#=
Vary the longwave configuration to narrow down which part of OneBandLongwave
causes CPU↔Reactant drift.

Usage: julia --project debug_drift_longwave_variant.jl <variant>

Variants:
  transparent    — OneBandLongwave with TransparentLongwaveTransmissivity (t≡1, no exp)
  constant       — OneBandLongwave with ConstantLongwaveTransmissivity    (single exp per layer, no σ^4)
  frierson       — OneBandLongwave with FriersonLongwaveTransmissivity    (default, baseline)
=#
cd(@__DIR__)
import Pkg
Pkg.activate(@__DIR__)

using SpeedyWeather
using Reactant
using LinearAlgebra
using CUDA
using Statistics: mean
import SpeedyWeather: ReactantDevice, first_timesteps!, later_timestep!

include("setup.jl")

const TRUNC = 31
const NSTEPS = 10
ModelType = PrimitiveWetModel

variant = length(ARGS) >= 1 ? ARGS[1] : "frierson"

function make_longwave(SG, variant)
    if variant == "transparent"
        return OneBandLongwave(SG; transmissivity = TransparentLongwaveTransmissivity(SG))
    elseif variant == "constant"
        return OneBandLongwave(SG; transmissivity = ConstantLongwaveTransmissivity(SG))
    elseif variant == "constant_low_t"
        # low transmissivity forces small per-layer t similar to Frierson deep layers
        return OneBandLongwave(SG; transmissivity = ConstantLongwaveTransmissivity(SG; transmissivity = 0.2))
    elseif variant == "frierson"
        return OneBandLongwave(SG)  # default = FriersonLongwaveTransmissivity
    else
        error("unknown variant $variant")
    end
end

println("="^70)
println("Longwave variant test — $variant")
println("="^70)

println("\n[1] Creating CPU model...")
arch_c = SpeedyWeather.CPU()
SG_c = SpectralGrid(; nlayers = 8, trunc = TRUNC, architecture = arch_c)
M_c = MatrixSpectralTransform(SG_c)
ic_c = InitialConditions(SG_c, ModelType)
model_c = ModelType(
    SG_c;
    spectral_transform = M_c,
    feedback = nothing,
    dynamics = false,
    dynamics_only = false,
    initial_conditions = ic_c,
    convection = nothing,
    longwave_radiation = make_longwave(SG_c, variant),
)
sim_c = initialize!(model_c)

println("[2] Creating Reactant model...")
arch_r = SpeedyWeather.ReactantDevice()
SG_r = SpectralGrid(; architecture = arch_r, nlayers = 8, trunc = TRUNC)
M_r = MatrixSpectralTransform(SG_r)
ic_r = InitialConditions(SG_r, ModelType)
model_r = ModelType(
    SG_r;
    spectral_transform = M_r,
    feedback = nothing,
    dynamics = false,
    dynamics_only = false,
    initial_conditions = ic_r,
    convection = nothing,
    longwave_radiation = make_longwave(SG_r, variant),
)
sim_r = initialize!(model_r)

println("[3] Spin-up CPU 1 day, sync to Reactant...")
run!(sim_c; period = Day(1))
SpeedyWeather.copy!(sim_r.variables, sim_c.variables)

println("[4] Compiling Reactant kernels...")
initialize!(sim_r; steps = NSTEPS)
r_first! = @compile first_timesteps!(sim_r)
r_later! = @compile later_timestep!(sim_r)

println("[5] Resync CPU from Reactant post-compile...")
SpeedyWeather.copy!(sim_c.variables, sim_r.variables)

println("[6] Initialize for $NSTEPS steps...")
initialize!(sim_c; steps = NSTEPS)
initialize!(sim_r; steps = NSTEPS)

println("[7] Running CPU...")
SpeedyWeather.run!(sim_c; steps = NSTEPS)

println("[8] Running Reactant...")
SpeedyWeather.time_stepping!(sim_r, r_first!, r_later!)
SpeedyWeather.finalize!(sim_r)
SpeedyWeather.finalize!(sim_c)

println("[9] Comparing...")
T_c = Array(sim_c.variables.grid.temperature)
T_r = Array(sim_r.variables.grid.temperature)
q_c = Array(sim_c.variables.grid.humidity)
q_r = Array(sim_r.variables.grid.humidity)

T_max = maximum(abs.(T_c .- T_r))
q_max = maximum(abs.(q_c .- q_r))

# DEBUG: compare transmissivity and flux snapshots (debug_t/U/D) layer by layer
t_c = Array(sim_c.variables.scratch.grid.debug_t)
t_r = Array(sim_r.variables.scratch.grid.debug_t)
U_c = Array(sim_c.variables.scratch.grid.debug_U)
U_r = Array(sim_r.variables.scratch.grid.debug_U)
D_c = Array(sim_c.variables.scratch.grid.debug_D)
D_r = Array(sim_r.variables.scratch.grid.debug_D)

function compare_layer(name, A_c, A_r)
    diff = abs.(A_c .- A_r)
    println("\n--- debug_$(name) comparison ---  shape=", size(A_c))
    println("  global max|Δ|: ", maximum(diff))
    for k in 1:size(A_c, 2)
        dk = @view diff[:, k]
        ij_max = argmax(dk)
        rel = iszero(A_c[ij_max, k]) ? 0.0 : abs(dk[ij_max] / A_c[ij_max, k])
        println("  layer k=$k  max|Δ|=", dk[ij_max], "  (rel=", rel, ")  at ij=$ij_max  (CPU=", A_c[ij_max, k], "  Reactant=", A_r[ij_max, k], ")")
    end
    return maximum(diff)
end

t_max = compare_layer("t", t_c, t_r)
U_max = compare_layer("U", U_c, U_r)
D_max = compare_layer("D", D_c, D_r)

# Also probe ij=1 (north pole) per-layer
println("\n--- ij=1 (north pole) per-layer ---")
for k in 1:size(t_c, 2)
    println("  k=$k  t_Δ=", t_c[1, k] - t_r[1, k],
            "  U_Δ=", U_c[1, k] - U_r[1, k],
            "  D_Δ=", D_c[1, k] - D_r[1, k])
end

println("\n" * "="^70)
println("RESULT  longwave_variant=$variant")
println("  T_max_diff: $T_max")
println("  q_max_diff: $q_max")
println("  debug_t max|Δ|: ", t_max)
println("  debug_U max|Δ|: ", U_max)
println("  debug_D max|Δ|: ", D_max)
println("="^70)
