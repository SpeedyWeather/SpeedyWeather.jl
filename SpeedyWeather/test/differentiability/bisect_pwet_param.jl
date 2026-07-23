# =============================================================================
# Parameter-AD bisection driver for the Julia 1.12 Enzyme codegen assertion
#   "unhandled accumulate with partial sizes"  (DiffeGradientUtils.cpp:433)
# =============================================================================
# full_diff_CI.jl (full PrimitiveWet time_step! WITH physics, Duplicated(model))
# ABORTS on x86 Julia 1.12 (exit 134) but PASSES on 1.10. State AD (Const model)
# was already fixed by the view-preserving make_zero(::Variables). The remaining
# trigger is the PARAMETER-AD (Duplicated(model)) reverse of some component.
#
# This driver runs EXACTLY ONE autodiff call, selected by ARGS[1], so that an
# Enzyme codegen abort (which kills the process) isolates that single call.
# It reuses the existing test harness (test_utils.jl: ADSimulation/ADseed) so the
# variants line up with the Const(model) tests in primitivewet.jl — the ONLY
# difference is Const(model) -> Duplicated(model) (parameter AD).
#
# The abort is a COMPILE-time event (CreatePrimalAndGradient), so it is
# independent of seed values / resolution — trunc=5,nlayers=1 (== full_diff_CI)
# reproduces it. An aborted variant = the culprit component; a variant that prints
# "RESULT: SUCCESS" compiled its reverse fine.
#
# Run one variant per process (an abort kills the process, isolating that variant).
# Set the env up ONCE first with setup_diff_env.jl. Then either:
#   - all variants, sequential:  run_bisect_local.sh
#   - SLURM array (extensive):   sbatch bisect_pwet_param.slurm
#   - single variant, directly:
#       julia --project=SpeedyWeather/test/differentiability \
#             SpeedyWeather/test/differentiability/bisect_pwet_param.jl <variant>
#
# Variants (see BISECTION_x86.md for the decision tree):
#   A1_state       full time_step!, Const(model)            [state AD, expect PASS]
#   A2_param       full time_step!, Duplicated(model)       [param AD, expect ABORT]  == full_diff_CI
#   B1_dynonly     full time_step! dynamics_only, Duplicated(model)  [physics vs dyn-core]
#   C1_dyntend     dynamics_tendencies!,        Duplicated(model)
#   C2_paramtend   parameterization_tendencies!, Duplicated(model)   [physics]
#   C3_implicit    implicit_correction!,        Duplicated(model)
#   C4_diffusion   horizontal_diffusion!,       Duplicated(model)
#   C5_updateprog  update_prognostic!,          Duplicated(model)
#   C6_transform   transform!,                  Duplicated(model)
#   C7_ocean       ocean_timestep!,             Duplicated(model)   [physics]
#   C8_land        land_timestep!,              Duplicated(model)   [physics]
# =============================================================================
import Pkg
Pkg.activate(@__DIR__)
using SpeedyWeather, Enzyme, Test
include("test_utils.jl")

const VARIANT = get(ARGS, 1, "A2_param")
let v = try; string(pkgversion(Enzyme)); catch; "?"; end
    println("=== BISECT variant=", VARIANT, "   Julia ", VERSION, "   Enzyme ", v, " ===")
end
flush(stdout)

const RA = set_runtime_activity(Reverse)
dup(x) = Duplicated(x, make_zero(x))

# physics ON for everything except the dynamics_only variant
const DYN_ONLY = VARIANT == "B1_dynonly"
spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
model = PrimitiveWetModel(; spectral_grid, dynamics_only = DYN_ONLY)
simulation = initialize!(model)
initialize!(simulation)
run!(simulation, period = Hour(6))          # cheap spin-up for nonzero fields
adsim = ADSimulation(simulation)
m = adsim.model
dm = make_zero(m)                            # single model shadow -> preserves aliasing

function go(variant)
    if variant == "A1_state"
        vars, dvars = ADseed(adsim, :prognostic)
        autodiff(RA, SpeedyWeather.time_step!, Const,
            Duplicated(vars, dvars), Const(m.time_stepping), Const(m))

    elseif variant == "A2_param" || variant == "B1_dynonly"
        vars, dvars = ADseed(adsim, :prognostic)
        autodiff(RA, SpeedyWeather.time_step!, Const,
            Duplicated(vars, dvars),
            Duplicated(m.time_stepping, dm.time_stepping),
            Duplicated(m, dm))

    elseif variant == "C1_dyntend"
        vars, dvars = ADseed(adsim, :tendencies)
        autodiff(RA, SpeedyWeather.dynamics_tendencies!, Const,
            Duplicated(vars, dvars), Duplicated(m, dm))

    elseif variant == "C2_paramtend"
        vars, dvars = ADseed(adsim, :tendencies)
        autodiff(RA, SpeedyWeather.parameterization_tendencies!, Const,
            Duplicated(vars, dvars), Duplicated(m, dm))

    elseif variant == "C3_implicit"
        vars, dvars = ADseed(adsim, :tendencies)
        autodiff(RA, SpeedyWeather.implicit_correction!, Const,
            Duplicated(vars, dvars),
            Duplicated(m.implicit, dm.implicit),
            Duplicated(m.time_stepping, dm.time_stepping),
            Duplicated(m, dm))

    elseif variant == "C4_diffusion"
        vars, dvars = ADseed(adsim, :tendencies)
        autodiff(RA, SpeedyWeather.horizontal_diffusion!, Const,
            Duplicated(vars, dvars), Duplicated(m, dm))

    elseif variant == "C5_updateprog"
        vars, dvars = ADseed(adsim, :prognostic)
        autodiff(RA, SpeedyWeather.update_prognostic!, Const,
            Duplicated(vars, dvars), Duplicated(m, dm))

    elseif variant == "C6_transform"
        vars, dvars = ADseed(adsim, :grid)
        autodiff(RA, SpeedyWeather.transform!, Const,
            Duplicated(vars, dvars), Duplicated(m, dm))

    elseif variant == "C7_ocean"
        vars, dvars = ADseed(adsim, :prognostic)
        autodiff(RA, SpeedyWeather.ocean_timestep!, Const,
            Duplicated(vars, dvars), Duplicated(m, dm))

    elseif variant == "C8_land"
        vars, dvars = ADseed(adsim, :prognostic)
        autodiff(RA, SpeedyWeather.land_timestep!, Const,
            Duplicated(vars, dvars), Duplicated(m, dm))

    else
        error("unknown variant: $variant")
    end
    return nothing
end

println(">>> entering autodiff for variant ", VARIANT); flush(stdout)
go(VARIANT)
println("RESULT: SUCCESS — variant ", VARIANT, " compiled + autodiff returned")
