# Tier 1 MWE: crashes inside `initialize!(model)` alone, no `run!`, no CPU
# comparison, no vertical-integration checks. This is the entire
# `vertical_integration.jl` failure boiled down to the one call that trips it.
#
# Usage:
#   julia --project=SpeedyWeather docs/dev/2026-08/repro/tier1_initialize_model.jl
#
# Expected on a healthy toolchain: prints "OK" and exits 0.
# Expected on the failing runner: InvalidIRError while compiling
# `set_field_3d_kernel!` (reported by CI as `gpu_set_field_3d_kernel!`).

using AMDGPU
using SpeedyWeather

arch = SpeedyWeather.GPU()
spectral_grid = SpectralGrid(architecture = arch)
model = PrimitiveWetModel(spectral_grid, dynamics_only = true)

# crash happens inside initialize! -> initialize!(vars, model.initial_conditions, model)
# -> initialize!(vars, ::ZonalWind, model) -> set!(vars, model; vorticity=JablonowskiVorticity(...), ...)
simulation = initialize!(model)

println("OK")
