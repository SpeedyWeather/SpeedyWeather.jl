# Tier 2 MWE: isolates the exact failing kernel launch (`set!` on a 3D Field
# with a `JablonowskiVorticity`/`JablonowskiDivergence` functor) without building
# a full model, running `initialize!`'s ocean/land/greenhouse-gas machinery, or
# touching orography/pressure/temperature ICs. Mirrors the pattern already used
# in `SpeedyWeather/test/GPU/set.jl`, just with the Jablonowski functor instead
# of `zfac`/an anonymous function.
#
# Usage:
#   julia --project=SpeedyWeather docs/dev/2026-08/repro/tier2_set_kernel.jl
#
# If this crashes the same way as tier1, the full model/initial-conditions
# stack is irrelevant -- the bug is purely in `set!`'s kernel launch path with
# this specific functor. If this does NOT crash but tier1 does, something else
# in the `initialize!` call chain (ocean/land init, greenhouse gases, orography,
# JablonowskiTemperature, ConstantRelativeHumidity, ...) is implicated instead.

using AMDGPU
using SpeedyWeather

arch = SpeedyWeather.GPU()
spectral_grid = SpectralGrid(architecture = arch)
geometry = SpeedyWeather.Geometry(spectral_grid)

field = zeros(Float32, spectral_grid.grid, spectral_grid.nlayers)

# same parameter construction as ZonalWind's default initialize!, see
# SpeedyWeather/src/dynamics/initial_conditions.jl:429-449
perturb_lat, perturb_lon, perturb_uₚ, perturb_radius = 40f0, 20f0, 1f0, 0.1f0
u₀, η₀ = 35f0, 0.252f0
radius = 6.371f6
sinφc, cosφc = sind(perturb_lat), cosd(perturb_lat)
R = radius * perturb_radius

J = SpeedyWeather.JablonowskiVorticity(sinφc, cosφc, perturb_lon, radius, u₀, η₀, perturb_uₚ, R)

set!(field, J, geometry; static_func = true)

println("OK")
