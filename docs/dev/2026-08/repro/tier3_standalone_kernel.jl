# Tier 3 MWE: NO SpeedyWeather dependency at all. Verbatim copy of the
# callable-struct functor (SpeedyWeather/src/dynamics/initial_conditions.jl:457-491,
# JablonowskiVorticity) and the kernel it's called from
# (SpeedyWeather/src/variables/set.jl:230-234, set_field_3d_kernel!), launched
# directly with KernelAbstractions against AMDGPU. Grid sizes below mimic a
# default T31-ish octahedral Gaussian grid (npoints=3168) with 8 layers; exact
# sizes shouldn't matter for an InvalidIRError, they're just illustrative.
#
# Usage:
#   julia --project=<env-with-AMDGPU-and-KernelAbstractions> docs/dev/2026-08/repro/tier3_standalone_kernel.jl
#
# Purpose: if this crashes, it's a pure AMDGPU.jl/GPUCompiler bug triggered by
# this exact functor+kernel shape, independent of anything SpeedyWeather-specific,
# and worth reporting upstream to JuliaGPU/AMDGPU.jl directly. If it does NOT
# crash, something about how SpeedyWeather actually invokes this (Adapt-ed
# argument types, launch! sizing/workgroup choice, closures, etc.) matters --
# see tier2 for the SpeedyWeather-level equivalent.

using AMDGPU
using KernelAbstractions
using Adapt

struct JablonowskiVorticity{NF} <: Function
    sinφc::NF
    cosφc::NF
    λc::NF
    radius::NF
    u₀::NF
    η₀::NF
    perturb_uₚ::NF
    R::NF
end

Adapt.@adapt_structure JablonowskiVorticity

@inline function (J::JablonowskiVorticity)(λ, φ, η)
    (; sinφc, cosφc, λc, radius, u₀, η₀, perturb_uₚ, R) = J

    X = clamp(sinφc * sind(φ) + cosφc * cosd(φ) * cosd(λ - λc), 0, 1)
    r = radius * acos(X)

    cosηᵥ = cos((η - η₀) * π / 2)
    ζ = -4 * u₀ / radius * (cosηᵥ * sqrt(cosηᵥ)) * sind(φ) * cosd(φ) * (2 - 5sind(φ)^2)

    perturbation = perturb_uₚ / radius * exp(-(r / R)^2) *
        (tand(φ) - 2 * (radius / R)^2 * acos(X) * (sinφc * cosd(φ) - cosφc * sind(φ) * cosd(λ - λc)) / sqrt(1 - X^2))

    return ζ + perturbation
end

@kernel function set_field_3d_kernel!(var, londs, latds, σ_levels_full, f, kernel_func)
    ij, k = @index(Global, NTuple)
    var[ij, k] = kernel_func(var[ij, k], f(londs[ij], latds[ij], σ_levels_full[k]))
end

npoints, nlayers = 3168, 8

var = ROCArray(zeros(Float32, npoints, nlayers))
londs = ROCArray(Float32.(360 .* rand(npoints)))
latds = ROCArray(Float32.(180 .* rand(npoints) .- 90))
σ_levels_full = ROCArray(Float32.((1:nlayers) ./ nlayers))

radius = 6.371f6
J = JablonowskiVorticity(sind(40f0), cosd(40f0), 20f0, radius, 35f0, 0.252f0, 1f0, radius / 10)
kernel_func = (a, b) -> b   # add = false

backend = KernelAbstractions.get_backend(var)
loop! = set_field_3d_kernel!(backend)
loop!(var, londs, latds, σ_levels_full, J, kernel_func; ndrange = (npoints, nlayers))
KernelAbstractions.synchronize(backend)

println("OK, sample value: ", Array(var)[1, 1])
