# TEMPORARY debugging file for the AMDGPU `InvalidIRError` seen in
# `vertical_integration.jl` (compiling `set_field_3d_kernel!` for
# `JablonowskiVorticity`/`JablonowskiDivergence`). See
# docs/dev/2026-08/amdgpu-vertical-integration-ci-failure.md for full context.
#
# `SpeedyWeather/test/GPU/runtests.jl` currently includes ONLY this file (all
# other GPU test includes are commented out there) so buildkite runs stay fast
# and a crash in one tier doesn't prevent the others from running in the same
# job. Revert that once the crash is root-caused.
#
# All three tiers are nested inside one outer @testset on purpose: a *nested*
# @testset that errors just records the error in its parent and execution
# continues with the next sibling; only the outermost @testset throws (after
# all children have run), so we get all three results from a single CI run.
@testset "AMDGPU crash MRE" begin

    @testset "tier1: initialize! on default PrimitiveWetModel" begin
        arch = SpeedyWeather.GPU()
        spectral_grid = SpectralGrid(architecture = arch)
        model = PrimitiveWetModel(spectral_grid, dynamics_only = true)
        simulation = initialize!(model)
        @test true   # reaching here means initialize! didn't crash
    end

    @testset "tier2: set! with JablonowskiVorticity kernel" begin
        arch = SpeedyWeather.GPU()
        spectral_grid = SpectralGrid(architecture = arch)
        geometry = SpeedyWeather.Geometry(spectral_grid)
        field = zeros(Float32, spectral_grid.grid, spectral_grid.nlayers)

        # same parameter construction as ZonalWind's default initialize!, see
        # SpeedyWeather/src/dynamics/initial_conditions.jl:429-449
        perturb_lat, perturb_lon, perturb_uₚ, perturb_radius = 40.0f0, 20.0f0, 1.0f0, 0.1f0
        u₀, η₀ = 35.0f0, 0.252f0
        radius = 6.371f6
        sinφc, cosφc = sind(perturb_lat), cosd(perturb_lat)
        R = radius * perturb_radius

        J = SpeedyWeather.JablonowskiVorticity(sinφc, cosφc, perturb_lon, radius, u₀, η₀, perturb_uₚ, R)
        set!(field, J, geometry; static_func = true)
        @test true
    end

    @testset "tier3: standalone kernel, no SpeedyWeather init path" begin
        # verbatim copy of the functor (initial_conditions.jl:457-491) and the
        # kernel (set.jl:230-234), launched directly via KernelAbstractions, to
        # isolate whether this is a pure AMDGPU.jl/GPUCompiler bug independent
        # of anything else SpeedyWeather does in initialize!/set!.
        struct MRE_JablonowskiVorticity{NF} <: Function
            sinφc::NF
            cosφc::NF
            λc::NF
            radius::NF
            u₀::NF
            η₀::NF
            perturb_uₚ::NF
            R::NF
        end

        Adapt.@adapt_structure MRE_JablonowskiVorticity

        @inline function (J::MRE_JablonowskiVorticity)(λ, φ, η)
            (; sinφc, cosφc, λc, radius, u₀, η₀, perturb_uₚ, R) = J

            X = clamp(sinφc * sind(φ) + cosφc * cosd(φ) * cosd(λ - λc), 0, 1)
            r = radius * acos(X)

            cosηᵥ = cos((η - η₀) * π / 2)
            ζ = -4 * u₀ / radius * (cosηᵥ * sqrt(cosηᵥ)) * sind(φ) * cosd(φ) * (2 - 5sind(φ)^2)

            perturbation = perturb_uₚ / radius * exp(-(r / R)^2) *
                (tand(φ) - 2 * (radius / R)^2 * acos(X) * (sinφc * cosd(φ) - cosφc * sind(φ) * cosd(λ - λc)) / sqrt(1 - X^2))

            return ζ + perturbation
        end

        @kernel function mre_set_field_3d_kernel!(var, londs, latds, σ_levels_full, f, kernel_func)
            ij, k = @index(Global, NTuple)
            var[ij, k] = kernel_func(var[ij, k], f(londs[ij], latds[ij], σ_levels_full[k]))
        end

        arch = SpeedyWeather.GPU()
        npoints, nlayers = 3168, 8

        var = on_architecture(arch, zeros(Float32, npoints, nlayers))
        londs = on_architecture(arch, Float32.(360 .* rand(npoints)))
        latds = on_architecture(arch, Float32.(180 .* rand(npoints) .- 90))
        σ_levels_full = on_architecture(arch, Float32.((1:nlayers) ./ nlayers))

        radius = 6.371f6
        J = MRE_JablonowskiVorticity(sind(40.0f0), cosd(40.0f0), 20.0f0, radius, 35.0f0, 0.252f0, 1.0f0, radius / 10)
        kernel_func = (a, b) -> b   # add = false

        backend = KernelAbstractions.get_backend(var)
        loop! = mre_set_field_3d_kernel!(backend)
        loop!(var, londs, latds, σ_levels_full, J, kernel_func; ndrange = (npoints, nlayers))
        KernelAbstractions.synchronize(backend)
        @test true
    end
end
