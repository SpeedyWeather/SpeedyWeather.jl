# TEMPORARY debugging file, companion to `mre_amdgpu_crash.jl`. Does NOT
# replace or modify anything in mre_amdgpu_crash.jl; its tier1-3 stay as-is.
# This file adds further tiers to bisect *why* the crash survives.
#
# Why this file exists: PR #1196 (Julia 1.12) ran the *full*
# `vertical_integration.jl` test on the AMDGPU runner -- i.e. with the merged
# `@fastmath` fix for `acos`/`sqrt` in `JablonowskiVorticity`/
# `JablonowskiDivergence` (initial_conditions.jl) and `inbounds=true` on
# `set_field_3d_kernel!` (set.jl) both already applied -- and still hit the
# same `InvalidIRError` at the same `set_field_3d_kernel!`/`JablonowskiVorticity`
# call site. So the existing fix is necessary but not sufficient.
#
# Local analysis (no GPU needed) found a concrete, previously-unconsidered
# culprit:
#   - `@macroexpand @fastmath sind(1.0)` returns `sind(1.0)`, unchanged --
#     `@fastmath` does NOT rewrite `sind`/`cosd`/`tand`. They are not in
#     `Base.FastMath`'s rewrite table, unlike `sin`/`cos`/`tan`/`acos`/`sqrt`,
#     which the existing fix relies on being rewritten.
#   - `@code_lowered sind(1.0f0)` shows `sind` itself lowers to
#     `isinf(x) && throw(DomainError(x, "sind(x) is only defined for finite x."))`
#     plus an `isnan` check and a multi-branch range reduction -- i.e. `sind`
#     carries its own throw path that the existing fix cannot touch at all.
#   - `JablonowskiVorticity`/`JablonowskiDivergence` call `sind`/`cosd`/`tand`
#     several times each, entirely outside the `@fastmath` blocks that were
#     added (those only wrap the `acos`/`sqrt` sub-expressions).
#
# This file bisects the functor into isolated single-expression pieces (with
# and without `@fastmath`, where applicable) to find out which construct(s)
# actually retain a hostcall/throw path on AMDGPU, and tests a candidate fix
# (tier8) that additionally replaces `sind`/`cosd`/`tand` with
# `sin`/`cos`/`tan` composed with `deg2rad` -- which ARE fastmath-rewritable --
# inside one `@fastmath begin ... end` block.
#
# All testsets nested under one outer @testset, same reasoning as
# mre_amdgpu_crash.jl: a crash in one tier shouldn't stop the others from
# running in the same CI job.
@testset "AMDGPU crash bisection" begin

    # This is a compile-time IR error, not a runtime numerical one -- grid
    # size can't affect whether GPUCompiler proves the IR valid, so keep it
    # tiny for fast iteration.
    npoints, nlayers = 8, 2

    # One shared kernel for every tier below, exactly mirroring set.jl's
    # `set_field_3d_kernel!` (single kernel, generic over the functor `f`).
    @kernel inbounds = true function bisect_kernel!(var, londs, latds, σ_levels_full, f, kernel_func)
        ij, k = @index(Global, NTuple)
        var[ij, k] = kernel_func(var[ij, k], f(londs[ij], latds[ij], σ_levels_full[k]))
    end

    function run_bisect_tier(J)
        arch = SpeedyWeather.GPU()
        var = on_architecture(arch, zeros(Float32, npoints, nlayers))
        londs = on_architecture(arch, Float32.(360 .* rand(npoints)))
        latds = on_architecture(arch, Float32.(180 .* rand(npoints) .- 90))
        σ_levels_full = on_architecture(arch, Float32.((1:nlayers) ./ nlayers))
        kernel_func = (a, b) -> b   # add = false, same as real code

        backend = KernelAbstractions.get_backend(var)
        loop! = bisect_kernel!(backend)
        loop!(var, londs, latds, σ_levels_full, J, kernel_func; ndrange = (npoints, nlayers))
        KernelAbstractions.synchronize(backend)
        return true
    end

    # ---- tier4: current (fixed) functor verbatim, standalone --------------
    # Exactly initial_conditions.jl's JablonowskiVorticity as it stands today
    # (with @fastmath already applied to acos/sqrt), launched with no
    # SpeedyWeather model/init involved. Predicts a crash based on PR #1196's
    # hardware run; confirms whether that's reproducible standalone too, not
    # just via the full initialize!/set! path.
    #
    # Adapt is deliberately NOT used/imported here: set.jl:198-224 never calls
    # adapt() on `f` before passing it into launch!, so tier3's
    # `Adapt.@adapt_structure` was dead code for this repro shape.
    @testset "tier4: current fastmath fix, standalone" begin
        struct MRE4_JablonowskiVorticity{NF} <: Function
            sinφc::NF
            cosφc::NF
            λc::NF
            radius::NF
            u₀::NF
            η₀::NF
            perturb_uₚ::NF
            R::NF
        end

        @inline function (J::MRE4_JablonowskiVorticity)(λ, φ, η)
            (; sinφc, cosφc, λc, radius, u₀, η₀, perturb_uₚ, R) = J
            X = clamp(sinφc * sind(φ) + cosφc * cosd(φ) * cosd(λ - λc), 0, 1)
            r = @fastmath radius * acos(X)
            cosηᵥ = cos((η - η₀) * π / 2)
            ζ = @fastmath -4 * u₀ / radius * (cosηᵥ * sqrt(cosηᵥ)) * sind(φ) * cosd(φ) * (2 - 5sind(φ)^2)
            perturbation = @fastmath perturb_uₚ / radius * exp(-(r / R)^2) *
                (tand(φ) - 2 * (radius / R)^2 * acos(X) * (sinφc * cosd(φ) - cosφc * sind(φ) * cosd(λ - λc)) / sqrt(1 - X^2))
            return ζ + perturbation
        end

        radius = 6.371f6
        J = MRE4_JablonowskiVorticity(sind(40.0f0), cosd(40.0f0), 20.0f0, radius, 35.0f0, 0.252f0, 1.0f0, radius / 10)
        @test run_bisect_tier(J)
    end

    # ---- tier5a/b: isolate sind/cosd/tand alone -----------------------------
    # No acos, no sqrt at all -- tests whether sind/cosd/tand's own internal
    # throw(DomainError)/isnan branches are sufficient by themselves to
    # trigger the hostcall/InvalidIRError, independent of anything the
    # existing fix touches.
    @testset "tier5a: sind/cosd/tand only, unwrapped" begin
        struct MRE5a_DegTrig{NF} <: Function
            a::NF
            b::NF
        end
        @inline function (J::MRE5a_DegTrig)(λ, φ, η)
            (; a, b) = J
            return a * sind(φ) + b * cosd(φ) * cosd(λ) + tand(φ)
        end

        J = MRE5a_DegTrig(1.0f0, 1.0f0)
        @test run_bisect_tier(J)
    end

    @testset "tier5b: sind/cosd/tand only, candidate fix (sin/cos/tan ∘ deg2rad, fastmath)" begin
        # deg2rad is a plain multiplication (no throw path), and sin/cos/tan
        # ARE in Base.FastMath's rewrite table (unlike sind/cosd/tand), so
        # this replaces the untouchable degree-trig calls with
        # fastmath-eligible equivalents. If this tier passes where tier5a
        # fails, this is (part of) the actual fix.
        struct MRE5b_DegTrig{NF} <: Function
            a::NF
            b::NF
        end
        @inline function (J::MRE5b_DegTrig)(λ, φ, η)
            (; a, b) = J
            return @fastmath a * sin(deg2rad(φ)) + b * cos(deg2rad(φ)) * cos(deg2rad(λ)) + tan(deg2rad(φ))
        end

        J = MRE5b_DegTrig(1.0f0, 1.0f0)
        @test run_bisect_tier(J)
    end

    # ---- tier6a/b: isolate acos(clamp(...)) alone --------------------------
    # No sind/cosd/tand at all (φ, η already lie in usable ranges from the
    # harness's londs/latds/σ_levels_full, no degree-trig needed to build a
    # clamp-then-acos chain). Sanity-checks that @fastmath does what the
    # existing fix assumes it does, in isolation from everything else.
    @testset "tier6a: acos(clamp(...)) only, unwrapped" begin
        struct MRE6a_AcosOnly{NF} <: Function
            radius::NF
        end
        @inline function (J::MRE6a_AcosOnly)(λ, φ, η)
            X = clamp(η, 0.0f0, 1.0f0)
            return J.radius * acos(X)
        end

        J = MRE6a_AcosOnly(6.371f6)
        @test run_bisect_tier(J)
    end

    @testset "tier6b: acos(clamp(...)) only, fastmath" begin
        struct MRE6b_AcosOnly{NF} <: Function
            radius::NF
        end
        @inline function (J::MRE6b_AcosOnly)(λ, φ, η)
            X = clamp(η, 0.0f0, 1.0f0)
            return @fastmath J.radius * acos(X)
        end

        J = MRE6b_AcosOnly(6.371f6)
        @test run_bisect_tier(J)
    end

    # ---- tier7a/b: isolate sqrt alone --------------------------------------
    @testset "tier7a: sqrt only, unwrapped" begin
        struct MRE7a_SqrtOnly{NF} <: Function
            η₀::NF
        end
        @inline function (J::MRE7a_SqrtOnly)(λ, φ, η)
            cosηᵥ = cos((η - J.η₀) * π / 2)
            return cosηᵥ * sqrt(cosηᵥ)
        end

        J = MRE7a_SqrtOnly(0.252f0)
        @test run_bisect_tier(J)
    end

    @testset "tier7b: sqrt only, fastmath" begin
        struct MRE7b_SqrtOnly{NF} <: Function
            η₀::NF
        end
        @inline function (J::MRE7b_SqrtOnly)(λ, φ, η)
            cosηᵥ = cos((η - J.η₀) * π / 2)
            return @fastmath cosηᵥ * sqrt(cosηᵥ)
        end

        J = MRE7b_SqrtOnly(0.252f0)
        @test run_bisect_tier(J)
    end

    # ---- tier8: full functor, complete candidate fix -----------------------
    # Combines tier5b + tier6b + tier7b's fixes into the full expression: one
    # @fastmath begin...end block wrapping everything, with sind/cosd/tand
    # replaced by sin/cos/tan ∘ deg2rad so the whole block is actually
    # fastmath-rewritable. If tier4 crashes but tier8 doesn't, this is the fix
    # to port back into initial_conditions.jl for real.
    @testset "tier8: full functor, degree-trig + acos/sqrt all fastmath'd" begin
        struct MRE8_JablonowskiVorticityFixed{NF} <: Function
            sinφc::NF
            cosφc::NF
            λc::NF
            radius::NF
            u₀::NF
            η₀::NF
            perturb_uₚ::NF
            R::NF
        end
        @inline function (J::MRE8_JablonowskiVorticityFixed)(λ, φ, η)
            (; sinφc, cosφc, λc, radius, u₀, η₀, perturb_uₚ, R) = J
            @fastmath begin
                sinφ, cosφ = sin(deg2rad(φ)), cos(deg2rad(φ))
                cosΔλ = cos(deg2rad(λ - λc))
                X = clamp(sinφc * sinφ + cosφc * cosφ * cosΔλ, 0, 1)
                r = radius * acos(X)
                cosηᵥ = cos((η - η₀) * π / 2)
                ζ = -4 * u₀ / radius * (cosηᵥ * sqrt(cosηᵥ)) * sinφ * cosφ * (2 - 5sinφ^2)
                perturbation = perturb_uₚ / radius * exp(-(r / R)^2) *
                    (tan(deg2rad(φ)) - 2 * (radius / R)^2 * acos(X) * (sinφc * cosφ - cosφc * sinφ * cosΔλ) / sqrt(1 - X^2))
                return ζ + perturbation
            end
        end

        radius = 6.371f6
        J = MRE8_JablonowskiVorticityFixed(sind(40.0f0), cosd(40.0f0), 20.0f0, radius, 35.0f0, 0.252f0, 1.0f0, radius / 10)
        @test run_bisect_tier(J)
    end
end
