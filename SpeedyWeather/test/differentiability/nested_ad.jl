# Nested Enzyme AD (reverse-over-forward) through the spectral transform.
#
# A Hutchinson trace estimate vᵀJv is a forward-mode JVP contracted to a scalar; differentiating it
# w.r.t. a parameter is therefore a second derivative — reverse over forward.
# This test checks that nesting now COMPILES and is CORRECT.

# Top-level (non-closure) function that we differentiate
function _nested_scalar!(out, coeffs, p, spec, grid, S)
    @inbounds for k in eachindex(parent(spec))
        parent(spec)[k] = p[1] * complex(coeffs[k], zero(eltype(coeffs)))
    end
    transform!(grid, spec, S)               # spec→grid; grid is linear in p
    out[1] = sum(abs2, grid)
    return nothing
end

# inner forward-mode JVP: directional derivative of `out` along `v`, at parameter `p`
function _nested_vJv(coeffs, p, v, spec, grid, S)
    o = zero(coeffs[1:1]); Jv = zero(coeffs[1:1])
    autodiff(
        set_runtime_activity(Forward), _nested_scalar!, Const,
        Duplicated(o, Jv), Duplicated(coeffs, v), Duplicated(p, make_zero(p)),
        Duplicated(spec, make_zero(spec)), Duplicated(grid, make_zero(grid)), Const(S),
    )
    return Jv[1]
end

@testset "Differentiability: nested reverse-over-forward through transform!" begin
    NF = Float64
    spectral_grid = SpectralGrid(truncation = 9, nlayers = 1, NF = NF)
    S = SpectralTransform(spectral_grid)

    spec0 = zeros(Complex{NF}, spectral_grid.spectrum)
    grid0 = zeros(NF, spectral_grid.grid)
    nc = length(parent(spec0))
    coeffs = randn(Random.MersenneTwister(2), NF, nc)
    v = ones(NF, nc)                        # forward-mode probe direction

    # OUTER reverse pass over the parameter p of the inner JVP — the nested case that used to fail.
    # spec/grid/S are passed as explicit outer arguments (not captured).
    pB = [1.0]
    dp = zeros(1)
    autodiff(
        set_runtime_activity(Reverse), _nested_vJv, Active,
        Const(coeffs), Duplicated(pB, dp), Const(v),
        Duplicated(spec0, make_zero(spec0)), Duplicated(grid0, make_zero(grid0)), Const(S),
    )
    @test all(isfinite, dp)
    @test dp[1] != 0

    # correctness: the outer reverse derivative must match a finite difference of the inner JVP w.r.t p
    jvp(p) = _nested_vJv(coeffs, [p], v, zeros(Complex{NF}, spectral_grid.spectrum), zeros(NF, spectral_grid.grid), S)
    fd = (jvp(1.0 + 1.0e-4) - jvp(1.0 - 1.0e-4)) / 2.0e-4
    @test isapprox(dp[1], fd; rtol = 1.0e-6)
end
