# Tests for SpeedyTransforms.ScaledLegendre only (step 1 of the recompute-on-the-fly Legendre
# polynomials feature, see docs/dev/2026-09/recompute-legendre-polynomials.md). No
# `SpectralTransform` integration here yet, that is a separate step; these tests exercise the
# module's own math against AssociatedLegendrePolynomials.jl's Float64 reference.

const SL = SpeedyTransforms.ScaledLegendre
const Legendre = SpeedyTransforms.Legendre

# equally spaced colatitudes standing in for a Gaussian ring set (as in the numerical prototype),
# always includes a ring close to the pole via the extra argument below
gaussian_colat(nlat_half) = [(j - 0.5) * π / (2nlat_half) for j in 1:nlat_half]

"Float64 reference lower-triangular table of λ_l^m(cos(colat[j])), via AssociatedLegendrePolynomials."
function reference_legendre(lmax, mmax, cos_colat)
    nlat_half = length(cos_colat)
    n = SL.count_nonzeros(lmax, mmax)
    Λ = zeros(Float64, n, nlat_half)
    Λj = zeros(Float64, lmax, mmax)
    for j in 1:nlat_half
        Legendre.λlm!(Λj, lmax - 1, mmax - 1, cos_colat[j])
        for m in 1:mmax, l in m:lmax
            Λ[SL.lm2i(l, m, lmax), j] = Λj[l, m]
        end
    end
    return Λ
end

"Recompute the full triangle with ScaledLegendre in number format NF, for comparison against reference_legendre."
function recomputed_legendre(::Type{NF}, lmax, mmax, cos_colat) where {NF}
    nlat_half = length(cos_colat)
    n = SL.count_nonzeros(lmax, mmax)
    αhi, αlo, βhi, βlo = SL.recursion_coefficients(NF, lmax, mmax)
    sechi, seclo, scale = SL.sectoral_modes(NF, cos_colat, mmax)
    Λ = zeros(NF, n, nlat_half)
    col = zeros(NF, lmax)
    for j in 1:nlat_half, m in 1:mmax
        ncolumn = lmax - m + 1
        lm_offset = SL.lm2i(m, m, lmax) - 1
        SL.legendre_column!(
            view(col, 1:ncolumn), cos_colat[j], αhi, αlo, βhi, βlo, lm_offset, ncolumn,
            sechi[j, m], seclo[j, m], scale[j, m],
        )
        Λ[(lm_offset + 1):(lm_offset + ncolumn), j] .= col[1:ncolumn]
    end
    return Λ
end

@testset "ScaledLegendre: coefficient identities" begin
    # the one-term sectoral -> first-off-diagonal step is the same two-term recursion with
    # α(m+1, m) == ν(m+1) and β(m+1, m) == 0 (0-based l, m); this is what makes a single
    # branch-free recursion valid for the whole column, see legendre_column!
    for m in 0:60
        l = m + 1
        @test SL.coeff_α(Float64, l, m) ≈ SL.coeff_ν(Float64, l)
        @test SL.coeff_β(Float64, l, m) == 0
    end

    # β(l,m) ≈ α(l,m) / α(l-1,m) for l > m+1 (noted in the design doc as a possible future
    # simplification, verified numerically here)
    for m in 0:10:60, l in (m + 2):5:(m + 80)
        αl = SL.coeff_α(Float64, l, m)
        αlm1 = SL.coeff_α(Float64, l - 1, m)
        βl = SL.coeff_β(Float64, l, m)
        @test βl ≈ αl / αlm1 rtol = 1e-10
    end
end

@testset "ScaledLegendre: legendre_column! vs AssociatedLegendrePolynomials" begin
    for trunc in (31, 63)
        lmax = mmax = trunc + 1
        nlat_half = trunc + 1
        colat = vcat(gaussian_colat(nlat_half), deg2rad(0.05))  # + one ring very close to the pole
        cos_colat = cos.(colat)

        ref = reference_legendre(lmax, mmax, cos_colat)

        for NF in (Float32, Float64)
            rec = recomputed_legendre(NF, lmax, mmax, cos_colat)
            err = maximum(abs.(Float64.(rec) .- ref))
            @test isfinite(err)
            if NF == Float32
                @test err < 2e-6
            else
                # `scale_shift(Float64) == 61` puts the flush-to-zero threshold at
                # 2^-61 ≈ 4.3e-19, well below eps(Float64), so the error floor here is the
                # recursion's own roundoff (which grows like eps*lmax^2 near the poles),
                # not the extended-exponent scaling.
                @test err < 1e-12
            end
        end
    end
end

@testset "ScaledLegendre: extreme underflow near the pole" begin
    # a ring extremely close to the pole (~0.05°) at high order m: the sectoral mode underflows
    # long before the recursion has produced anything, so this exercises the extended-exponent
    # scaling's ability to carry the recursion through the underflow region without flushing to
    # zero and without ever overflowing/NaN-ing along the way
    trunc = 255
    lmax = mmax = trunc + 1
    cos_colat = [cos(deg2rad(0.05))]

    ref = reference_legendre(lmax, mmax, cos_colat)

    for NF in (Float32, Float64)
        rec = recomputed_legendre(NF, lmax, mmax, cos_colat)
        @test all(isfinite, rec)
        err = maximum(abs.(Float64.(rec) .- ref))
        @test isfinite(err)
        @test err < (NF == Float32 ? 2e-6 : 3e-10)
    end
end
