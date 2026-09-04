# Where does the zonal-mean (m = 0) error of an analysed product come from on HEALPix?
# Per ring: (sampled mean of the Nyquist-truncated product on the ring's nlon points)
#         - (exact ring mean of the full product, Σ_m F_m conj(G_m))
# split into the truncation loss (orders above the ring's mmax dropped) and the fold
# (k = n·nlon content landing on m = 0 when the truncated product is sampled).
using SpeedyWeather
using SpeedyWeather: SpeedyTransforms, LowerTriangularArrays, RingGrids
using SpeedyWeather.SpeedyTransforms: EqualAreaQuadrature, PerOrderQuadrature
using SpeedyWeather.LowerTriangularArrays: get_lm_range
using Printf, Random, LinearAlgebra, Logging, Statistics

const NF = Float64
const T = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 128
quiet(f) = with_logger(f, ConsoleLogger(stderr, Logging.Error))

function random_field(spectrum, slope; seed)
    Random.seed!(seed)
    a = randn(LowerTriangularMatrix{Complex{NF}}, spectrum)
    for m in 1:spectrum.mmax
        range = get_lm_range(m, spectrum.lmax - 1)
        for (offset, lm) in enumerate(range)
            l = m + offset - 2
            a[lm] *= (1 + l)^(-slope)
        end
    end
    m0 = get_lm_range(1, spectrum.lmax - 1)
    a[m0] = complex.(real.(a[m0]))
    return a
end

# per-ring Fourier coefficients F[m+1, j] (north) and south, all orders m = 0..T, no truncation
function ring_fourier(S, a)
    Λ = S.legendre_polynomials
    lmax, mmax = S.spectrum.lmax, S.spectrum.mmax
    nlat_half = S.grid.nlat_half
    Fn = zeros(Complex{NF}, mmax, nlat_half); Fs = zeros(Complex{NF}, mmax, nlat_half)
    for m in 1:mmax
        rng = get_lm_range(m, lmax - 1)
        for j in 1:nlat_half
            n = zero(Complex{NF}); s = zero(Complex{NF})
            for (off, lm) in enumerate(rng)
                v = a[lm] * Λ[lm, j]
                n += v
                s += iseven(off - 1) ? v : -v
            end
            Fn[m, j] = n; Fs[m, j] = s
        end
    end
    return Fn, Fs
end

series(F, M, φ) = real(F[1]) + 2 * sum((real(F[m + 1] * cis(m * φ)) for m in 1:M); init = 0.0)

# errors of the ring mean of the product f·g when the ring has longitudes `lond` and the fields on
# it are truncated at order M: total, truncation loss, fold
function ring_errors(F, G, lond, M)
    φ = deg2rad.(lond)
    f = [series(F, M, p) for p in φ]
    g = [series(G, M, p) for p in φ]
    sampled = mean(f .* g)
    exact = real(F[1] * G[1]) + 2 * sum(real(F[m + 1] * conj(G[m + 1])) for m in 1:(length(F) - 1))
    trunc = real(F[1] * G[1]) + 2 * sum((real(F[m + 1] * conj(G[m + 1])) for m in 1:M); init = 0.0)
    return sampled - exact, trunc - exact, sampled - trunc, exact
end

function analyse(Grid, dealiasing, slope; same = false, pad = 0)
    nlat_half = SpeedyTransforms.get_nlat_half(T, dealiasing)
    S = quiet() do
        SpectralTransform(Spectrum(T), Grid(nlat_half); NF, Quadrature = PerOrderQuadrature)
    end
    grid = S.grid
    a = random_field(S.spectrum, slope; seed = 11)
    b = same ? a : random_field(S.spectrum, slope; seed = 23)
    Fn, Fs = ring_fourier(S, a)
    Gn, Gs = ring_fourier(S, b)

    # validate against the real transform on the first belt ring (nothing truncated there)
    fa = transform(a, S)
    rings = RingGrids.eachring(grid)
    jv = findfirst(j -> S.mmax_truncation[j] == S.spectrum.mmax - 1, 1:nlat_half)
    lond = RingGrids.get_lond_per_ring(typeof(grid), nlat_half, jv)
    mine = [series(Fn[:, jv], S.mmax_truncation[jv], deg2rad(p)) for p in lond]
    theirs = fa[rings[jv]]
    @printf("  [validation vs transform, ring %d: max abs diff %.1e, field rms %.2f]\n",
        jv, maximum(abs.(mine - theirs)), sqrt(mean(theirs .^ 2)))

    nlat = RingGrids.get_nlat(grid)
    solid = RingGrids.get_solid_angles(grid)      # per ring? check length
    ring_area = length(solid) == nlat ? solid : [sum(solid[r]) for r in rings]
    w = [ring_area[j] * RingGrids.get_nlon_per_ring(grid, j) for j in 1:nlat] ./ (4π)   # ring weight
    if length(solid) == nlat
        w = [solid[j] * RingGrids.get_nlon_per_ring(grid, j) for j in 1:nlat] ./ (4π)
    end

    # per ring (north hemisphere, j = 1..nlat_half) errors; south by symmetry index nlat+1-j
    tot = zeros(nlat); trn = zeros(nlat); fld = zeros(nlat); ex = zeros(nlat)
    for j in 1:nlat_half
        nlon = RingGrids.get_nlon_per_ring(grid, j)
        lond = RingGrids.get_lond_per_ring(typeof(grid), nlat_half, j)
        M = S.mmax_truncation[j]
        if pad > 0 && nlon < 3T + 1      # hypothetical: pad this ring by `pad` longitudes
            nlon2 = nlon + pad
            lond = [lond[1] + 360 / nlon2 * (p - 1) for p in 1:nlon2]
            M = min(S.spectrum.mmax - 1, (nlon2 - 1) ÷ 2)
        end
        for (hemi, F, G) in ((j, Fn, Gn), (nlat + 1 - j, Fs, Gs))
            t, r, f, e = ring_errors(F[:, j], G[:, j], lond, M)
            tot[hemi], trn[hemi], fld[hemi], ex[hemi] = t, r, f, e
        end
    end
    return (; S, grid, nlat, nlat_half, w, tot, trn, fld, ex, Fn, Gn)
end

function report(label, R; N = (0, 1, 2, 4, 8, 16, 32))
    (; nlat, nlat_half, w, tot, trn, fld, ex) = R
    scale = sqrt(sum(w .* ex .^ 2))            # rms of the exact ring means
    println("\n== $label  (nlat_half = $nlat_half, rms exact ring mean = $(@sprintf("%.3e", scale)))")
    println("  ring   nlon   lat      total err   trunc loss    fold      |rel|")
    for j in vcat(1:8, 12, 16, 24, 32, 48, 64, 96, nlat_half)
        j > nlat_half && continue
        @printf("  %4d  %5d  %6.2f  %+.3e  %+.3e  %+.3e  %.1e\n", j,
            RingGrids.get_nlon_per_ring(typeof(R.grid), R.nlat_half, j), RingGrids.get_latd(typeof(R.grid), R.nlat_half)[j],
            tot[j], trn[j], fld[j], abs(tot[j]) / scale)
    end
    # global-mean (l = 0) error, and where it comes from: cumulative over the innermost rings
    g = sum(w .* tot)
    @printf("  global-mean error Σ w_j e_j = %+.3e  (relative to global mean of product %.3e: %.1e)\n",
        g, sum(w .* ex), abs(g / sum(w .* ex)))
    println("  zonal error left if the innermost N rings (each hemisphere) had NO error:")
    for n in N
        keep = [j for j in 1:nlat if j > n && j <= nlat - n]
        e2 = sqrt(sum(w[keep] .* tot[keep] .^ 2))
        e2all = sqrt(sum(w .* tot .^ 2))
        @printf("    N = %2d:  rms zonal-mean error %.3e  (%.1f%% of all-ring value %.3e), |Σ w e| = %.2e\n",
            n, e2, 100 * e2 / e2all, e2all, abs(sum(w[keep] .* tot[keep])))
    end
end


using SpeedyWeather.SpeedyTransforms: ring_solid_angles

"""Analyse the per-ring error as the model would: a zonally uniform grid field carrying e_j on ring j,
pushed through the HEALPix transform's own analysis; returns the m = 0 coefficients."""
function zonal_error_spectrum(R, N)
    (; S, grid, nlat, tot) = R
    field = zeros(NF, grid)
    for (j, r) in enumerate(RingGrids.eachring(grid))
        (j <= N || j > nlat - N) && continue
        field[r] .= tot[j]
    end
    spec = transform(field, S)
    return spec[get_lm_range(1, S.spectrum.lmax - 1)]
end

function truth_zonal(slope, Sref)
    a = random_field(Sref.spectrum, slope; seed = 11)
    b = random_field(Sref.spectrum, slope; seed = 23)
    p = transform(transform(a, Sref) .* transform(b, Sref), Sref)
    return p[get_lm_range(1, Sref.spectrum.lmax - 1)]
end

"""1/σmin of the m = 0 analysis block (physical metric) when the innermost N rings are excluded."""
function zonal_conditioning(R, N)
    (; S, grid, nlat_half, nlat) = R
    nlons = [RingGrids.get_nlon_per_ring(grid, j) for j in 1:nlat_half]
    g0 = ring_solid_angles(nlons, nlat, nlat_half, RingGrids.get_solid_angles(grid))
    rng = get_lm_range(1, S.spectrum.lmax - 1)
    keep = (N + 1):nlat_half
    out = Float64[]
    for parity in 0:1
        rows = [lm for (off, lm) in enumerate(rng) if (off - 1) % 2 == parity && off - 1 <= T]
        B = [S.legendre_polynomials[lm, j] * sqrt(g0[j]) for lm in rows, j in keep]
        push!(out, 1 / minimum(svdvals(B)))
    end
    return out
end

println("T = $T, Float64. Ring-mean errors of the analysed product of two band-limited fields.\n")
Sref = quiet() do
    nl = SpeedyTransforms.get_nlat_half(T, 6.0)
    SpectralTransform(Spectrum(T), FullGaussianGrid(nl); NF, Quadrature = EqualAreaQuadrature)
end
for slope in (1.0,)
    println("\n######## coefficient slope (1+l)^-$slope ########")
    tz = truth_zonal(slope, Sref)
    low = 1:11
    for (label, Grid, d) in (("HEALPixGrid d3.5", HEALPixGrid, 3.5),
                             ("OctaHEALPixGrid d3.5", OctaHEALPixGrid, 3.5),
                             ("OctahedralGaussianGrid d2", OctahedralGaussianGrid, 2.0),
                             ("OctaminimalGaussianGrid d2", OctaminimalGaussianGrid, 2.0),
                             ("OctaminimalGaussianGrid d3.5", OctaminimalGaussianGrid, 3.5),
                             ("HEALPixPaddedGrid d3.5", HEALPixPaddedGrid, 3.5))
        R = analyse(Grid, d, slope)
        report(label, R)
        println("  through the transform: relative error of the zonal-mean coefficients, all l / l <= 10,")
        println("  with the innermost N rings (per hemisphere) removed from the analysis; and 1/sigma_min of the")
        println("  m = 0 analysis block (even, odd parity) on the remaining rings:")
        for N in (0, 2, 4, 6, 8, 12, 16)
            ez = zonal_error_spectrum(R, N)
            c = zonal_conditioning(R, N)
            @printf("    N = %2d:  all l %.2e   l<=10 %.2e   1/sigma_min = %.3f / %.3f\n", N,
                norm(ez) / norm(tz), norm(ez[low]) / norm(tz[low]), c[1], c[2])
        end
    end
    println("\n-- hypothetical: same HEALPix latitudes, every ring below 3T+1 longitudes padded by +16 --")
    R = analyse(HEALPixGrid, 3.5, slope; pad = 16)
    report("HEALPixGrid d3.5, rings +16 longitudes", R)
    println("\n-- hypothetical: +8 longitudes --")
    R = analyse(HEALPixGrid, 3.5, slope; pad = 8)
    report("HEALPixGrid d3.5, rings +8 longitudes", R; N = (0, 4, 8))
    println("\n-- same field twice (u·u-like, sign of the error is systematic) --")
    R = analyse(HEALPixGrid, 3.5, slope; same = true)
    report("HEALPixGrid d3.5, a·a", R; N = (0, 4, 8))
end
