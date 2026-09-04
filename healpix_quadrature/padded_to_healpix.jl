# Experiment 3 of the polar-truncation ladder: can a field integrated on `HEALPixPaddedGrid` be
# delivered on true HEALPix pixels without loss?
#
# The padded grid shares every latitude, every ring area and the belt with `HEALPixGrid`; the only
# difference is that the polar-cap ring `j` carries `4j + 16` longitudes instead of `4j`. Two rings
# at the same latitude with different equidistant longitude samplings are related by a Fourier
# resampling along the ring, which is *exact* for content below the shorter ring's Nyquist limit.
# So the question is not whether the resampling is accurate — it is exactly what a truncation to
# the target ring's Nyquist does — but how much the delivered HEALPix map differs from what a
# HEALPix-native run would have produced, and whether the difference is the discarded polar
# content (expected, and the whole point) or something worse.
#
# Three deliveries are compared, all producing the identical `HEALPixGrid` pixel set:
#
#   1. Fourier ring resampling  (this file, `resample_ring`)   -- exact up to the target Nyquist
#   2. Spectral synthesis onto HEALPixGrid                     -- the reference, what output would do
#   3. The model's own AnvilInterpolator                       -- what NetCDFOutput uses today
#
#   julia --project=healpix_quadrature healpix_quadrature/padded_to_healpix.jl [truncation]
using SpeedyWeather
using SpeedyWeather: SpeedyTransforms, RingGrids, LowerTriangularArrays
using SpeedyWeather.LowerTriangularArrays: get_lm_range
using Printf, Random, LinearAlgebra, Logging, FFTW, Statistics, Test

const NF = Float64
const T = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 128
quiet(f) = with_logger(f, ConsoleLogger(stderr, Logging.Error))

"""A random real field truncated at `T` with an atmosphere-like `(1+l)^-1` spectrum."""
function random_field(spectrum; slope = 1.0, seed = 11)
    Random.seed!(seed)
    a = randn(LowerTriangularMatrix{Complex{NF}}, spectrum)
    for m in 1:spectrum.mmax
        rng = get_lm_range(m, spectrum.lmax - 1)
        for (off, lm) in enumerate(rng)
            a[lm] *= (1 + (m + off - 2))^(-slope)
        end
    end
    m0 = get_lm_range(1, spectrum.lmax - 1)
    a[m0] = complex.(real.(a[m0]))
    return a
end

"""Resample one ring of `n_in` equidistant samples starting at longitude `lon0_in` onto `n_out`
equidistant samples starting at `lon0_out`. Exact for content below `min(n_in, n_out)/2`.

`rfft` puts its phase origin on the *first sample*, so `F[m]` is the coefficient about `lon0_in`
and the coefficient about longitude 0 is `F[m]·cis(-m·lon0_in)`. Re-expressing it about
`lon0_out` gives the shift below. Verified exact to roundoff on signals both samplings resolve,
in both the up- and down-sampling direction."""
function resample_ring(values::AbstractVector, lon0_in, n_out, lon0_out)
    n_in = length(values)
    F = FFTW.rfft(values) ./ n_in                    # coefficients about lon0_in
    nkeep = min(length(F), n_out ÷ 2 + 1)
    G = zeros(ComplexF64, n_out ÷ 2 + 1)
    for k in 1:nkeep
        m = k - 1
        G[k] = F[k] * cis(m * deg2rad(lon0_out - lon0_in))
    end
    # the Nyquist bin of an even-length output carries only the real part
    iseven(n_out) && (G[end] = complex(real(G[end])))
    return FFTW.irfft(G .* n_out, n_out)
end

"""Map a `HEALPixPaddedField` onto the `HEALPixGrid` pixel set by per-ring Fourier resampling."""
function padded_to_healpix(field_padded, grid_padded, grid_healpix)
    out = zeros(NF, grid_healpix)
    nlat = RingGrids.get_nlat(grid_healpix)
    rings_in = RingGrids.eachring(grid_padded)
    rings_out = RingGrids.eachring(grid_healpix)
    for j in 1:nlat
        vin = field_padded[rings_in[j]]
        n_out = RingGrids.get_nlon_per_ring(grid_healpix, j)
        lon0_in = first(RingGrids.get_lond_per_ring(typeof(grid_padded), grid_padded.nlat_half, j))
        lon0_out = first(RingGrids.get_lond_per_ring(typeof(grid_healpix), grid_healpix.nlat_half, j))
        out[rings_out[j]] .= resample_ring(vin, lon0_in, n_out, lon0_out)
    end
    return out
end

"""The resampler must be exact whenever both samplings resolve the signal — that is the property
the whole delivery argument rests on, so it is asserted rather than assumed."""
function selftest()
    f(x) = 2.0 + 0.7 * cosd(x - 33.0)               # m = 0 and 1 only
    for (n_in, n_out) in ((20, 4), (24, 8), (32, 16), (20, 20), (4, 20), (288, 288))
        lon0_in, lon0_out = 360 / n_in * 0.5, 360 / n_out * 0.5
        got = resample_ring(f.([360 / n_in * (i - 0.5) for i in 1:n_in]), lon0_in, n_out, lon0_out)
        want = f.([360 / n_out * (i - 0.5) for i in 1:n_out])
        @assert maximum(abs, got - want) < 1.0e-12 "resampler wrong for $n_in -> $n_out"
    end
    return println("resampler self-test passed (exact to 1e-12 up- and down-sampling)\n")
end

function main()
    selftest()
    nlat_half = SpeedyTransforms.get_nlat_half(T, 3.5)
    grid_padded = HEALPixPaddedGrid(nlat_half)
    grid_healpix = HEALPixGrid(nlat_half)
    S_padded = quiet() do
        SpectralTransform(Spectrum(T), grid_padded; NF)
    end
    S_healpix = quiet() do
        SpectralTransform(Spectrum(T), grid_healpix; NF)
    end

    println("T$T, nlat_half = $nlat_half, Float64")
    println("  padded  npoints = ", RingGrids.get_npoints(HEALPixPaddedGrid, nlat_half))
    println("  healpix npoints = ", RingGrids.get_npoints(HEALPixGrid, nlat_half))
    println()

    a = random_field(S_padded.spectrum)

    # the reference delivery: synthesise the same spectral state directly onto HEALPix pixels
    reference = transform(a, S_healpix)

    # delivery 1: synthesise on the padded grid (what the model integrates on), then resample rings
    on_padded = transform(a, S_padded)
    resampled = padded_to_healpix(on_padded, grid_padded, grid_healpix)

    # delivery 3: the model's own interpolator, padded -> healpix
    interpolated = zeros(NF, grid_healpix)
    RingGrids.interpolate!(interpolated, on_padded)

    scale = sqrt(mean(reference .^ 2))
    err(x) = sqrt(mean((x .- reference) .^ 2)) / scale

    @printf("relative rms difference from direct spectral synthesis onto HEALPixGrid:\n")
    @printf("  Fourier ring resampling   %.3e\n", err(resampled))
    @printf("  AnvilInterpolator         %.3e\n", err(interpolated))
    println()

    # where the resampling difference sits: it should be exactly the content above each cap ring's
    # own Nyquist limit, i.e. what a HEALPix-native run could never have represented anyway
    println("per-ring relative difference (Fourier resampling), cap rings:")
    println("  ring  nlon_pad  nlon_hp   rel diff")
    rings_out = RingGrids.eachring(grid_healpix)
    for j in vcat(1:8, 12, 16, 24, 32, 48)
        r = rings_out[j]
        d = sqrt(mean((resampled[r] .- reference[r]) .^ 2)) / sqrt(mean(reference[r] .^ 2))
        @printf(
            "  %4d  %8d  %7d   %.3e\n", j,
            RingGrids.get_nlon_per_ring(grid_padded, j),
            RingGrids.get_nlon_per_ring(grid_healpix, j), d
        )
    end
    println()

    # And the decisive check for output fidelity: does the delivered HEALPix map carry the same
    # *large-scale* information? Analyse both back to spectral on the HEALPix grid and compare the
    # low degrees, which is what any downstream user looks at.
    back_resampled = transform(resampled, S_healpix)
    back_reference = transform(reference, S_healpix)
    low = Int[]
    for m in 1:S_healpix.spectrum.mmax
        rng = get_lm_range(m, S_healpix.spectrum.lmax - 1)
        for (off, lm) in enumerate(rng)
            (m + off - 2) <= 20 && push!(low, lm)
        end
    end
    @printf(
        "re-analysed on HEALPixGrid, relative difference: all degrees %.3e, degree <= 20 %.3e\n",
        norm(back_resampled - back_reference) / norm(back_reference),
        norm(back_resampled[low] - back_reference[low]) / norm(back_reference[low])
    )
    return
end

main()
