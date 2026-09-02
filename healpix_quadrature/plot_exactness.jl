# Transform exactness of the HEALPix grids across truncation and dealiasing, for each quadrature
# scheme: equal area (the status quo), per-ring (classical HEALPix weights, the healpy/ducc/cuHPX
# baseline), per-order (exact on the retained band) and contractive (inexact by construction, but
# non-expansive and fitted to be the most accurate weighting that is). One figure per number format.
#
#   julia --project=docs plot_exactness.jl
#
# Writes exactness_<NF>.png and analysis_norm.png.
using SpeedyTransforms, RingGrids, LowerTriangularArrays
using SpeedyTransforms: EqualAreaQuadrature, RingQuadrature, PerOrderQuadrature, ContractiveQuadrature
using CairoMakie, Printf, Logging, LinearAlgebra, Random

const DEALIASINGS = (1.5, 2.0, 2.5, 3.0, 3.5)
const TRUNCATIONS = (32, 64, 128, 256)
const GRIDS = (HEALPixGrid, OctaHEALPixGrid)
const QUADRATURES = (EqualAreaQuadrature, RingQuadrature, PerOrderQuadrature, ContractiveQuadrature)
const NUMBER_FORMATS = (Float64, Float32)

const QUADRATURE_LABELS = Dict(
    EqualAreaQuadrature => "equal area (status quo)",
    RingQuadrature => "per ring (classical HEALPix / cuHPX)",
    PerOrderQuadrature => "per order (exact on the band)",
    ContractiveQuadrature => "contractive (fitted, ‖A‖ ≤ 1)",
)

# Level at which the round trip is limited by the arithmetic rather than by the quadrature, and the
# y-range to plot. Both are empirical: an exact transform lands at ~5·10⁻¹⁶ in Float64 and ~2·10⁻⁷
# in Float32 in this metric, the latter set by accumulating over ~nlat rings.
const ROUNDOFF = Dict(Float64 => 1.0e-14, Float32 => 5.0e-7)
const YLIMS = Dict(Float64 => (1.0e-17, 3.0), Float32 => (1.0e-9, 3.0))

"""Round-trip error of the spectral transform: spectral → grid → spectral applied to a random
band-limited field, as a relative L2 and L∞ over the coefficients. Also returns the round-trip
energy gain `‖a2‖/‖a‖` — the number that decides stability, as opposed to the error norms which
only say how far from the identity the round trip is — and how many orders `m` no ring carries at
all, which no weighting can repair."""
function roundtrip_error(Grid, truncation, dealiasing, Quadrature, NF; seed = 7)
    nlat_half = SpeedyTransforms.get_nlat_half(truncation, dealiasing)
    S = with_logger(ConsoleLogger(stderr, Logging.Error)) do
        SpectralTransform(Spectrum(truncation), Grid(nlat_half); NF, Quadrature)
    end
    Random.seed!(seed)
    a = randn(LowerTriangularMatrix{Complex{NF}}, S.spectrum)
    m0 = LowerTriangularArrays.get_lm_range(1, S.spectrum.lmax - 1)
    a[m0] = complex.(real.(a[m0]))              # a real field has real m = 0 coefficients
    a2 = transform(transform(a, S), S)
    dead = count(m -> all(j -> S.mmax_truncation[j] + 1 < m, 1:nlat_half), 1:S.spectrum.mmax)
    return (;
        nlat_half, dead, gain = norm(a2) / norm(a),
        L2 = norm(a2 - a) / norm(a), L∞ = maximum(abs, a2 - a) / maximum(abs, a),
    )
end

results = Dict()
println("Round-trip error, spectral → grid → spectral of a random band-limited field\n")
@printf(
    "%-9s %-18s %-30s %-7s %-9s %-8s %-7s %-11s %-11s %-11s %s\n",
    "NF", "grid", "quadrature", "trunc", "dealias", "nlat½", "slack", "L2 rel", "L∞ rel", "gain", "dead m"
)
for NF in NUMBER_FORMATS, Grid in GRIDS, Q in QUADRATURES
    for truncation in TRUNCATIONS, dealiasing in DEALIASINGS
        r = roundtrip_error(Grid, truncation, dealiasing, Q, NF)
        results[(NF, Grid, Q, truncation, dealiasing)] = r
        @printf(
            "%-9s %-18s %-30s T%-6d %-9.1f %-8d %-+7d %-11.3e %-11.3e %-11.8f %d\n",
            NF, nameof(Grid), QUADRATURE_LABELS[Q], truncation, dealiasing, r.nlat_half,
            r.nlat_half - truncation, r.L2, r.L∞, r.gain, r.dead
        )
    end
end

function make_figure(NF)
    nrow, ncol = length(GRIDS), length(QUADRATURES)
    fig = Figure(size = (380 * ncol + 60, 340 * nrow + 150))
    Label(
        fig[1, 1:ncol], "HEALPix spectral transform: round-trip error vs truncation\n" *
            "$NF, spectral → grid → spectral of a random band-limited field",
        fontsize = 22, font = :bold
    )

    colors = cgrad(:viridis, length(DEALIASINGS), categorical = true)
    # distinct markers so exact curves, which collapse onto the roundoff floor, stay distinguishable
    markers = [:circle, :rect, :utriangle, :diamond, :xcross, :star5]
    lines_for_legend = Any[]

    for (row, Grid) in enumerate(GRIDS), (col, Q) in enumerate(QUADRATURES)
        ax = Axis(
            fig[row + 1, col],
            title = "$(nameof(Grid)) — $(QUADRATURE_LABELS[Q])",
            titlesize = 14,
            xlabel = row == nrow ? "truncation" : "",
            ylabel = col == 1 ? "relative L2 round-trip error" : "",
            xscale = log2, yscale = log10,
            xticks = (collect(TRUNCATIONS), ["T$t" for t in TRUNCATIONS]),
        )

        for (i, dealiasing) in enumerate(DEALIASINGS)
            errs = [
                max(results[(NF, Grid, Q, t, dealiasing)].L2, YLIMS[NF][1])
                    for t in TRUNCATIONS
            ]
            l = scatterlines!(
                ax, collect(TRUNCATIONS), errs,
                color = colors[i], marker = markers[i], markersize = 10, linewidth = 2
            )
            row == 1 && col == 1 && push!(lines_for_legend, l)
        end

        # anything at this level is exact; the rest is quadrature error
        hlines!(ax, [ROUNDOFF[NF]], color = (:black, 0.45), linestyle = :dash)
        text!(
            ax, 0.03, 0.03, space = :relative, text = "$NF roundoff",
            fontsize = 10, color = (:black, 0.6)
        )
        ylims!(ax, YLIMS[NF]...)
    end

    Legend(
        fig[nrow + 2, 1:ncol], lines_for_legend,
        [@sprintf("dealiasing %.1f", d) for d in DEALIASINGS],
        orientation = :horizontal, framevisible = false, labelsize = 14, nbanks = 1
    )

    path = joinpath(@__DIR__, "exactness_$NF.png")
    save(path, fig)
    return path
end

for NF in NUMBER_FORMATS
    println("\nwrote ", make_figure(NF))
end


# ---------------------------------------------------------------------------------------------
# The round-trip error above is measured on a band-limited field, which is the one input for which
# `PerOrderQuadrature` is exact by construction. A nonlinear model never sends analysis a
# band-limited field, so the norm that decides whether the transform can inject energy is that of
# the analysis operator alone, over ALL grid fields:
#
#   A = Λ diag(g) diag(g⁰)^(-1/2)     per order m, per parity of l-m
#
# with Λ[l, j] = λ_lm(μ_j), g the quadrature weights in use and g⁰ the grid's geometric ring areas
# (the physical L2 metric on grid space). σmax(A) ≤ 1 means no grid field of any wavenumber can be
# analysed into more spectral energy than it had.

"""Largest singular value of the analysis operator in the physical grid metric, over all orders
`m` and both parities of `l - m`. Exactly 1 for a grid with an exact quadrature rule."""
function analysis_norm(Grid, truncation, dealiasing, Quadrature)
    nlat_half = SpeedyTransforms.get_nlat_half(truncation, dealiasing)
    S = with_logger(ConsoleLogger(stderr, Logging.Error)) do
        SpectralTransform(Spectrum(truncation), Grid(nlat_half); NF = Float64, Quadrature)
    end
    (; nlat, nlons, mmax_truncation) = S
    (; lmax, mmax) = S.spectrum
    λ = Array(S.legendre_polynomials.data)
    ΔΩ = abs.(Array(S.solid_angles_rotated))
    geometric = RingGrids.get_solid_angles(S.grid)
    multiplicity(j) = (nlat - j + 1 == j) ? 1 : 2
    g⁰ = [multiplicity(j) * nlons[j] * Float64(geometric[j]) for j in 1:nlat_half]

    σmax = 0.0
    for m in 1:mmax
        rings = [j for j in 1:nlat_half if mmax_truncation[j] + 1 >= m]
        isempty(rings) && continue
        Λ = Float64.(λ[LowerTriangularArrays.get_lm_range(m, lmax - 1), rings])
        g = [multiplicity(j) * nlons[j] * Float64(ΔΩ[m, j]) for j in rings]
        A = Λ * Diagonal(g ./ sqrt.(g⁰[rings]))
        for parity in (0, 1)
            block = @view A[(1 + parity):2:end, :]
            size(block, 1) == 0 && continue
            σmax = max(σmax, opnorm(block))
        end
    end
    return σmax
end

function make_norm_figure()
    nrow, ncol = length(GRIDS), length(QUADRATURES)
    fig = Figure(size = (380 * ncol + 60, 340 * nrow + 150))
    Label(
        fig[1, 1:ncol], "HEALPix spectral transform: norm of the analysis operator\n" *
            "‖A‖ > 1 means aliased grid content can be transformed into MORE spectral energy",
        fontsize = 22, font = :bold
    )
    colors = cgrad(:viridis, length(DEALIASINGS), categorical = true)
    markers = [:circle, :rect, :utriangle, :diamond, :xcross, :star5]
    lines_for_legend = Any[]

    for (row, Grid) in enumerate(GRIDS), (col, Q) in enumerate(QUADRATURES)
        ax = Axis(
            fig[row + 1, col],
            title = "$(nameof(Grid)) — $(QUADRATURE_LABELS[Q])", titlesize = 14,
            xlabel = row == nrow ? "truncation" : "",
            ylabel = col == 1 ? "‖analysis‖ (physical grid metric)" : "",
            xscale = log2, xticks = (collect(TRUNCATIONS), ["T$t" for t in TRUNCATIONS]),
        )
        for (i, dealiasing) in enumerate(DEALIASINGS)
            ns = [analysis_norm(Grid, t, dealiasing, Q) for t in TRUNCATIONS]
            l = scatterlines!(
                ax, collect(TRUNCATIONS), ns,
                color = colors[i], marker = markers[i], markersize = 10, linewidth = 2
            )
            row == 1 && col == 1 && push!(lines_for_legend, l)
        end
        # the stability threshold: at or below this line the transform cannot inject energy
        hlines!(ax, [1.0], color = (:red, 0.7), linestyle = :dash, linewidth = 2)
        text!(
            ax, 0.03, 0.06, space = :relative, text = "‖A‖ = 1, no energy injected",
            fontsize = 10, color = (:red, 0.8)
        )
        ylims!(ax, 0.995, 1.09)
    end

    Legend(
        fig[nrow + 2, 1:ncol], lines_for_legend,
        [@sprintf("dealiasing %.1f", d) for d in DEALIASINGS],
        orientation = :horizontal, framevisible = false, labelsize = 14, nbanks = 1
    )
    path = joinpath(@__DIR__, "analysis_norm.png")
    save(path, fig)
    return path
end

println("\nwrote ", make_norm_figure())
