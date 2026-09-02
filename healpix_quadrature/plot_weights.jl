# The fitted quadrature weights `g^(m)_j` as a function of order `m`, one panel per scheme.
#
#   julia --project=healpix_quadrature healpix_quadrature/plot_weights.jl
#
# Each order is fitted independently, so nothing forces `g^(m)` and `g^(m+1)` to be close. Nonlinear
# terms couple orders, so a ragged weight spectrum could imprint a systematic pattern — this is the
# "worth plotting before running anything long" check in the plan's risk list. Also reports how
# much the ring totals `Σ_j g_j` depart from 4π at each order, which per-order weights no longer
# constrain except at m = 0.
using SpeedyWeather, SpeedyWeather.SpeedyTransforms
using CairoMakie, Printf, Statistics, Logging

const ST = SpeedyTransforms
const GRIDS = (HEALPixGrid, OctaHEALPixGrid)
const QUADRATURES = (
    ST.EqualAreaQuadrature, ST.RingQuadrature, ST.PerOrderQuadrature, ST.ContractiveQuadrature,
)
const LABELS = Dict(
    ST.EqualAreaQuadrature => "equal area",
    ST.RingQuadrature => "per ring",
    ST.PerOrderQuadrature => "per order",
    ST.ContractiveQuadrature => "contractive",
)

"""The weights `ΔΩ[m, j]` a transform ended up with, recovered from the fused
`ΔΩ conj(lon_offset)` (the rotation has unit modulus), together with the ring totals per order."""
function weights(Grid, truncation, Quadrature)
    spectral_grid = SpectralGrid(; truncation, nlayers = 1, Grid)
    transform = with_logger(ConsoleLogger(stderr, Logging.Error)) do
        SpectralTransform(spectral_grid; Quadrature)
    end
    ΔΩ = abs.(Array(transform.solid_angles_rotated))
    (; nlat, nlons) = transform
    nlat_half = transform.grid.nlat_half
    # `get_solid_angles` covers all `nlat` rings; the weights only span one hemisphere
    geometric = RingGrids.get_solid_angles(transform.grid)[1:nlat_half]
    ring_total = [
        sum(((nlat - j + 1 == j) ? 1 : 2) * nlons[j] * Float64(ΔΩ[m, j]) for j in 1:nlat_half)
            for m in axes(ΔΩ, 1)
    ]
    return ΔΩ ./ geometric', ring_total, spectral_grid
end

const TRUNCATION = 64

println("Ring totals Σ_j g_j per order m, relative to 4π (only m = 0 is constrained)\n")
@printf("%-18s %-14s %-13s %-13s %-13s %-13s\n",
        "grid", "quadrature", "m=0", "max |dev|", "min ratio", "max ratio")

figure = Figure(size = (420 * length(QUADRATURES) + 60, 340 * length(GRIDS) + 120))
Label(
    figure[0, 1:length(QUADRATURES)],
    "Quadrature weights relative to equal area, g^(m)_j / ΔΩ_j, at T$TRUNCATION\n" *
        "one line per ring j, against order m",
    fontsize = 20, font = :bold
)

for (row, Grid) in enumerate(GRIDS), (column, Quadrature) in enumerate(QUADRATURES)
    ratio, ring_total, spectral_grid = weights(Grid, TRUNCATION, Quadrature)
    deviation = ring_total ./ 4π .- 1
    @printf(
        "%-18s %-14s %-+13.2e %-13.2e %-13.4f %-13.4f\n",
        nameof(Grid), LABELS[Quadrature], deviation[1], maximum(abs, deviation),
        minimum(ratio), maximum(ratio)
    )

    axis = Axis(
        figure[row, column], xlabel = "order m", ylabel = "g^(m)_j / ΔΩ_j",
        title = "$(nameof(Grid)) — $(LABELS[Quadrature])"
    )
    for j in axes(ratio, 2)
        # rings the shortcut drops at high m keep the geometric weight; only plot where they live
        m_max = spectral_grid.spectrum.mmax
        lines!(axis, 0:(m_max - 1), view(ratio, 1:m_max, j), linewidth = 0.7, color = (:steelblue, 0.55))
    end
    hlines!(axis, [1.0], color = :black, linestyle = :dash, linewidth = 1)
    ylims!(axis, 0.75, 1.25)
end

filename = joinpath(@__DIR__, "weights_vs_order.png")
save(filename, figure)
println("\nwrote $filename")
