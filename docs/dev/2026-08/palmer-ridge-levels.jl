# Level lines of sigma and hybrid sigma-pressure coordinates over the 10 km "Palmer Ridge"
# (resting-atmosphere stress test, see PR #1236): ridge of ~9.7 km crest height (10 km with
# truncation), full width at half maximum 2Δx at T64 (Δx = πR/nlat ≈ 208 km), cross section.
using SpeedyWeather, CairoMakie
const SW = SpeedyWeather

nlayers = 31                                     # L31 as in the Palmer Ridge runs
sg = SpectralGrid(truncation = 64; nlayers)
σ_half = SW.sigma_half_spacing(nlayers, SW.frierson_profile)
coords = (
    "Sigma (Frierson spacing)" => SigmaCoordinates(sg, σ_half),
    "Ramp, default (0.2, 0.3)" => RampSigmaPressureCoordinates(sg, σ_half),
    "Ramp, higher start (0.05, 0.2)" =>
        RampSigmaPressureCoordinates(sg, σ_half; pressure_only_above = 0.05, ramp_width = 0.2),
)

# ridge: Gaussian cross section, FWHM = 2Δx at T64
Δx = π * 6371 / 96                               # [km], nlat = 96 at T64
fwhm = 2Δx
x = range(-1500, 1500, 801)                      # cross-ridge distance [km]
crest = 9.7                                      # crest height [km], 9.6-9.8 km after truncation
z = @. crest * exp(-4log(2) * x^2 / fwhm^2)      # orography [km]

# surface pressure of the hydrostatic US standard atmosphere (6.5 K/km to 11 km)
standard_pressure(z) = 1013.25e2 * (1 - 6.5 * z / 288.15)^5.2559   # z in km, p in Pa
pₛ = standard_pressure.(z)
p_crest = minimum(pₛ)

fig = Figure(size = (1500, 560), fontsize = 14)
for (i, (name, c)) in enumerate(coords)
    Δ(ps) = [SW.pressure_thickness(k, ps, c) for k in 1:nlayers]
    # thinnest layer at the crest relative to its sigma-coordinate thickness (p_crest/pₛ₀ * Δp at pₛ₀)
    ratio = minimum(Δ(p_crest) ./ (p_crest / 1e5 .* Δ(1e5)))
    crossing = ratio < 0
    ax = Axis(
        fig[1, i]; title = name,
        subtitle = "thinnest layer at crest ($(round(Int, p_crest / 100)) hPa): " *
            "$(round(ratio, digits = 2)) × σ thickness" * (crossing ? " (layers cross!)" : ""),
        subtitlecolor = crossing ? :firebrick : :black,
        xlabel = "cross-ridge distance [km]",
        ylabel = i == 1 ? "pressure [hPa]" : "", yreversed = true,
        xgridvisible = false, ygridvisible = false,
    )
    i > 1 && hideydecorations!(ax; ticks = false, ticklabels = false)
    band!(ax, x, pₛ ./ 100, fill(1020, length(x)); color = (:gray60, 0.6))
    P = [SW.pressure_half(k, ps, c) / 100 for ps in pₛ, k in 1:(nlayers + 1)]  # [hPa]
    for k in 1:(nlayers + 1)
        # interfaces that cross a neighbour somewhere (negative layer thickness) in red
        crosses = (k > 1 && any(P[:, k] .< P[:, k - 1])) || (k <= nlayers && any(P[:, k] .> P[:, k + 1]))
        lines!(ax, x, P[:, k]; color = crosses ? :firebrick : :black, linewidth = k == nlayers + 1 ? 2 : 0.8)
    end
    ylims!(ax, 1020, 0)
    xlims!(ax, extrema(x)...)
end
Label(
    fig[0, :], "Palmer Ridge: $(crest) km crest, FWHM 2Δx ≈ $(round(Int, fwhm)) km at T64, L$nlayers, " *
        "standard atmosphere surface pressure"; fontsize = 16, font = :bold,
)
save(joinpath(@__DIR__, "palmer-ridge-levels.png"), fig, px_per_unit = 2)
