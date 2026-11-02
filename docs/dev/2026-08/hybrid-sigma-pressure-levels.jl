using SpeedyWeather, CairoMakie
const SW = SpeedyWeather
nlayers = 16
sg = SpectralGrid(truncation = 31; nlayers)
σ_half = SW.sigma_half_spacing(nlayers, SW.frierson_profile)
coords = (
    "Sigma (Frierson spacing)" => SigmaCoordinates(sg, σ_half),
    "Cubic transition 0.2–0.8" => CubicSigmaPressureCoordinates(sg, σ_half; pressure_only_above = 0.2, σ_only_below = 0.8),
    "Ramp transition (new)" => RampSigmaPressureCoordinates(sg, σ_half),
)
x = range(-3, 3, 400)                           # horizontal distance [arbitrary]
pₛ = @. 1e5 - 5e4 * exp(-x^2)                   # Gaussian mountain, surface pressure down to 500 hPa

fig = Figure(size = (1200, 480), fontsize = 14)
for (i, (name, c)) in enumerate(coords)
    Δ(ps) = [SW.pressure_thickness(k, ps, c) for k in 1:nlayers]
    ratio = minimum(Δ(5e4) ./ (0.5 .* Δ(1e5)))
    ax = Axis(fig[1, i]; title = name, subtitle = "thinnest layer at 500 hPa: $(round(ratio, digits = 2)) × σ thickness",
        xlabel = "horizontal distance",
        ylabel = i == 1 ? "pressure [hPa]" : "", yreversed = true,
        xgridvisible = false, ygridvisible = false)
    i > 1 && hideydecorations!(ax; ticks = false, ticklabels = false)
    band!(ax, x, pₛ ./ 100, fill(1000, length(x)); color = (:gray60, 0.6))
    for k in 1:(nlayers + 1)
        p = [SW.pressure_half(k, ps, c) for ps in pₛ] ./ 100
        lines!(ax, x, p; color = :black, linewidth = k == nlayers + 1 ? 2 : 1)
    end
    ylims!(ax, 1010, 0)
end
save(joinpath(@__DIR__, "hybrid-sigma-pressure-levels.png"), fig, px_per_unit = 2)
