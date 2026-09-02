# Summarise the quadrature stability matrix produced by `run_case.jl`.
#
#   julia --project=healpix_quadrature healpix_quadrature/analyse_runs.jl [runs_directory]
#
# Prints one table row per run — did it survive, when did it fail, how far did the conserved
# quantities drift, how did the spectral tail move — and writes `stability_<grid>_T<n>.png` with
# the trajectories and the vorticity spectra.

using JLD2, Printf, Statistics, Dates
using CairoMakie

const RUNS = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "runs")

const CASE_LABELS = Dict(
    "A" => "A  dealias 3.0, equal area",
    "B" => "B  dealias 3.5, equal area",
    "C" => "C  dealias 3.5, per ring",
    "D" => "D  dealias 3.5, per order",
    "E" => "E  dealias 3.5, contractive",
)

"""Everything `run_case.jl` stored, as a NamedTuple."""
function load_run(path)
    return jldopen(path, "r") do file
        NamedTuple(Symbol(k) => file[k] for k in keys(file))
    end
end

"""Fractional drift of `series` from its first value, over the finite part of the record."""
function drift(series)
    finite = filter(isfinite, series)
    (isempty(finite) || iszero(first(finite))) && return NaN
    return last(finite) / first(finite) - 1
end

"""Ratio of the vorticity power in the top decade of degrees to the total, first and last record.
A tail that rises relative to the same-grid, old-weights case is what the plan flags as a problem."""
function tail_fraction(spectrum)
    size(spectrum, 2) == 0 && return (NaN, NaN)
    lmax = size(spectrum, 1)
    tail = max(1, round(Int, 0.9lmax)):lmax
    fraction(column) = sum(view(column, tail)) / sum(column)
    finite_columns = [j for j in axes(spectrum, 2) if all(isfinite, view(spectrum, :, j))]
    isempty(finite_columns) && return (NaN, NaN)
    return fraction(view(spectrum, :, first(finite_columns))),
        fraction(view(spectrum, :, last(finite_columns)))
end

paths = sort(filter(p -> endswith(p, ".jld2"), readdir(RUNS, join = true)))
isempty(paths) && error("no .jld2 runs found in $RUNS")
runs = load_run.(paths)

println("Quadrature stability matrix, $RUNS\n")
@printf(
    "%-28s %-6s %-6s %-7s %-9s %-12s %-11s %-11s %-11s %-11s %-9s\n",
    "case", "grid", "T", "diff", "survived", "failed at", "Δ lnpₛ00", "Δ mass", "Δ energy",
    "Δ enstroph", "tail"
)
for run in sort(runs, by = r -> (r.grid, r.truncation, r.diffusion_hours, r.case))
    tail_first, tail_last = tail_fraction(run.spectrum)
    @printf(
        "%-28s %-6s T%-5d %-7s %-9s %-12s %-+11.2e %-+11.2e %-+11.2e %-+11.2e %.2e→%.2e\n",
        CASE_LABELS[run.case], replace(run.grid, "Grid" => ""), run.truncation,
        string(run.diffusion_hours) * "h", run.finished ? "yes" : "NO",
        isnothing(run.diverged_at) ? "-" : Dates.format(run.diverged_at, "yyyy-mm-dd"),
        drift(run.mean_lnpressure), drift(run.mean_pressure), drift(run.energy),
        drift(run.enstrophy), tail_first, tail_last
    )
end

# ---------------------------------------------------------------------------------------------
# figures, one per (grid, truncation, diffusion) group

for group_key in unique([(r.grid, r.truncation, r.diffusion_hours) for r in runs])
    grid, truncation, diffusion_hours = group_key
    group = sort(
        filter(r -> (r.grid, r.truncation, r.diffusion_hours) == group_key, runs), by = r -> r.case
    )
    colors = Makie.wong_colors()

    figure = Figure(size = (1500, 900))
    Label(
        figure[0, 1:3],
        "$grid T$truncation, hyperdiffusion $(diffusion_hours) h — quadrature schemes compared",
        fontsize = 22, font = :bold
    )

    panels = (
        (:max_vorticity, "max |ζ| [1/s]", log10),
        (:enstrophy, "enstrophy [1/s²]", log10),
        (:energy, "kinetic energy [m²/s²]", identity),
        (:angular_momentum, "relative angular momentum", identity),
        (:mean_pressure, "global mean pₛ [Pa]", identity),
        (:max_vertical_velocity, "max |w|", log10),
    )

    for (index, (field, label, scale)) in enumerate(panels)
        row, column = divrem(index - 1, 3)
        axis = Axis(
            figure[row + 1, column + 1], xlabel = "years", ylabel = label,
            yscale = scale === log10 ? log10 : identity
        )
        for (color_index, run) in enumerate(group)
            years = [(t - first(run.time)).value / (1000 * 3600 * 24 * 365.25) for t in run.time]
            series = getfield(run, field)
            keep = findall(i -> isfinite(series[i]) && (scale !== log10 || series[i] > 0), eachindex(series))
            lines!(
                axis, years[keep], series[keep], color = colors[color_index],
                label = CASE_LABELS[run.case], linewidth = 1.5
            )
        end
        index == 1 && axislegend(axis, position = :lt, labelsize = 10, framevisible = false)
    end

    # the spectral tail: first and last vorticity spectrum of every case
    axis = Axis(
        figure[3, 1:3], xlabel = "degree l", ylabel = "vorticity power",
        xscale = log10, yscale = log10, title = "vorticity power spectrum, first (dashed) and last (solid) record"
    )
    for (color_index, run) in enumerate(group)
        size(run.spectrum, 2) == 0 && continue
        finite = [j for j in axes(run.spectrum, 2) if all(isfinite, view(run.spectrum, :, j))]
        isempty(finite) && continue
        degrees = 1:size(run.spectrum, 1)
        for (column, style) in ((first(finite), :dash), (last(finite), :solid))
            power = view(run.spectrum, :, column)
            keep = findall(>(0), power)
            lines!(
                axis, degrees[keep], power[keep], color = colors[color_index],
                linestyle = style, linewidth = 1.5,
                label = style === :solid ? CASE_LABELS[run.case] : nothing
            )
        end
    end
    axislegend(axis, position = :lb, labelsize = 10, framevisible = false, merge = true)

    filename = joinpath(RUNS, "stability_$(grid)_T$(truncation)_diff$(diffusion_hours)h.png")
    save(filename, figure)
    println("\nwrote $filename")
end
