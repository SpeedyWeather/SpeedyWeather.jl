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

"""Fraction of the vorticity power sitting in the top 10% of degrees, at record `column`. A tail
that rises relative to the same-grid, old-weights case is what the plan flags as a problem."""
function tail_fraction(spectrum, column)
    (size(spectrum, 2) == 0 || column > size(spectrum, 2)) && return NaN
    power = view(spectrum, :, column)
    all(isfinite, power) || return NaN
    lmax = size(spectrum, 1)
    return sum(view(power, max(1, round(Int, 0.9lmax)):lmax)) / sum(power)
end

"""The record index closest to `years` after the start, or `nothing` if the run did not get there."""
function record_at(run, years)
    isempty(run.time) && return nothing
    target = first(run.time) + Millisecond(round(Int, years * 365.25 * 24 * 3600 * 1000))
    target > last(run.time) && return nothing
    return argmin(abs.(Dates.value.(run.time .- target)))
end

paths = sort(filter(p -> endswith(p, ".jld2"), readdir(RUNS, join = true)))
isempty(paths) && error("no .jld2 runs found in $RUNS")
runs = load_run.(paths)

# Everything is compared at a common point in the run rather than at each run's last record:
# the runs fail at different times, and a diagnostic read off just before failure says more about
# how close that run was to failing than about the quadrature.
const COMPARISON_YEARS = 2.0

println("Quadrature stability matrix, $RUNS")
println("drifts and the spectral tail are measured at year $COMPARISON_YEARS, common to all runs\n")
@printf(
    "%-28s %-6s %-6s %-6s %-5s %-9s %-12s %-7s %-11s %-11s %-11s %-11s %-9s\n",
    "case", "grid", "T", "diff", "seed", "survived", "failed at", "years", "Δ lnpₛ00", "Δ mass",
    "Δ energy", "Δ enstroph", "tail@2y"
)
for run in sort(runs, by = r -> (r.grid, r.truncation, r.diffusion_hours, r.case, get(r, :seed, 0)))
    at = record_at(run, COMPARISON_YEARS)
    window = isnothing(at) ? Colon() : (1:at)
    failure_years = isnothing(run.diverged_at) ? NaN :
        (run.diverged_at - first(run.time)).value / (1000 * 3600 * 24 * 365.25)
    @printf(
        "%-28s %-6s T%-5d %-6s %-5d %-9s %-12s %-7.2f %-+11.2e %-+11.2e %-+11.2e %-+11.2e %.3e\n",
        CASE_LABELS[run.case], replace(run.grid, "Grid" => ""), run.truncation,
        string(run.diffusion_hours) * "h", get(run, :seed, 0), run.finished ? "yes" : "NO",
        isnothing(run.diverged_at) ? "-" : Dates.format(run.diverged_at, "yyyy-mm-dd"),
        failure_years,
        drift(run.mean_lnpressure[window]), drift(run.mean_pressure[window]),
        drift(run.energy[window]), drift(run.enstrophy[window]),
        isnothing(at) ? NaN : tail_fraction(run.spectrum, at)
    )
end

# ---------------------------------------------------------------------------------------------
# figures, one per (grid, truncation, diffusion) group

for group_key in unique([(r.grid, r.truncation, r.diffusion_hours) for r in runs])
    grid, truncation, diffusion_hours = group_key
    group = sort(
        filter(r -> (r.grid, r.truncation, r.diffusion_hours) == group_key, runs),
        by = r -> (r.case, get(r, :seed, 0))
    )
    # colour by case, so ensemble members of the same case share a colour and only the
    # unperturbed member carries the legend entry
    palette = Makie.wong_colors()
    case_color = Dict(c => palette[i] for (i, c) in enumerate(sort(unique(r.case for r in group))))
    color_of(run) = (case_color[run.case], get(run, :seed, 0) == 0 ? 1.0 : 0.45)
    label_of(run) = get(run, :seed, 0) == 0 ? CASE_LABELS[run.case] : nothing

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
        for run in group
            years = [(t - first(run.time)).value / (1000 * 3600 * 24 * 365.25) for t in run.time]
            series = getfield(run, field)
            keep = findall(i -> isfinite(series[i]) && (scale !== log10 || series[i] > 0), eachindex(series))
            lines!(
                axis, years[keep], series[keep], color = color_of(run),
                label = label_of(run), linewidth = 1.5
            )
        end
        index == 1 && axislegend(axis, position = :lt, labelsize = 10, framevisible = false)
    end

    # the spectral tail: first and last vorticity spectrum of every case
    axis = Axis(
        figure[3, 1:3], xlabel = "degree l", ylabel = "vorticity power",
        xscale = log10, yscale = log10, title = "vorticity power spectrum, first (dashed) and last (solid) record"
    )
    for run in group
        get(run, :seed, 0) == 0 || continue      # one spectrum per case keeps the panel readable
        size(run.spectrum, 2) == 0 && continue
        finite = [j for j in axes(run.spectrum, 2) if all(isfinite, view(run.spectrum, :, j))]
        isempty(finite) && continue
        degrees = 1:size(run.spectrum, 1)
        for (column, style) in ((first(finite), :dash), (last(finite), :solid))
            power = view(run.spectrum, :, column)
            keep = findall(>(0), power)
            lines!(
                axis, degrees[keep], power[keep], color = case_color[run.case],
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
