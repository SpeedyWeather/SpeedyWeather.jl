#=
Generate `docs/src/benchmarks.md` from the benchmark results JSON produced by
`SpeedyWeather/benchmark/manual_benchmarking.jl`. Called from `docs/make.jl`
before `makedocs` so the page is a regular static markdown file by the time
Documenter sees it.

The page is structured as:
  - introduction
  - cross-architecture overview table (PrimitiveWet resolution sweep)
  - one comparison figure per L value (SYPD vs T, one line per arch)
  - one section per architecture mirroring the per-arch markdown stored in JSON

If the JSON is missing or empty, a placeholder page is written telling the
reader how to populate it. This keeps `make.jl` runnable on a fresh checkout.
=#

using JSON3, CairoMakie, SpeedyWeather

const BENCHMARK_RESULTS_PATH = joinpath(pkgdir(SpeedyWeather), "benchmark", "assets", "benchmark_results.json")
const DOCS_ASSETS_DIR = joinpath(@__DIR__, "src", "assets", "benchmarks")
const DOCS_PAGE_PATH = joinpath(@__DIR__, "src", "benchmarks.md")

# Stable arch ordering — matches manual_benchmarking.jl
const ARCH_ORDER = ["cpu-arm", "cpu-x86", "gpu-nvidia", "gpu-amd", "reactant-cpu", "reactant-nvidia-gpu", "reactant-amd-gpu", "reactant-tpu"]

function sorted_arch_labels(results)
    known = filter(in(keys(results)), ARCH_ORDER)
    extra = sort!([k for k in keys(results) if !(k in ARCH_ORDER)])
    return vcat(known, extra)
end

function load_results()
    isfile(BENCHMARK_RESULTS_PATH) || return Dict{String, Any}()
    try
        return Dict{String, Any}(JSON3.read(read(BENCHMARK_RESULTS_PATH, String), Dict{String, Any}))
    catch err
        @warn "Could not parse $BENCHMARK_RESULTS_PATH; rendering empty benchmarks page" exception = err
        return Dict{String, Any}()
    end
end

function overview_nlayers(all_results, labels)
    s = Set{Int}()
    for label in labels
        ov = get(all_results[label], "overview", nothing)
        ov === nothing && continue
        for l in ov["nlayers"]
            push!(s, Int(l))
        end
    end
    return sort!(collect(s))
end

# Average resolution in degrees (as `SpectralGrid` shows it) of the grid with `nlat` rings,
# the resolution sweep runs on the default grid. Inverts `get_nlat`: grids with an odd number
# of rings have nlat_half = (nlat + 1) ÷ 2 (Equator included), otherwise nlat_half = nlat ÷ 2
function resolution_degrees(nlat::Integer)
    Grid = SpeedyWeather.DEFAULT_GRID
    nlat_half = (nlat + RingGrids.nlat_odd(Grid)) ÷ 2
    return RingGrids.average_resolution(Grid, nlat_half)
end

# two significant digits below 10°, e.g. "3.5°", "0.7°", otherwise whole degrees, e.g. "11°"
format_degrees(degrees) = degrees < 10 ? "$(round(degrees, sigdigits = 2))°" : "$(round(Int, degrees))°"

# Render one SYPD-vs-T comparison figure for a given `nlayers_target`.
# Each architecture contributes two lines: standard Legendre+FFT (LT+FFT, solid
# circles) and the matrix transform (MT, dashed diamonds). Lines for the same
# architecture share a colour. The points of the top line (the one with the most
# points, of those the fastest) are labelled with their resolution in degrees.
# Returns the path saved to, or `nothing` if no architecture had data at that L.
function write_overview_figure(all_results, labels, nlayers_target, png_path)
    fig = Figure(size = (980, 480))
    ax = Axis(
        fig[1, 1];
        xlabel = "Spectral truncation T",
        ylabel = "SYPD (simulated years per wallclock day)",
        title = "PrimitiveWetModel, L=$nlayers_target",
        yscale = log10,
    )
    palette = Makie.wong_colors()
    plotted_anything = false
    # Build legend entries manually so each series shows its actual line style
    # AND marker (scatterlines' auto-legend entry only renders as a solid line).
    legend_elems = Any[]
    legend_labels = String[]
    top_line = (xs = Int[], ys = Float64[], nlats = Int[])
    xmax, ymax = 0, 0.0
    for (j, label) in enumerate(labels)
        ov = get(all_results[label], "overview", nothing)
        ov === nothing && continue
        color = palette[mod1(j, length(palette))]
        transforms = get(ov, "spectral_transform", fill("default", length(ov["truncation"])))

        for (kind, short, marker, linestyle) in (
                ("default", "LT+FFT", :circle, :solid),
                ("matrix", "MT", :diamond, :dash),
            )
            xs = Int[]
            ys = Float64[]
            nlats = Int[]
            for i in eachindex(ov["truncation"])
                Int(ov["nlayers"][i]) == nlayers_target || continue
                String(transforms[i]) == kind || continue
                s = ov["sypd"][i]
                (s isa Number && isfinite(s) && s > 0) || continue
                push!(xs, Int(ov["truncation"][i]))
                push!(ys, Float64(s))
                push!(nlats, Int(ov["nlat"][i]))
            end
            isempty(xs) && continue
            xmax, ymax = max(xmax, maximum(xs)), max(ymax, maximum(ys))
            order = sortperm(xs)
            scatterlines!(
                ax, xs[order], ys[order];
                color = color, marker = marker, linestyle = linestyle,
                linewidth = 2, markersize = 10,
            )
            if (length(xs), sum(log, ys)) > (length(top_line.xs), sum(log, top_line.ys; init = -Inf))
                top_line = (xs = xs[order], ys = ys[order], nlats = nlats[order])
            end
            push!(
                legend_elems, [
                    LineElement(color = color, linestyle = linestyle, linewidth = 2),
                    MarkerElement(color = color, marker = marker, markersize = 10),
                ]
            )
            push!(legend_labels, "$label ($short)")
            plotted_anything = true
        end
    end
    plotted_anything || return nothing
    # resolution labels above and right of each point of the top line, where the line has
    # already descended (SYPD decreases with T); the axis is extended up and right to fit them
    for (x, y, nlat) in zip(top_line.xs, top_line.ys, top_line.nlats)
        text!(
            ax, x, y; text = format_degrees(resolution_degrees(nlat)),
            align = (:left, :bottom), offset = (4, 4), fontsize = 12, color = (:black, 0.6),
        )
    end
    ylims!(ax, nothing, max(2.5 * maximum(top_line.ys), 1.5 * ymax))
    xlims!(ax, nothing, 1.1 * xmax)
    # next to the axis, as it has many entries. `patchsize` widens the swatch so the dash pattern reads clearly.
    Legend(fig[1, 2], legend_elems, legend_labels; patchsize = (70, 15))
    save(png_path, fig)
    return png_path
end

format_sypd_cell(s) = s < 10 ? string(round(s; digits = 1)) : string(round(Int, s))
transform_short_label(s) = String(s) == "matrix" ? "MT" : "LT+FFT"

function write_overview_table(io, all_results, labels)
    rows = Tuple{Int, Int, String}[]
    for label in labels
        ov = get(all_results[label], "overview", nothing)
        ov === nothing && continue
        transforms = get(ov, "spectral_transform", fill("default", length(ov["truncation"])))
        for i in eachindex(ov["truncation"])
            r = (Int(ov["truncation"][i]), Int(ov["nlayers"][i]), String(transforms[i]))
            r in rows || push!(rows, r)
        end
    end
    # LT+FFT (default) before MT (matrix) within the same (T, L) group.
    sort!(rows; by = r -> (r[1], r[2], r[3] == "default" ? 0 : 1))

    header = "| T | L | Transform | " * join(labels, " | ") * " |"
    sep = "| --- | --- | --- | " * join(fill("---", length(labels)), " | ") * " |"
    write(io, header * "\n")
    write(io, sep * "\n")
    for (t, l, tr) in rows
        cells = String[]
        for label in labels
            ov = get(all_results[label], "overview", nothing)
            cell = "—"
            if ov !== nothing
                transforms = get(ov, "spectral_transform", fill("default", length(ov["truncation"])))
                for i in eachindex(ov["truncation"])
                    if Int(ov["truncation"][i]) == t && Int(ov["nlayers"][i]) == l && String(transforms[i]) == tr
                        s = ov["sypd"][i]
                        cell = (s isa Number && isfinite(s)) ? format_sypd_cell(s) : "—"
                        break
                    end
                end
            end
            push!(cells, cell)
        end
        write(io, "| $t | $l | $(transform_short_label(tr)) | " * join(cells, " | ") * " |\n")
    end
    write(io, "\n")
    return
end

# Suffix every `### ` and `#### ` heading line with " — <label>" so the slug
# is unique across the docs and Documenter's @ref resolver isn't confused by
# generic suite titles (e.g. "Grids" or "Models"). Also rewrites any
# single-dash table separators (`| - | - |`) into triple-dash form
# (`| --- | --- |`), which Julia's Markdown parser actually recognises as a
# table — older JSON dumps may still contain the single-dash form.
function rewrite_arch_markdown(md::AbstractString, label::AbstractString)
    io = IOBuffer()
    for line in eachline(IOBuffer(md), keep = true)
        stripped = chomp(line)
        if startswith(stripped, "### ") || startswith(stripped, "#### ")
            write(io, stripped, " — ", label, "\n")
        elseif is_table_separator(stripped)
            write(io, normalize_separator(stripped), "\n")
        else
            write(io, line)
        end
    end
    return String(take!(io))
end

is_table_separator(s::AbstractString) =
    startswith(s, "|") && !isempty(s) && all(c -> c in ('|', '-', ':', ' '), s)

# Expand each isolated `-` (between pipes/spaces/colons) to `---`.
function normalize_separator(s::AbstractString)
    return replace(s, r"(?<![-])-(?![-])" => "---")
end

function write_empty_page(path)
    open(path, "w") do io
        write(io, "# [Benchmarks](@id benchmarks_page)\n\n")
        write(io, "No benchmark results have been collected yet. To populate this page, run from `SpeedyWeather/benchmark`:\n\n")
        write(io, "```\n")
        write(io, "julia --project=. manual_benchmarking.jl                # CPU\n")
        write(io, "julia --project=. manual_benchmarking.jl gpu            # CUDA GPU\n")
        write(io, "```\n\n")
        write(io, "Results are persisted to `SpeedyWeather/benchmark/assets/benchmark_results.json` and read from there by the docs build.\n")
    end
    return
end

function generate_benchmarks_page()
    mkpath(DOCS_ASSETS_DIR)
    results = load_results()

    if isempty(results)
        @info "No benchmark JSON found; writing placeholder benchmarks page"
        write_empty_page(DOCS_PAGE_PATH)
        return
    end

    labels = sorted_arch_labels(results)

    # Render figures into docs/src/assets/benchmarks/
    figure_links = Tuple{Int, String}[]   # (L, asset path used in markdown)
    for L in overview_nlayers(results, labels)
        fname = "primitivewet_L$(L).png"
        out = write_overview_figure(results, labels, L, joinpath(DOCS_ASSETS_DIR, fname))
        out === nothing || push!(figure_links, (L, joinpath("assets", "benchmarks", fname)))
    end

    open(DOCS_PAGE_PATH, "w") do io
        write(io, "# [Benchmarks](@id benchmarks_page)\n\n")
        write(io, "Performance benchmarks for SpeedyWeather.jl across architectures. ")
        write(io, "This page is auto-generated at doc-build time from `SpeedyWeather/benchmark/assets/benchmark_results.json`, ")
        write(io, "which itself is updated by `SpeedyWeather/benchmark/manual_benchmarking.jl`. ")
        write(io, "Running the benchmark script on a new architecture adds a column to the overview table, a series to each comparison figure, and a per-architecture section to this page.\n\n")

        write(io, "All simulations are benchmarked over several seconds (wallclock time) without output. ")
        write(io, "Benchmarking excludes initialization. By default each configuration is timed up to 5 times (as long as all its timed runs together take at most 30 s), ")
        write(io, "each architecture's section states its settings. SYPD is the median and Spread shows how much the runs differ. ")
        write(io, "Timings still vary between processes and machines, so treat the numbers as rough rather than precise.\n\n")

        write(io, "## Overview: PrimitiveWet resolution across architectures\n\n")
        write(io, "Simulated years per wallclock day (SYPD) for the `PrimitiveWetModel` resolution sweep, one column per architecture. ")
        write(io, "Each (T, L) configuration is reported for both the standard Legendre transform and fast Fourier transform (LT+FFT) ")
        write(io, "and the single matrix transform (MT). ")
        write(io, "Empty cells (—) mean the architecture has either not been benchmarked yet or did not run that specific configuration. ")
        write(io, "The labels in the figures give the average grid spacing of each resolution in degrees.\n\n")

        for (L, link) in figure_links
            write(io, "![PrimitiveWet L=$L across architectures]($link)\n\n")
        end

        write_overview_table(io, results, labels)

        # Per-architecture sections — reuse the markdown blob stored in JSON.
        # Disambiguate the per-arch sub-headings so Documenter's @ref resolver
        # can't confuse them with top-level pages (e.g. the `### Grids` suite
        # heading would otherwise collide with the `# Grids` page).
        for label in labels
            record = results[label]
            meta = record["meta"]
            write(io, "## Architecture: `$label`\n\n")
            write(io, "Created for SpeedyWeather.jl v$(meta["speedyweather_version"]) on $(meta["timestamp"]).\n\n")
            write(io, "### Machine details — $label\n\n")
            write(io, meta["machine_info"])
            write(io, "\n")
            write(io, rewrite_arch_markdown(record["markdown"], label))
            write(io, "\n")
        end
    end

    @info "Generated $DOCS_PAGE_PATH ($(length(labels)) arch(s), $(length(figure_links)) figure(s))"
    return
end

generate_benchmarks_page()
