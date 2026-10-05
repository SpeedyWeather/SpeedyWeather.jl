const USAGE = """
Helpers for the benchmark-regression skill, operating on the JSON files written by

    julia manual_benchmarking.jl --debug --output=FILE

Usage: julia regression.jl <command> [options]

    benchmarked-commit [--arch cpu|gpu|amdgpu] [--arch-label LABEL] [--ref origin/main] [--repo .]
        first-parent commit that introduced the currently stored results of an architecture
    stamp FILE --tree TREE --revision SHA
        verify that FILE was produced by the packages in TREE and record the revision
    table --result NAME=FILE[,FILE...] ... [--candidate main] [--threshold 0.85]
        markdown table of SYPD per revision, geometric-mean ratios and a verdict per reference
    verdict --reference FILE[,FILE...] --candidate FILE[,FILE...] --cutoff X [--group all|default|matrix]
        good (exit 0) or bad (exit 1) classification against the reference, for git bisect

Each revision may be measured several times (comma-separated files); the best (max) SYPD per
configuration is used, as noise from other processes only ever slows a run down.
Errors exit with 2, which bisect_step.sh maps to "skip".
"""

# own environment next to this script, instantiated on first use
if !isfile(joinpath(@__DIR__, "Manifest.toml"))
    import Pkg
    Pkg.activate(@__DIR__; io = devnull)
    Pkg.instantiate(; io = devnull)
end
pushfirst!(LOAD_PATH, @__DIR__)

using JSON3, Printf

const RESULTS_JSON = "SpeedyWeather/benchmark/assets/benchmark_results.json"
const GROUPS = ("all", "default", "matrix")
const GROUP_LABELS = Dict("all" => "all", "default" => "LT+FFT", "matrix" => "MT")

struct UsageError <: Exception
    msg::String
end

# `--key value` or `--key=value` options (repeatable) and positional arguments
function parse_arguments(args)
    positional = String[]
    options = Dict{String, Vector{String}}()
    i = 1
    while i <= length(args)
        arg = args[i]
        if !startswith(arg, "--")
            push!(positional, arg)
        else
            if contains(arg, "=")
                key, value = split(arg[3:end], "=", limit = 2)
            else
                i < length(args) || throw(UsageError("missing value for $arg"))
                key, value = arg[3:end], args[i += 1]
            end
            push!(get!(options, String(key), String[]), String(value))
        end
        i += 1
    end
    return positional, options
end

option(options, key, default) = haskey(options, key) ? last(options[key]) : default
function option(options, key)
    haskey(options, key) || throw(UsageError("--$key is required"))
    return last(options[key])
end

# same labels as manual_benchmarking.jl
function arch_label(arch)
    arch == "gpu" && return "gpu-nvidia"
    arch == "amdgpu" && return "gpu-amd"
    arch == "cpu" || throw(UsageError("unknown architecture $arch, use cpu, gpu or amdgpu"))
    arch_str = String(Sys.ARCH)
    return (startswith(arch_str, "aarch") || arch_str == "arm64") ? "cpu-arm" : "cpu-x86"
end

git(repo, args...) = readchomp(pipeline(`git -C $repo $args`, stderr = devnull))

function stored_timestamp(commit, label, repo)
    data = try
        JSON3.read(git(repo, "show", "$commit:$RESULTS_JSON"))
    catch
        return nothing
    end
    record = get(data, Symbol(label), nothing)
    isnothing(record) && return nothing
    return get(get(record, :meta, Dict()), :timestamp, nothing)
end

function benchmarked_commit(positional, options)
    label = haskey(options, "arch-label") ? option(options, "arch-label") : arch_label(option(options, "arch", "cpu"))
    ref = option(options, "ref", "origin/main")
    repo = option(options, "repo", ".")
    commits = split(git(repo, "log", "--first-parent", "--format=%H", ref, "--", RESULTS_JSON))
    current = isempty(commits) ? nothing : stored_timestamp(commits[1], label, repo)
    isnothing(current) && error("no stored benchmark results for $label on $ref")
    found = commits[1]
    for commit in commits[2:end]    # walk back while the stored record is unchanged
        stored_timestamp(commit, label, repo) == current || break
        found = commit
    end
    println(found)
    return 0
end

function stamp(positional, options)
    length(positional) == 1 || throw(UsageError("stamp takes exactly one FILE"))
    file = only(positional)
    data = copy(JSON3.read(read(file, String)))     # mutable Dict{Symbol, Any}
    tree = realpath(option(options, "tree"))
    dirs = get(data[:meta], :package_dirs, Dict{Symbol, Any}())
    wrong = filter(((name, dir),) -> !startswith(realpath(dir), tree * "/"), dirs)
    isempty(dirs) && error("no package_dirs recorded in $file")
    isempty(wrong) || error("packages not loaded from $tree: $wrong")
    data[:meta][:revision] = option(options, "revision")
    open(io -> JSON3.pretty(io, data), file, "w")
    return 0
end

"""(meta, Dict((truncation, nlayers, transform) => best SYPD)) over the comma-separated `files`."""
function load(files)
    meta = nothing
    best = Dict{Tuple{Int, Int, String}, Float64}()
    for file in split(files, ",")
        data = JSON3.read(read(file, String))
        meta = something(meta, data.meta)
        (; truncation, nlayers, spectral_transform, sypd) = data.overview
        for (t, l, transform, s) in zip(truncation, nlayers, spectral_transform, sypd)
            isnothing(s) && continue
            key = (Int(t), Int(l), String(transform))
            best[key] = max(get(best, key, 0.0), Float64(s))
        end
    end
    return meta, best
end

in_group(key, group) = group == "all" || key[3] == group

function geomean_ratio(candidate, reference, group)
    common = [key for key in keys(candidate) if haskey(reference, key) && in_group(key, group)]
    isempty(common) && return nothing
    return exp(sum(log(candidate[key] / reference[key]) for key in common) / length(common))
end

format_sypd(sypd) = isnothing(sypd) ? "—" : sypd < 10 ? @sprintf("%.1f", sypd) : @sprintf("%.0f", sypd)

function format_change(ratio, threshold)
    isnothing(ratio) && return "—"
    change = replace(@sprintf("%+.0f%%", (ratio - 1) * 100), "-" => "−")
    return ratio < threshold ? "**$change**" : change
end

markdown_row(cells) = "| " * join(cells, " | ") * " |"

function table(positional, options)
    haskey(options, "result") || throw(UsageError("table needs at least one --result NAME=FILE[,FILE...]"))
    threshold = parse(Float64, option(options, "threshold", "0.85"))
    candidate_name = option(options, "candidate", "main")

    names = String[]
    results = Dict{String, Any}()   # name => (meta, best, number of runs)
    for spec in options["result"]
        contains(spec, "=") || throw(UsageError("--result must be NAME=FILE[,FILE...], got $spec"))
        name, files = split(spec, "=", limit = 2)
        push!(names, name)
        results[name] = (load(files)..., length(split(files, ",")))
    end
    candidate_name in names || throw(UsageError("--candidate $candidate_name is not one of the --result names $names"))
    candidate = results[candidate_name][2]
    references = filter(!=(candidate_name), names)

    configurations = collect(union((keys(results[name][2]) for name in names)...))
    sort!(configurations, by = key -> (key[3] != "default", key[1], key[2]))
    header = ["T", "L", "Transform", names..., ("$candidate_name vs $name" for name in references)...]
    println(markdown_row(header))
    println("|", " --- |"^length(header))
    for key in configurations
        sypds = [format_sypd(get(results[name][2], key, nothing)) for name in names]
        changes = map(references) do name
            reference = results[name][2]
            ratio = haskey(candidate, key) && haskey(reference, key) ? candidate[key] / reference[key] : nothing
            format_change(ratio, threshold)
        end
        transform = key[3] == "matrix" ? "MT" : "LT+FFT"
        println(markdown_row([key[1], key[2], transform, sypds..., changes...]))
    end
    for group in GROUPS
        changes = [format_change(geomean_ratio(candidate, results[name][2], group), threshold) for name in references]
        println(markdown_row(["**geomean $(GROUP_LABELS[group])**", "", "", fill("", length(names))..., changes...]))
    end

    println("\nRevisions (SYPD = simulated years per wallclock day, higher is better, best of N runs):")
    for name in names
        meta, _, nruns = results[name]
        revision = first(string(get(meta, :revision, "?")), 10)
        println("- $name: $revision, v$(meta.speedyweather_version), $(get(meta, :arch_label, "?")), N=$nruns")
    end

    println("\nVerdict (regression = geomean ratio < $threshold):")
    for name in references
        ratios = [group => geomean_ratio(candidate, results[name][2], group) for group in GROUPS]
        filter!(!isnothing ∘ last, ratios)
        summary = join(("$(GROUP_LABELS[group]) $(@sprintf("%.3f", ratio))" for (group, ratio) in ratios), ", ")
        regressed = filter(<(threshold) ∘ last, ratios)
        if isempty(regressed)
            println("- vs $name: ok ($summary)")
        else
            group, ratio = regressed[argmin(last.(regressed))]
            cutoff = sqrt(ratio)    # geometric midpoint between reference (1) and candidate
            println("- vs $name: REGRESSION ($summary); bisect with --group $group --cutoff $(@sprintf("%.3f", cutoff))")
        end
    end
    return 0
end

function verdict(positional, options)
    _, reference = load(option(options, "reference"))
    _, candidate = load(option(options, "candidate"))
    group = option(options, "group", "all")
    group in GROUPS || throw(UsageError("--group must be one of $(join(GROUPS, ", "))"))
    cutoff = parse(Float64, option(options, "cutoff"))
    ratio = geomean_ratio(candidate, reference, group)
    isnothing(ratio) && error("no common configurations")
    good = ratio >= cutoff
    @printf("geomean %s ratio %.3f %s cutoff %.3f → %s\n", GROUP_LABELS[group], ratio, good ? "≥" : "<", cutoff, good ? "good" : "bad")
    return good ? 0 : 1
end

const COMMANDS = Dict(
    "benchmarked-commit" => benchmarked_commit,
    "stamp" => stamp,
    "table" => table,
    "verdict" => verdict,
)

function main(args)
    if isempty(args) || !haskey(COMMANDS, args[1])
        print(stderr, USAGE)
        return 2
    end
    positional, options = parse_arguments(args[2:end])
    return COMMANDS[args[1]](positional, options)
end

status = try
    main(ARGS)
catch err
    println(stderr, err isa UsageError ? "usage error: $(err.msg)" : "error: $(sprint(showerror, err))")
    2
end
exit(status)
