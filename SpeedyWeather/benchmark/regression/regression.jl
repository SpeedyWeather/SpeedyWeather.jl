const USAGE = """
Check SpeedyWeather.jl for performance regression with the debug mode of the benchmark suite
(`manual_benchmarking.jl --debug`: PrimitiveWet resolution sweep, 8 layers, truncation ≤ 128).

Usage: julia SpeedyWeather/benchmark/regression/regression.jl <command> [options]

    check [--candidate REV] [--arch cpu|gpu|amdgpu] [--ref origin/main] [--threshold 0.85]
          [--references benchmarked,release,main] [--branch-name NAME] [--output-dir DIR]
          [--bisect] [--no-confirm] [--no-fetch] [--fail-on-regression]
        Benchmark the candidate (default: main, e.g. HEAD for a PR branch) and as references main,
        the latest release and the latest benchmarked revision of this arch (or only those given
        with --references), and write report.md and summary.json to DIR (default WORKDIR/results),
        naming the candidate in the headline by its git ref or by --branch-name. A regression (geometric mean of the SYPD
        ratios candidate / reference < threshold, for all configurations or one of the two
        transforms) is confirmed by running both again, and with --bisect traced to its first bad
        commit. Exits with 1 on a regression with --fail-on-regression.
    benchmark REVISION|TREE OUTPUT.json [--arch cpu]
        Debug benchmark of one git revision, or of an existing checkout as is.
    bisect --good REV --bad REV --reference FILE[,FILE...] --cutoff X [--group all|default|matrix]
        First-parent git bisect for the first revision whose geometric-mean SYPD ratio against
        the reference falls below the cutoff.
    table --result NAME=FILE[,FILE...] ... [--candidate main] [--threshold 0.85]
        Comparison table of result files, e.g. after re-running a revision by hand.
    benchmarked-commit [--arch cpu] [--arch-label LABEL] [--ref origin/main]
        First-parent commit that introduced the stored results of an arch in benchmark_results.json.

All commands accept --workdir WORKDIR (default \$SPEEDY_BENCH_WORKDIR, else
\$TMPDIR/speedyweather-benchmark-regression). Revisions are benchmarked in detached worktrees
under WORKDIR, always with the benchmark harness of this checkout, so all revisions are measured
identically and the working tree is never touched. Requires Julia ≥ 1.11. Errors exit with 2.
"""

# own environment next to this script, instantiated on first use
if !isfile(joinpath(@__DIR__, "Manifest.toml"))
    import Pkg
    Pkg.activate(@__DIR__; io = devnull)
    Pkg.instantiate(; io = devnull)
end
pushfirst!(LOAD_PATH, @__DIR__)

using JSON3, Printf, TOML

git(dir, args...) = readchomp(`git -C $dir $args`)
is_ancestor(a, b) = success(`git -C $REPO merge-base --is-ancestor $a $b`)

const REPO = git(@__DIR__, "rev-parse", "--show-toplevel")
const BENCHMARK_DIR = dirname(@__DIR__)     # the harness used for every revision
const HARNESS_FILES = ("manual_benchmarking.jl", "benchmark_suite.jl", "define_benchmarks.jl")
const RESULTS_JSON = "SpeedyWeather/benchmark/assets/benchmark_results.json"
const AMDGPU_UUID = "21141c5a-9bdb-4563-92ae-f87d6854732e"     # not a dependency of the benchmark project
const GROUPS = ("all", "default", "matrix")
const GROUP_LABELS = Dict("all" => "all", "default" => "LT+FFT", "matrix" => "MT")
const REFERENCES = ("benchmarked", "release", "main")

default_workdir() = get(ENV, "SPEEDY_BENCH_WORKDIR", joinpath(tempdir(), "speedyweather-benchmark-regression"))
short(sha) = first(sha, 10)

# same labels as manual_benchmarking.jl
function arch_label(arch)
    arch == "gpu" && return "gpu-nvidia"
    arch == "amdgpu" && return "gpu-amd"
    arch == "cpu" || throw(ArgumentError("unknown architecture $arch, use cpu, gpu or amdgpu"))
    arch_str = String(Sys.ARCH)
    return (startswith(arch_str, "aarch") || arch_str == "arm64") ? "cpu-arm" : "cpu-x86"
end

# REVISIONS

latest_release() = first(split(git(REPO, "tag", "--list", "v[0-9]*", "--sort=-v:refname")))

function stored_timestamp(commit, label)
    data = try
        JSON3.read(read(pipeline(`git -C $REPO show $commit:$RESULTS_JSON`, stderr = devnull), String))
    catch
        return nothing
    end
    record = get(data, Symbol(label), nothing)
    return isnothing(record) ? nothing : get(get(record, :meta, Dict()), :timestamp, nothing)
end

"""First-parent commit on `ref` that introduced the currently stored benchmark results of `label`."""
function benchmarked_commit(label; ref = "origin/main")
    commits = split(git(REPO, "log", "--first-parent", "--format=%H", ref, "--", RESULTS_JSON))
    current = isempty(commits) ? nothing : stored_timestamp(commits[1], label)
    isnothing(current) && error("no stored benchmark results for $label on $ref")
    found = commits[1]
    for commit in commits[2:end]    # walk back while the stored record is unchanged
        stored_timestamp(commit, label) == current || break
        found = commit
    end
    return String(found)
end

# BENCHMARKING

function worktree(revision, workdir)
    sha = git(REPO, "rev-parse", "--verify", "$revision^{commit}")
    tree = joinpath(workdir, "trees", short(sha))
    if !isdir(tree)
        git(REPO, "worktree", "prune")
        git(REPO, "worktree", "add", "--detach", "--quiet", tree, sha)
    end
    return tree
end

function remove_worktrees(workdir)
    trees = joinpath(workdir, "trees")
    for tree in [(isdir(trees) ? readdir(trees, join = true) : String[]); joinpath(workdir, "bisect")]
        isdir(tree) && run(pipeline(ignorestatus(`git -C $REPO worktree remove --force $tree`), stderr = devnull))
    end
    git(REPO, "worktree", "prune")
    return nothing
end

# environment with only what the debug mode loads, SpeedyWeather from `tree`
function write_project(env, tree, arch)
    deps = TOML.parsefile(joinpath(BENCHMARK_DIR, "Project.toml"))["deps"]
    project_deps = Dict(name => deps[name] for name in ("BenchmarkTools", "Dates", "JSON3", "Printf", "SpeedyWeather"))
    arch == "gpu" && (project_deps["CUDA"] = deps["CUDA"])
    arch == "amdgpu" && (project_deps["AMDGPU"] = AMDGPU_UUID)
    sources = Dict("SpeedyWeather" => Dict("path" => joinpath(tree, "SpeedyWeather")))
    open(io -> TOML.print(io, Dict("deps" => project_deps, "sources" => sources)), joinpath(env, "Project.toml"), "w")
    return nothing
end

# check that the tree's packages were benchmarked (not e.g. a registered release) and record the revision
function stamp!(output, tree, sha)
    data = copy(JSON3.read(read(output, String)))   # mutable Dict{Symbol, Any}
    dirs = get(data[:meta], :package_dirs, Dict{Symbol, Any}())
    isempty(dirs) && error("no package_dirs recorded in $output")
    wrong = filter(((name, dir),) -> !startswith(realpath(dir), realpath(tree) * "/"), dirs)
    isempty(wrong) || error("packages not loaded from $tree: $wrong")
    data[:meta][:revision] = sha
    open(io -> JSON3.pretty(io, data), output, "w")
    return nothing
end

"""Debug benchmark of `target`, a git revision (checked out in a detached worktree under `workdir`)
or an existing checkout, with the harness of this checkout. Writes `output` and `output.log`,
returns whether it succeeded."""
function benchmark(target, output; arch = "cpu", workdir = default_workdir())
    VERSION >= v"1.11" || error("Julia ≥ 1.11 is required ([sources] is ignored before), this is $VERSION")
    tree = isdir(target) ? realpath(target) : worktree(target, workdir)
    sha = git(tree, "rev-parse", "HEAD")
    env = joinpath(workdir, "envs", "$(short(sha))-$arch")
    mkpath(env)
    for file in HARNESS_FILES
        cp(joinpath(BENCHMARK_DIR, file), joinpath(env, file), force = true)
    end
    write_project(env, tree, arch)

    output = abspath(output)
    mkpath(dirname(output))
    @info "Benchmarking $(short(sha)) ($(git(tree, "log", "-1", "--format=%s"))) on $arch → $output"
    julia = `$(Base.julia_cmd()) --startup-file=no --project=$env`
    commands = (
        `$julia -e "using Pkg; Pkg.resolve(); Pkg.instantiate()"`,
        `$julia $(joinpath(env, "manual_benchmarking.jl")) $arch --debug --output=$output`,
    )
    succeeded = open(output * ".log", "w") do log
        all(command -> success(pipeline(command, stdout = log, stderr = log)), commands)
    end
    if !succeeded
        @warn "Benchmark of $(short(sha)) failed, last lines of $output.log:\n" * join(last(readlines(output * ".log"), 30), "\n")
        return false
    end
    try
        stamp!(output, tree, sha)
    catch err
        @warn "Benchmark of $(short(sha)) cannot be trusted" exception = err
        return false
    end
    return true
end

# COMPARISON

"""Results of one revision: SYPD per (truncation, nlayers, transform), the best (max) over all
`files`, as noise from other processes only ever slows a run down."""
struct Measurement
    name::String
    files::Vector{String}
    meta::Any
    sypd::Dict{Tuple{Int, Int, String}, Float64}
end

function Measurement(name, files)
    meta = nothing
    sypd = Dict{Tuple{Int, Int, String}, Float64}()
    for file in files
        data = JSON3.read(read(file, String))
        meta = something(meta, data.meta)
        (; truncation, nlayers, spectral_transform) = data.overview
        for (t, l, transform, s) in zip(truncation, nlayers, spectral_transform, data.overview.sypd)
            isnothing(s) && continue
            key = (Int(t), Int(l), String(transform))
            sypd[key] = max(get(sypd, key, 0.0), Float64(s))
        end
    end
    return Measurement(name, collect(files), meta, sypd)
end

revision(m::Measurement) = string(get(m.meta, :revision, "?"))

in_group(key, group) = group == "all" || key[3] == group

function geomean_ratio(candidate::Measurement, reference::Measurement, group)
    common = [key for key in keys(candidate.sypd) if haskey(reference.sypd, key) && in_group(key, group)]
    isempty(common) && return nothing
    return exp(sum(log(candidate.sypd[key] / reference.sypd[key]) for key in common) / length(common))
end

"""Geometric-mean ratios of `candidate` vs `reference` per group; on a regression (any ratio
below `threshold`) the worst group and a bisect cutoff halfway (geometrically) between both."""
function compare(candidate::Measurement, reference::Measurement, threshold)
    ratios = Pair{String, Float64}[]
    for group in GROUPS
        ratio = geomean_ratio(candidate, reference, group)
        isnothing(ratio) || push!(ratios, group => ratio)
    end
    regressed = filter(<(threshold) ∘ last, ratios)
    group, ratio = isempty(regressed) ? ("all", 1.0) : regressed[argmin(last.(regressed))]
    return (; reference, ratios, regressed = !isempty(regressed), group, ratio, cutoff = sqrt(ratio))
end

format_sypd(sypd) = isnothing(sypd) ? "—" : sypd < 10 ? @sprintf("%.1f", sypd) : @sprintf("%.0f", sypd)
format_percent(ratio) = replace(@sprintf("%+.0f%%", (ratio - 1) * 100), "-" => "−")

function format_change(ratio, threshold)
    isnothing(ratio) && return "—"
    return ratio < threshold ? "**$(format_percent(ratio))**" : format_percent(ratio)
end

markdown_row(cells) = "| " * join(cells, " | ") * " |\n"

"""Markdown table of SYPD per configuration and the change of `candidate` vs every reference."""
function comparison_table(candidate::Measurement, references, threshold)
    measurements = [references; candidate]
    configurations = collect(union((keys(m.sypd) for m in measurements)...))
    sort!(configurations, by = key -> (key[3] != "default", key[1], key[2]))
    layers = length(unique(key[2] for key in configurations)) > 1     # column only if it varies
    io = IOBuffer()
    header = ["T", (layers ? ["L"] : [])..., "Transform", (m.name for m in measurements)..., ("$(candidate.name) vs $(m.name)" for m in references)...]
    print(io, markdown_row(header), "|", " --- |"^length(header), "\n")
    for key in configurations
        sypds = [format_sypd(get(m.sypd, key, nothing)) for m in measurements]
        changes = map(references) do reference
            both = haskey(candidate.sypd, key) && haskey(reference.sypd, key)
            format_change(both ? candidate.sypd[key] / reference.sypd[key] : nothing, threshold)
        end
        transform = key[3] == "matrix" ? "MT" : "LT+FFT"
        print(io, markdown_row([key[1], (layers ? [key[2]] : [])..., transform, sypds..., changes...]))
    end
    for group in GROUPS
        changes = [format_change(geomean_ratio(candidate, reference, group), threshold) for reference in references]
        empty_cells = fill("", length(measurements) + (layers ? 2 : 1))
        print(io, markdown_row(["**geomean $(GROUP_LABELS[group])**", empty_cells..., changes...]))
    end

    println(io, "\nSYPD = simulated years per wallclock day, higher is better; best of N runs.")
    for m in measurements
        println(io, "- $(m.name): $(short(revision(m))), v$(m.meta.speedyweather_version), $(get(m.meta, :arch_label, "?")), N=$(length(m.files))")
    end
    return String(take!(io))
end

function verdict_lines(comparisons, threshold)
    io = IOBuffer()
    for c in comparisons
        summary = join(("$(GROUP_LABELS[group]) $(@sprintf("%.3f", ratio))" for (group, ratio) in c.ratios), ", ")
        if c.regressed
            println(io, "- vs $(c.reference.name): **regression** ($summary), bisect on $(c.group) with cutoff $(@sprintf("%.3f", c.cutoff))")
        else
            println(io, "- vs $(c.reference.name): ok ($summary)")
        end
    end
    return String(take!(io))
end

# BISECTION

const BISECT_LINE = r"^# (good|bad|skip|first bad commit|possible first bad commit): \[([0-9a-f]{40})\] (.*)$"m

"""First-parent `git bisect` between `good` and `bad` in a worktree under `workdir`. Every step is
benchmarked into `output_dir/bisect-<sha>.json` and is bad if its geometric-mean SYPD ratio (of
`group`) against the `reference` files is below `cutoff`. Returns the culprit (or `nothing`), the
remaining candidates if only skipped commits are left, and the tested steps."""
function bisect(;
        good, bad, reference, cutoff, group = "all", arch = "cpu",
        workdir = default_workdir(), output_dir = joinpath(workdir, "results"),
    )
    tree = joinpath(workdir, "bisect")
    isdir(tree) && run(ignorestatus(`git -C $REPO worktree remove --force $tree`))
    git(REPO, "worktree", "prune")
    git(REPO, "worktree", "add", "--detach", "--quiet", tree, bad)
    step = `$(Base.julia_cmd()) --startup-file=no $(@__FILE__) bisect-step`
    step = `$step --reference $(join(reference, ",")) --cutoff $cutoff --group $group`
    step = `$step --arch $arch --workdir $workdir --output-dir $output_dir`
    log = try
        run(`git -C $tree bisect start --first-parent $bad $good`)
        run(ignorestatus(`git -C $tree bisect run $step`))
        git(tree, "bisect", "log")
    finally
        run(pipeline(ignorestatus(`git -C $tree bisect reset`), stdout = devnull, stderr = devnull))
        run(ignorestatus(`git -C $REPO worktree remove --force $tree`))
    end
    entries = [(; kind = m[1], sha = m[2], subject = m[3]) for m in eachmatch(BISECT_LINE, log)]
    culprit = findfirst(e -> e.kind == "first bad commit", entries)
    return (;
        good, bad, group, cutoff, log,
        culprit = isnothing(culprit) ? nothing : entries[culprit],
        possible = filter(e -> e.kind == "possible first bad commit", entries),
        steps = filter(e -> e.kind in ("good", "bad", "skip"), entries),
    )
end

# run by `git bisect run` in the bisect worktree: exit 0 = good, 1 = bad, 125 = skip
function bisect_step(; reference, cutoff, group, arch, workdir, output_dir)
    tree = git(pwd(), "rev-parse", "--show-toplevel")
    sha = git(tree, "rev-parse", "HEAD")
    output = joinpath(output_dir, "bisect-$(short(sha)).json")
    isfile(output) || benchmark(tree, output; arch, workdir) || return 125
    ratio = geomean_ratio(Measurement(short(sha), [output]), Measurement("reference", reference), group)
    isnothing(ratio) && return 125
    @printf("%s: geomean %s ratio %.3f vs cutoff %.3f → %s\n", short(sha), GROUP_LABELS[group], ratio, cutoff, ratio >= cutoff ? "good" : "bad")
    return ratio >= cutoff ? 0 : 1
end

function bisect_report(result, measurements, output_dir)
    (; reference) = result
    io = IOBuffer()
    println(io, "\n### Bisection\n")
    cutoff = @sprintf("%.3f", result.cutoff)
    println(io, "First parent of $(result.candidate_name) from $(short(result.good)) (good) to $(short(result.bad)) (bad); a revision is bad if its ")
    println(io, "geomean $(GROUP_LABELS[result.group]) SYPD ratio vs $(reference.name) is below $cutoff.\n")
    print(io, markdown_row(["commit", "ratio", "verdict"]), "| --- | --- | --- |\n")
    files = Dict(revision(m) => m.files for m in measurements)
    for step in result.steps
        file = get(files, step.sha, [joinpath(output_dir, "bisect-$(short(step.sha)).json")])
        ratio = all(isfile, file) ? geomean_ratio(Measurement(step.sha, file), reference, result.group) : nothing
        print(io, markdown_row(["$(short(step.sha)) $(step.subject)", isnothing(ratio) ? "—" : @sprintf("%.3f", ratio), step.kind]))
    end
    if !isnothing(result.culprit)
        println(io, "\n**First bad commit:** $(short(result.culprit.sha)) $(result.culprit.subject)")
    elseif !isempty(result.possible)
        println(io, "\n**Only skipped commits left**, the first bad commit is one of:")
        foreach(e -> println(io, "- $(short(e.sha)) $(e.subject)"), result.possible)
    else
        println(io, "\n**Bisection inconclusive**, see the bisect log in summary.json.")
    end
    return String(take!(io))
end

# CHECK

"""Benchmark `candidate` (default: main) and as references main, the latest release and the latest
benchmarked revision, compare, confirm a regression with a second run and optionally bisect it.
Writes report.md and summary.json to `output_dir` and returns whether the candidate regressed."""
function check(;
        arch = "cpu", ref = "origin/main", candidate = ref, threshold = 0.85, workdir = default_workdir(),
        output_dir = joinpath(workdir, "results"), confirm = true, bisect_regression = false, fetch = true,
        reference_names = REFERENCES, branch = nothing,
    )
    unknown = setdiff(reference_names, REFERENCES)
    isempty(unknown) || throw(ArgumentError("unknown references $(join(unknown, ", ")), use $(join(REFERENCES, ", "))"))
    fetch && git(REPO, "fetch", "origin", "--tags", "--quiet")
    label = arch_label(arch)
    main_sha = git(REPO, "rev-parse", "$ref^{commit}")
    candidate_sha = git(REPO, "rev-parse", "$candidate^{commit}")
    candidate_name = candidate_sha == main_sha ? "main" : "branch"
    # for the headline, e.g. the PR branch
    branch = isnothing(branch) ? git(REPO, "rev-parse", "--abbrev-ref", candidate) : branch
    notes = String[]

    named = Pair{String, String}[]
    if "benchmarked" in reference_names
        try
            push!(named, "benchmarked" => benchmarked_commit(label; ref))
        catch err
            push!(notes, sprint(showerror, err))     # e.g. no stored results for this arch yet
        end
    end
    if "release" in reference_names
        release = latest_release()
        push!(named, "release $release" => git(REPO, "rev-parse", "$release^{commit}"))
    end
    "main" in reference_names && push!(named, "main" => main_sha)

    # references, merged if they are the same commit, and dropped if they are the candidate
    references = Pair{String, String}[]
    for (name, sha) in named
        if sha == candidate_sha
            name == "main" || push!(notes, "$name is the same commit as $candidate_name ($(short(sha))).")
            continue
        end
        i = findfirst(r -> r.second == sha, references)
        isnothing(i) ? push!(references, name => sha) : (references[i] = "$(references[i].first) = $name" => sha)
    end

    file(name, attempt) = joinpath(output_dir, replace(name, r"[^\w.-]+" => "_") * "-$attempt.json")
    try
        measured = Measurement[]
        for (name, sha) in references
            if benchmark(sha, file(name, 1); arch, workdir)
                push!(measured, Measurement(name, [file(name, 1)]))
            else
                push!(notes, "$name ($(short(sha))) could not be benchmarked with the current harness, see $(file(name, 1)).log")
            end
        end
        benchmark(candidate_sha, file(candidate_name, 1); arch, workdir) ||
            error("benchmarking $candidate_name failed, see $(file(candidate_name, 1)).log")
        tested = Measurement(candidate_name, [file(candidate_name, 1)])
        comparisons = [compare(tested, reference, threshold) for reference in measured]

        if confirm && any(c -> c.regressed, comparisons)
            @info "Possible regression, benchmarking $candidate_name and the regressed references again to confirm"
            for m in [tested; [c.reference for c in comparisons if c.regressed]]
                benchmark(revision(m), file(m.name, 2); arch, workdir) && push!(m.files, file(m.name, 2))
            end
            measured = [Measurement(m.name, m.files) for m in measured]
            tested = Measurement(tested.name, tested.files)
            comparisons = [compare(tested, reference, threshold) for reference in measured]
        end
        regressed = filter(c -> c.regressed, comparisons)

        bisection = nothing
        if bisect_regression && !isempty(regressed)
            # bisect from the most recent regressed reference (or where the candidate branched off it)
            newest = reduce((a, b) -> is_ancestor(revision(a.reference), revision(b.reference)) ? b : a, regressed)
            good = git(REPO, "merge-base", revision(newest.reference), candidate_sha)
            # where the candidate branched off the reference may already be slow, e.g. if the
            # reference fixed a regression later; bisecting from there would blame the wrong commit
            base_file = joinpath(output_dir, "bisect-$(short(good)).json")
            base_is_slow = good != revision(newest.reference) && benchmark(good, base_file; arch, workdir) &&
                geomean_ratio(Measurement("base", [base_file]), newest.reference, newest.group) < newest.cutoff
            if good == candidate_sha
                push!(notes, "Nothing to bisect, $candidate_name is an ancestor of $(newest.reference.name).")
            elseif base_is_slow
                push!(notes, "Not bisected: $(short(good)), where $candidate_name branched off $(newest.reference.name), is already slow. Merge or rebase onto $(newest.reference.name) first.")
            else
                bisection = bisect(;
                    good, bad = candidate_sha, reference = newest.reference.files, newest.cutoff, newest.group,
                    arch, workdir, output_dir,
                )
                bisection = (; bisection..., reference = newest.reference, candidate_name)
            end
        end

        headline = if isempty(measured)
            "nothing to compare $candidate_name with"
        elseif isempty(regressed)
            "no significant regression"
        else
            worst = regressed[argmin([c.ratio for c in regressed])]
            slower = @sprintf("%.0f%%", (1 - worst.ratio) * 100)
            "**regression**, $candidate_name is $slower slower than $(worst.reference.name) (geomean $(GROUP_LABELS[worst.group]))"
        end
        of = candidate_name == "main" ? "" : " of `$branch`"
        report = IOBuffer()
        println(report, "## Benchmark regression check$of: $headline\n")
        date = Libc.strftime("%Y-%m-%d", time())
        println(report, "`manual_benchmarking.jl --debug` (PrimitiveWet, L8, T ≤ 128) on $label, $(Sys.cpu_info()[1].model), Julia $VERSION, $date.\n")
        print(report, comparison_table(tested, measured, threshold))
        println(report, "\nRegression = geometric-mean SYPD ratio $candidate_name / reference < $threshold:")
        print(report, verdict_lines(comparisons, threshold))
        foreach(note -> println(report, "\nNote: ", note), notes)
        isnothing(bisection) || print(report, bisect_report(bisection, [measured; tested], output_dir))
        report = String(take!(report))

        mkpath(output_dir)
        write(joinpath(output_dir, "report.md"), report)
        summary = Dict(
            "regression" => !isempty(regressed),
            "candidate" => Dict("name" => candidate_name, "ref" => branch, "revision" => candidate_sha),
            "arch_label" => label,
            "threshold" => threshold,
            "measurements" => [Dict("name" => m.name, "revision" => revision(m), "files" => m.files) for m in [measured; tested]],
            "comparisons" => [
                Dict("reference" => c.reference.name, "ratios" => Dict(c.ratios), "regressed" => c.regressed, "group" => c.group, "cutoff" => c.cutoff)
                    for c in comparisons
            ],
            "bisect" => isnothing(bisection) ? nothing : Dict(
                    "good" => bisection.good, "bad" => bisection.bad, "group" => bisection.group, "cutoff" => bisection.cutoff,
                    "culprit" => isnothing(bisection.culprit) ? nothing : bisection.culprit.sha,
                    "possible" => [e.sha for e in bisection.possible], "log" => bisection.log,
                ),
            "notes" => notes,
        )
        open(io -> JSON3.pretty(io, summary), joinpath(output_dir, "summary.json"), "w")
        println(report)
        @info "Wrote $(joinpath(output_dir, "report.md")) and summary.json"
        return !isempty(regressed)
    finally
        remove_worktrees(workdir)
    end
end

# COMMAND LINE

const SWITCHES = ("bisect", "no-confirm", "no-fetch", "fail-on-regression")

# `--key value` or `--key=value` (repeatable), `--switch` and positional arguments
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
            elseif arg[3:end] in SWITCHES
                key, value = arg[3:end], "true"
            else
                i < length(args) || throw(ArgumentError("missing value for $arg"))
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
    haskey(options, key) || throw(ArgumentError("--$key is required"))
    return last(options[key])
end
comma_separated(value) = String.(split(value, ","))

function main(command, args)
    positional, options = parse_arguments(args)
    workdir = abspath(option(options, "workdir", default_workdir()))
    output_dir = abspath(option(options, "output-dir", joinpath(workdir, "results")))
    arch = option(options, "arch", "cpu")
    threshold = parse(Float64, option(options, "threshold", "0.85"))

    if command == "check"
        regressed = check(;
            arch, workdir, output_dir, threshold, ref = option(options, "ref", "origin/main"),
            candidate = option(options, "candidate", option(options, "ref", "origin/main")),
            confirm = !haskey(options, "no-confirm"), bisect_regression = haskey(options, "bisect"),
            fetch = !haskey(options, "no-fetch"), branch = option(options, "branch-name", nothing),
            reference_names = comma_separated(option(options, "references", join(REFERENCES, ","))),
        )
        return regressed && haskey(options, "fail-on-regression") ? 1 : 0

    elseif command == "benchmark"
        length(positional) == 2 || throw(ArgumentError("benchmark takes REVISION|TREE and OUTPUT.json"))
        try
            return benchmark(positional...; arch, workdir) ? 0 : 2
        finally
            remove_worktrees(workdir)
        end

    elseif command == "bisect"
        reference = Measurement("reference", comma_separated(option(options, "reference")))
        result = bisect(;
            good = git(REPO, "rev-parse", option(options, "good")), bad = git(REPO, "rev-parse", option(options, "bad")),
            reference = reference.files, cutoff = parse(Float64, option(options, "cutoff")),
            group = option(options, "group", "all"), arch, workdir, output_dir,
        )
        print(bisect_report((; result..., reference, candidate_name = option(options, "bad")), Measurement[], output_dir))
        return 0

    elseif command == "bisect-step"
        return bisect_step(;
            reference = comma_separated(option(options, "reference")), cutoff = parse(Float64, option(options, "cutoff")),
            group = option(options, "group", "all"), arch, workdir, output_dir,
        )

    elseif command == "table"
        haskey(options, "result") || throw(ArgumentError("table needs --result NAME=FILE[,FILE...]"))
        measurements = map(options["result"]) do spec
            contains(spec, "=") || throw(ArgumentError("--result must be NAME=FILE[,FILE...], got $spec"))
            name, files = split(spec, "=", limit = 2)
            Measurement(name, comma_separated(files))
        end
        candidate_name = option(options, "candidate", "main")
        i = findfirst(m -> m.name == candidate_name, measurements)
        isnothing(i) && throw(ArgumentError("--candidate $candidate_name is not one of the --result names"))
        candidate = measurements[i]
        references = deleteat!(copy(measurements), i)
        print(comparison_table(candidate, references, threshold))
        println("\nRegression = geometric-mean SYPD ratio $(candidate.name) / reference < $threshold:")
        print(verdict_lines([compare(candidate, reference, threshold) for reference in references], threshold))
        return 0

    elseif command == "benchmarked-commit"
        label = haskey(options, "arch-label") ? option(options, "arch-label") : arch_label(arch)
        println(benchmarked_commit(label; ref = option(options, "ref", "origin/main")))
        return 0
    end
    print(stderr, USAGE)
    return 2
end

if abspath(PROGRAM_FILE) == @__FILE__
    command = isempty(ARGS) ? "" : ARGS[1]
    status = try
        main(command, ARGS[2:end])
    catch err
        println(stderr, "error: ", sprint(showerror, err))
        command == "bisect-step" ? 125 : 2      # an error must not mark a bisect step as bad
    end
    exit(status)
end
