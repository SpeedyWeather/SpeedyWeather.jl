# Cost of SubcycledVerticalAdvection relative to WENOVerticalAdvection.
# usage: julia --project=<env with SpeedyWeather [+ CUDA]> subcycled_vertical_advection.jl [cpu|gpu] [spinup_hours]
# Spins up PrimitiveWetModel with SubcycledVerticalAdvection, then times vertical_advection!(vars, model)
# for each scheme on the same state, and reports the distribution of substeps over columns.
using SpeedyWeather, Printf, Statistics
arch_name = get(ARGS, 1, "cpu")
spinup_hours = parse(Int, get(ARGS, 2, "24"))
if arch_name == "gpu"
    using CUDA
    const ARCH = SpeedyWeather.GPU()
    sync() = CUDA.synchronize()
else
    const ARCH = SpeedyWeather.CPU()
    sync() = nothing
end

function timeit(f, n)
    f(); sync()                       # compile / warm up
    t = time_ns()
    for _ in 1:n
        f()
    end
    sync()
    return (time_ns() - t) / n / 1.0e6  # ms per call
end

function substep_distribution(vars, model, adv)
    w = Array(vars.dynamics.w.data)
    Δσ = Array(model.geometry.σ_levels_thick)
    interval = SpeedyWeather.advection_interval(adv, vars, model)
    nl = size(w, 2)
    C = [maximum(abs(w[ij, k]) * interval / min(Δσ[k], Δσ[k + 1]) for k in 1:(nl - 1)) for ij in 1:size(w, 1)]
    n = clamp.(ceil.(Int, C ./ adv.courant_target), 1, adv.substeps)
    return n, maximum(C)
end

println("arch = $arch_name, spin-up $spinup_hours h, Julia $(VERSION), threads $(Threads.nthreads())")
@printf("%5s %4s %8s | %9s %9s %9s %9s | %6s %6s %6s %6s | %7s %8s\n", "T", "L", "Δt [s]",
    "WENO", "sub ad.", "sub n=1", "sub n=4", "n=1 %", "mean n", "max n", "C max", "step", "sub/step")
for trunc in (31, 127), nlayers in (8, 100)
    spectral_grid = SpectralGrid(; truncation = trunc, nlayers, architecture = ARCH)
    sub = SubcycledVerticalAdvection(spectral_grid)
    model = PrimitiveWetModel(spectral_grid; vertical_advection = sub)
    model.feedback.verbose = false
    simulation = initialize!(model)
    ok = true
    try
        run!(simulation, period = Hour(spinup_hours))
    catch e
        ok = false
        @warn "spin-up failed at T$trunc L$nlayers" e
    end
    vars = simulation.variables
    ok = ok && all(isfinite, Array(vars.grid.temperature.data))

    schemes = (
        WENOVerticalAdvection(spectral_grid), sub,
        SubcycledVerticalAdvection(spectral_grid; courant_target = nothing, substeps = 1),
        SubcycledVerticalAdvection(spectral_grid; courant_target = nothing, substeps = 4),
    )
    nrep = trunc > 100 ? 10 : 50
    times = map(schemes) do adv
        m = PrimitiveWetModel(spectral_grid; vertical_advection = adv)   # same grid, variables are reused
        timeit(() -> SpeedyWeather.vertical_advection!(vars, m), nrep)
    end
    n, Cmax = substep_distribution(vars, model, sub)

    # full time step with the subcycled scheme
    Δt = SpeedyWeather.default_time_step(model.time_stepping) / 2
    tstep = timeit(() -> run!(simulation, steps = 10), 2) / 10
    @printf("%5d %4d %8.0f | %9.3f %9.3f %9.3f %9.3f | %6.1f %6.2f %6d %6.2f | %7.1f %7.1f%%  %s\n",
        trunc, nlayers, Δt, times..., 100 * mean(n .== 1), mean(n), maximum(n), Cmax, tstep,
        100 * times[2] / tstep, ok ? "" : "(state not finite)")
end
println("times in ms per vertical_advection! call (u, v, T, q); step = ms per full time step with the subcycled scheme")
