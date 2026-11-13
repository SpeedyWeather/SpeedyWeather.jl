include("gw.jl")
N = parse(Int, ARGS[1]); θ = Float32.(parse.(Float64, split(ARGS[2], ","))); Δt = Minute(parse(Int, ARGS[3]))
r = try rungw(; stepper=:ncycle, steps=N, θ, Δt) catch e; (ratio_end=NaN, rmax=NaN) end
@printf("N=%d θ=%s Δt=%3d: var=%.4f max=%.2f\n", N, ARGS[2], Δt.value, r.ratio_end, r.rmax)
