include("gw.jl")
# args: stepper N μ α Δt_min [days] [variant]
a = ARGS
stepper = Symbol(a[1]); N = parse(Int, a[2]); μ = a[3] == "-" ? nothing : parse(Float32, a[3]); α = parse(Float64, a[4]); Δt = Minute(parse(Int, a[5]))
days = length(a) > 5 ? parse(Int, a[6]) : 2
variant = length(a) > 6 ? eval(Symbol(a[7]))() : NCycleLorenzA()
r = try rungw(; stepper, steps=N, μ, centering=α, Δt, days, variant) catch e; (ratio_end=NaN, rmax=NaN, r=[NaN]) end
@printf("%s N=%d μ=%s α=%.2f Δt=%3d %s: var=%.4f max=%.2f\n", stepper, N, a[3], α, Δt.value, length(a)>6 ? a[7] : "A", r.ratio_end, r.rmax)
