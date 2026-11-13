include("linear_analysis_common.jl")
include("linear_analysis_formulations.jl")
using Base.Threads
println("LF β<0: ", [round(maximum(lf(fill(-β*im*l,1,1), fill(im*l,1,1)) for l in range(0,40,length=400)),digits=3) for β in (-0.05,-0.1,-0.25)])
Ls = [range(0,3,length=30); range(3,40,length=30)]
function score(θs, ws)
    N = length(ws)
    gwmax = maximum(abs(cycR(im*y, ws, θs)) for y in range(0.05,40,length=80))^(1/N)
    bet = maximum(amp(fill(-β*im*l,1,1), fill(im*l,1,1), ws, θs) for β in (0.1,0.3,0.5,0.7,0.9), l in Ls)
    dop = maximum(amp(fill(im*e,1,1), fill(im*l,1,1), ws, θs) for e in (0.1,0.2,0.3,0.4), l in Ls)
    fpl = maximum(amp(ComplexF64[0 -f 0; f 0 0; 0 0 0], ComplexF64[0 0 0;0 0 l;0 -l 0], ws, θs) for f in (0.1,0.35), l in Ls)
    damp = minimum(abs(cycR(im*y, ws, θs)) for y in range(0.05,1.5,length=15))^(1/N)
    m = max(gwmax, bet, dop, fpl)
    return -damp + (m > 1+1e-6 ? 1000*(m-1) : 0.0), (; damp, gwmax, bet, dop, fpl)
end
function opt(ws, init; iters=2000)
    best = copy(init); bs, = score(best, ws); step = 0.4
    for it in 1:iters
        c = best .+ step .* randn(length(ws)); s, = score(c, ws)
        s < bs && ((best, bs) = (c, s))
        it % 400 == 0 && (step *= 0.5)
    end
    best, score(best, ws)
end
for N in (3, 4)
    ws = [wA(k,N) for k in 0:N-1]
    println("N=$N old α=0.5: ", map(x->round(x,digits=4), score(fill(0.5,N), ws)[2]))
    res = Vector{Any}(undef, 28)
    @threads for t in 1:28
        res[t] = opt(ws, 0.6 .+ 0.8 .* randn(N))
    end
    sort!(res, by = r -> r[2][1])
    for r in res[1:4]
        println("N=$N θ=", round.(r[1],digits=4), " ", map(x->round(x,digits=5), r[2][2]))
    end
end
