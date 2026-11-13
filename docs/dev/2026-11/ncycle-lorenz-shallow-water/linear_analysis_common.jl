using LinearAlgebra
wA(k,N) = k==0 ? 1.0 : N/(N-k)
# code family with per-substep θ_k: G_k = w(E x + L(x + θ_k (x'-x))) + (1-w)G ; x' = x + G
function cyc_theta(E, L, ws, θs)
    n = size(E,1); I_ = Matrix{ComplexF64}(I,n,n)
    T = Matrix{ComplexF64}(I,2n,2n)
    for (w,θ) in zip(ws,θs)
        S = inv(I_ - w*θ*L)
        dx = S*w*(E+L); dG = S*(1-w)
        T = [I_+dx dG; dx dG]*T
    end
    T
end
amp(E,L,ws,θs) = maximum(abs, eigvals(cyc_theta(E,L,ws,θs)))^(1/length(ws))
cycR(z, ws, θs) = cyc_theta(zeros(ComplexF64,1,1), fill(z,1,1), ws, θs)[1,1]  # x->x map (G erased by w0=1)
ys = range(0.05, 30, length=200)
loss(θs, ws) = sum((abs(cycR(im*y, ws, θs)) - 1)^2 for y in ys)
# crude optimization: Nelder-Mead by hand via random search + coordinate descent
function opt(ws; iters=4000)
    N = length(ws); best = fill(0.5, N); bl = loss(best, ws)
    step = 0.2
    for it in 1:iters
        c = best .+ step .* randn(N)
        l = loss(c, ws)
        if l < bl; best, bl = c, l; end
        it % 500 == 0 && (step *= 0.5)
    end
    best, bl
end
using Random; Random.seed!(1)
for N in (3,4)
    ws = [wA(k,N) for k in 0:N-1]
    for trial in 1:3
        θ, l = opt(ws)
        println("N=$N θ=", round.(θ,digits=4), " loss=", l)
        println("   |R| at y=0.5,1,3,10: ", [round(abs(cycR(im*y,ws,θ)),digits=5) for y in (0.5,1,3,10)], "  phase err y=0.5: ", angle(cycR(0.5im,ws,θ)) - N*0.5)
    end
end
