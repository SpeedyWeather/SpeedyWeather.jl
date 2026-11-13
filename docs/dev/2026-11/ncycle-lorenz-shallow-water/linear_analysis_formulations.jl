using LinearAlgebra
wA(k,N) = k==0 ? 1.0 : N/(N-k)
# general matrix version; E,L matrices (in units Δt=1)
function cyc(E, L, ws, scheme; α=0.5)
    n = size(E,1); I_ = Matrix{ComplexF64}(I,n,n); Z = zeros(ComplexF64,n,n)
    T = Matrix{ComplexF64}(I,2n,2n)
    S = inv(I_ - α*L); R = S*(I_ + (1-α)*L)
    for w in ws
        if scheme == :lawson   # G̃ = wEx + (1-w)G ; x = R(x+G̃) ; G = R G̃
            gx = w*E; gG = (1-w)*I_
            M = [R*(I_+gx) R*gG; R*gx R*gG]
        elseif scheme == :docs
            gx = w*E; gG = (1-w)*I_
            M = [I_+S*(gx+L) S*gG; gx gG]
        elseif scheme == :code
            Sw = inv(I_ - w*α*L)
            dx = Sw*w*(E+L); dG = Sw*(1-w)
            M = [I_+dx dG; dx dG]
        end
        T = M*T
    end
    maximum(abs, eigvals(T))^(1/length(ws))
end
function lf(E,L; α=0.5, ν=0.1, a=0.53)
    n=size(E,1); I_=Matrix{ComplexF64}(I,n,n)
    S = inv(I_-2α*L); nm = S*(I_+2(1-α)*L); n0 = 2S*E
    dm = ν/2*(I_+nm); d0 = ν/2*(-2I_+n0)
    M = [a*dm I_+a*d0; nm-(1-a)*dm n0-(1-a)*d0]
    maximum(abs,eigvals(M))
end
A4=[wA(k,4) for k in 0:3]; A3=[wA(k,3) for k in 0:2]
sc(oE,oL) = (fill(im*oE,1,1), fill(im*oL,1,1))
println("scalar Doppler: worst over ωL∈[0,40] for ωE")
for oE in (0.0,0.2,0.4,0.6,0.7,0.8,0.9)
    r(s,ws) = round(maximum(cyc(sc(oE,l)...,ws,s) for l in range(0,40,length=400)),digits=4)
    println(oE, " LF ", round(maximum(lf(sc(oE,l)...) for l in range(0,40,length=400)),digits=4),
     " lawsonA3 ", r(:lawson,A3), " lawsonA4 ", r(:lawson,A4), " docsA4 ", r(:docs,A4), " codeA4 ", r(:code,A4))
end
# f-plane: x=(ζ,δ,η) ; E = Coriolis; L = GW with ωL = c k Δt
println("f-plane fΔt, worst over ωL")
for f in (0.05,0.2,0.5)
    E = ComplexF64[0 -f 0; f 0 0; 0 0 0]
    r(s,ws) = round(maximum(cyc(E,ComplexF64[0 0 0;0 0 l;0 -l 0],ws,s) for l in range(0,40,length=400)),digits=4)
    println(f," LF ",round(maximum(lf(E,ComplexF64[0 0 0;0 0 l;0 -l 0]) for l in range(0,40,length=400)),digits=4),
      " lawsonA4 ",r(:lawson,A4)," lawsonA3 ",r(:lawson,A3)," docsA4 ",r(:docs,A4)," codeA4 ",r(:code,A4))
end
println("pure GW damping per step (lawson = 1 by construction):", cyc(sc(0,0.7)...,A4,:lawson))
