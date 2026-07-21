# MWE — EnzymeNoTypeError from a runtime-length ntuple of wrapped views passed to a
# KernelAbstractions kernel (Julia >= 1.11 only).
#
# Reduced from SpeedyWeather.jl's leapfrog step `update_prognostic!`
# (SpeedyWeather/src/time_stepping/steppers/leapfrog.jl + time_stepping/steps.jl `get_steps`).
# No SpeedyWeather dependency: only Enzyme + KernelAbstractions.
#
#   julia +1.12 --project=<env with Enzyme+KA> MWE_enzyme_runtime_ntuple.jl
#     -> EnzymeNoTypeError: "Enzyme cannot statically prove the type of a value ..."
#        (also mentions "very large sized registers that exceed the maximum size of
#         Enzyme's type analysis")
#   julia +1.10 (same script)  -> SUCCESS, correct gradient
#
# Versions: Enzyme v0.13.173, KernelAbstractions v0.9, Julia 1.12.6 (fails) / 1.10.11 (works).
#
# The three load-bearing ingredients (each verified necessary by ablation):
#   1. `ntuple(s -> view3(var, s), size(var, 3))` with a RUNTIME length: ntuple's runtime-n
#      branch ladder returns a small Union of tuple types of large view aggregates.
#      With a literal length (`ntuple(f, 2)`) the error disappears.
#   2. The views must be passed to a differentiated KernelAbstractions CPU kernel.
#      A plain Julia loop with the identical body differentiates fine.
#   3. The array wrapper must carry a second, inactive nested-struct field (here MySpectrum
#      with Ints + Vector{Int}s). A data-only wrapper, or a bare Vector{Int} second field,
#      differentiates fine.
# Everything else was ablated away: no runtime step INDEX is needed (all view offsets are
# literals below), no boundschecks, no mutable Const structs, no big Const args.
#
# Workarounds verified: Enzyme.API.looseTypeAnalysis!(true)  (correct gradients),
# or making the ntuple length compile-time. Enzyme.API.maxtypeoffset!/maxtypedepth! do NOT help
# (tested up to maxtypeoffset 2^24).

using Enzyme
using KernelAbstractions

println("Julia ", VERSION, ", Enzyme ", pkgversion(Enzyme))

# spectrum-like inactive metadata struct (ingredient 3: required)
struct MySpectrum{V <: AbstractVector}
    lmax::Int
    mmax::Int
    l_indices::V
    m_indices::V
end

# LowerTriangularArray-like wrapper: data + inactive metadata
struct MyLTA{T, N, A <: AbstractArray{T, N}, S <: MySpectrum} <: AbstractArray{T, N}
    data::A
    spectrum::S
end
Base.size(L::MyLTA) = size(L.data)
Base.@propagate_inbounds Base.getindex(L::MyLTA, i::Int) = getindex(L.data, i)
Base.@propagate_inbounds Base.setindex!(L::MyLTA, x, i::Int) = setindex!(L.data, x, i)
Base.IndexStyle(::Type{<:MyLTA}) = IndexLinear()

# view on the last (step) dimension, wrapped back (as SpeedyWeather's lta_view/get_step)
view3(L::MyLTA, s) = MyLTA(view(L.data, :, :, s), L.spectrum)

# ingredient 1: RUNTIME-length ntuple — size(var, 3) is not a compile-time constant.
# Replacing `size(var, 3)` with a literal `2` makes the error disappear.
steps_of(var) = ntuple(s -> view3(var, s), size(var, 3))

# ingredient 2: a KernelAbstractions kernel (leapfrog-shaped: reads + writes the views).
# The same body as a plain Julia loop differentiates fine.
@kernel inbounds = true function leap!(var_old, var_new, var_lf, tendency, Δt, w1, w2)
    lmk = @index(Global, Linear)
    old = var_old[lmk]
    new = old + Δt * tendency[lmk]
    update = old - 2var_lf[lmk] + new
    var_old[lmk] = var_lf[lmk] + w1 * update
    var_new[lmk] = new - w2 * update
end

function step!(var, tend)
    var_old, var_new = steps_of(var)     # runtime-length ntuple of wrapped views
    var_lf = view3(var, 2)               # literal offsets everywhere — no runtime index needed
    var_tend = view3(tend, 1)
    loop = leap!(KernelAbstractions.CPU(), (65, 1), size(var_tend))
    loop(var_old, var_new, var_lf, var_tend, 0.001f0, 0.01f0, 0.005f0)
    return nothing
end

spec = MySpectrum(11, 10, ones(Int, 65), ones(Int, 65))
var = MyLTA(rand(ComplexF32, 65, 1, 2), spec)    # (coeffs, layers, 2 time steps)
tend = MyLTA(rand(ComplexF32, 65, 1, 1), spec)   # (coeffs, layers, 1 step)
dvar = Enzyme.make_zero(var)
dvar.data[:, :, 2] .= 1                          # seed
dtend = Enzyme.make_zero(tend)

try
    autodiff(set_runtime_activity(Reverse), step!, Const,
        Duplicated(var, dvar), Duplicated(tend, dtend))
    println("RESULT: SUCCESS  (dtend extrema = ", extrema(real, dtend.data), ")")
catch e
    println("RESULT: FAILED -> ", first(sprint(showerror, e), 200))
end
