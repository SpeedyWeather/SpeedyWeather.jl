# MWE: EnzymeNoTypeError on Julia >= 1.11 (works on 1.10). Enzyme + Base only.
#
# Reduced from SpeedyWeather.jl's spectral transform (parameter AD w.r.t. the model); same
# type-analysis failure class as https://github.com/EnzymeAD/Enzyme.jl/issues/3275, second shape.
#
# Ingredients (each ablation-verified necessary):
#  1. The struct is loaded from a MUTABLE Duplicated wrapper inside the differentiated function.
#     Identical content behind an immutable wrapper succeeds (case 2), as does Const(wrapper).
#  2. The struct carries a large nested INACTIVE field: a NamedTuple of NSPECTRA structs, each
#     holding just TWO Ints + ONE Vector{Int}. It is never read except for ONE Int (the bound).
#     Every part of that shape is load-bearing:
#       - count: the threshold is EXACTLY 21 -> 22 (21 such structs pass, 22 fail; bisected on
#         Julia 1.12.6 / Enzyme 0.13.173 — bump NSPECTRA for margin if a future version shifts
#         it; with 4 Vector fields per struct instead, 8 pass and 16 fail);
#       - BOTH leading Ints: with a single leading Int (Vector at byte offset 8 instead of 16 —
#         the exact size of the un-deducible copy) the same shape passes even at double count;
#       - the struct nesting itself: 128 FLAT Vector fields (~3 KB, no struct wrapper) pass,
#         so it is NOT a raw-byte threshold.
#  3. The differentiated loop's trip count derives from that Int, LOADED out of the Duplicated
#     chain at runtime; with a compile-time-constant bound the identical body succeeds.
#     The body is a plain nested loop accumulating a Matrix field into the Duplicated output
#     (a Vector field instead of the Matrix, or a single-element body, is below the floor).
#
# Verified (Enzyme 0.13.173):
#   Julia 1.12.6:  mutable wrapper -> EnzymeNoTypeError ("cannot deduce type of copy",
#                  24-byte Vector aggregate); immutable wrapper -> exact gradient (40).
#   Julia 1.10.11: both succeed with the exact gradient.
using Enzyme

println("Julia ", VERSION, ", Enzyme ", pkgversion(Enzyme))

struct Spectrum{V}             # mimics SpeedyWeather's Spectrum: 2 Ints + Vector field(s)
    lmax::Int
    mmax::Int
    l_indices::V
end
struct Transform{L, G}
    legendre::L                # Matrix{Float32}: the data the loop differentiates through
    gradients::G               # NamedTuple of NSPECTRA Spectrum — inactive, unread except g1.lmax
end
mutable struct MutModel{T}     # FAILS on Julia >= 1.11
    transform::T
end
struct ImmModel{T}             # identical content, SUCCEEDS
    transform::T
end

# exact threshold (bisected, Julia 1.12.6 / Enzyme 0.13.173): 21 spectra pass, 22 fail
const NSPECTRA = 22

spectrum() = Spectrum(10, 10, collect(1:8))
transform() = Transform(rand(Float32, 8, 8), NamedTuple{ntuple(i -> Symbol(:g, i), NSPECTRA)}(ntuple(_ -> spectrum(), NSPECTRA)))

function walk!(out, S)
    for j in 1:(S.gradients.g1.lmax - 2)    # trip count loaded from the Duplicated chain (= 8)
        for l in 1:5
            out[l] += S.legendre[l, j]
        end
    end
    return nothing
end
loss!(out, w) = walk!(out, w.transform)

function try_case(name, w)
    out = zeros(Float32, 8)
    dout = make_zero(out)
    dout .= 1
    dw = make_zero(w)
    try
        autodiff(Reverse, loss!, Const, Duplicated(out, dout), Duplicated(w, dw))
        println(rpad(name, 20), ": SUCCESS  (sum |d legendre| = ", sum(abs, dw.transform.legendre), ", expected 40.0)")
    catch e
        println(rpad(name, 20), ": FAILED -> ", first(sprint(showerror, e), 100))
    end
    return flush(stdout)
end

try_case("mutable wrapper", MutModel(transform()))
try_case("immutable wrapper", ImmModel(transform()))
