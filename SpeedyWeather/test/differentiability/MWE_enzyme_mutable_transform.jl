# MWE — error #2: EnzymeNoTypeError from loading a SpectralTransform out of a MUTABLE
# Duplicated struct and passing it through transform! (Julia >= 1.11, Enzyme 0.13.173).
#
# This is the reduced form of differentiating SpeedyWeather w.r.t. the MODEL
# (`Duplicated(model)`, parameter sensitivities): models are `mutable struct`s, and
# `vorticity_flux_curldiv!` does `S = model.spectral_transform; transform!(..., S)`.
#
#   julia +1.12 --project=SpeedyWeather --check-bounds=yes <this file>
#     m1 : mutable wrapper + transform!    -> EnzymeNoTypeError ("cannot deduce type of copy",
#                                             a 24-byte Vector aggregate inside S)
#     m1b: mutable wrapper + broadcast     -> SUCCESS (plain use of S is fine)
#     i1 : immutable wrapper + transform!  -> SUCCESS (mutability is the trigger; even an
#          immutable wrapper carrying ALL 18 fields of the real model passes)
#   On Julia 1.10 all three succeed.
#
# Both transform! halves fail independently under m1: _fourier! (custom EnzymeRules, FFTW)
# AND _legendre! (plain KernelAbstractions kernels) — the custom rules are NOT required.
# Workaround: Enzyme.API.looseTypeAnalysis!(true) (enabled by SpeedyWeatherEnzymeExt on
# Julia >= 1.11; gradients verified against FD and Julia 1.10).
# Related upstream issue (same type-analysis class): https://github.com/EnzymeAD/Enzyme.jl/issues/3275
using SpeedyWeather, Enzyme

Enzyme.API.looseTypeAnalysis!(false)    # undo the ext's workaround to expose the bug

spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)
model = BarotropicModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid))
sim = initialize!(model)
scratch0 = sim.variables.scratch.transform_memory
S0 = model.spectral_transform

println("Julia ", VERSION, ", Enzyme ", pkgversion(Enzyme))

mutable struct MutModel{ST}
    spectral_transform::ST
end
struct ImmModel{ST}
    spectral_transform::ST
end

ft!(spec, grid, scratch, w) = (S = w.spectral_transform; SpeedyWeather.transform!(spec, grid, scratch, S); nothing)
fb!(spec, grid, scratch, w) = (S = w.spectral_transform; spec.data[1:8] .+= complex.(S.gradients.eigenvalues[1:8]); nothing)

function try_case(name, fn, w)
    spec = rand(ComplexF32, spectral_grid.spectrum)
    grid = rand(Float32, spectral_grid.grid)
    dspec = make_zero(spec); dspec .= 1 + im
    dgrid = make_zero(grid)
    scratch = deepcopy(scratch0)
    try
        autodiff(
            set_runtime_activity(Reverse), fn, Const,
            Duplicated(spec, dspec), Duplicated(grid, dgrid),
            Duplicated(scratch, make_zero(scratch)), Duplicated(w, make_zero(w)),
        )
        println(rpad(name, 38), ": SUCCESS")
    catch e
        println(rpad(name, 38), ": FAILED -> ", first(sprint(showerror, e), 120))
    end
end

try_case("m1 : mutable wrapper + transform!", ft!, MutModel(S0))
try_case("m1b: mutable wrapper + broadcast", fb!, MutModel(S0))
try_case("i1 : immutable wrapper + transform!", ft!, ImmModel(S0))
