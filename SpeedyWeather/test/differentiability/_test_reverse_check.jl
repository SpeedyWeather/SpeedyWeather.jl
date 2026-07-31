# Focused check: does test_reverse (EnzymeTestUtils) work for the chunked transform! rule via
# nothing-returning wrappers? (The direct transform! errored in ETU's FD — returned array.)
# Run: julia --project=SpeedyWeather/test/differentiability [--check-bounds=yes] <file>
using SpeedyWeather
using EnzymeTestUtils, Enzyme, FiniteDifferences, Test
import EnzymeTestUtils: test_approx, test_reverse
import AbstractFFTs
import SpeedyWeather.SpeedyTransforms as ST
import SpeedyWeather.SpeedyTransforms: SpectralTransform, transform!, Spectrum

println("Julia ", VERSION, ", Enzyme ", pkgversion(Enzyme))

# FFT-plan test_approx override (copied from spectral_transform_ad_rules.jl)
function EnzymeTestUtils.test_approx(x::AbstractFFTs.Plan, y::AbstractFFTs.Plan, msg; kwargs...)
    EnzymeTestUtils.@test_msg "$msg: types must match" typeof(x) == typeof(y)
    names = fieldnames(typeof(x))[1:(end - 1)]
    if isempty(names)
        EnzymeTestUtils.@test_msg msg x == y
    else
        for k in names
            k isa Symbol && hasproperty(x, k) || continue
            EnzymeTestUtils.test_approx(getfield(x, k), getfield(y, k), "$msg.$k"; kwargs...)
        end
    end
    return nothing
end

@testset "test_reverse chunked transform" begin
    trunc = 5
    spectrum = Spectrum(trunc, one_degree_more = true)
    grid = FullGaussianGrid(ST.get_nlat_half(trunc, 2))
    NL = 4
    S = SpectralTransform(spectrum, grid; NF = Float32, nlayers = NL, transform_batch = [1])

    transform_s2g!(field, coeffs, scratch, S) = (transform!(field, coeffs, scratch, S); nothing)
    transform_g2s!(coeffs, field, scratch, S) = (transform!(coeffs, field, scratch, S); nothing)

    let field = zeros(Float32, grid, NL), coeffs = rand(ComplexF32, spectrum, NL)
        test_reverse(
            transform_s2g!, Const,
            (field, Duplicated), (coeffs, Duplicated),
            (deepcopy(S.scratch_memory), Const), (S, Const);
            fdm = FiniteDifferences.central_fdm(5, 1), rtol = 1.0e-2, atol = 1.0e-2,
        )
    end
    let coeffs = zeros(ComplexF32, spectrum, NL), field = rand(Float32, grid, NL)
        test_reverse(
            transform_g2s!, Const,
            (coeffs, Duplicated), (field, Duplicated),
            (deepcopy(S.scratch_memory), Const), (S, Const);
            fdm = FiniteDifferences.central_fdm(5, 1), rtol = 1.0e-2, atol = 1.0e-2,
        )
    end
end
