# Standalone run of the new "chunked transform Enzyme rules" testset content, in the
# differentiability project (which resolves the deps). Confirms it passes on 1.10 AND 1.12.
# Run: julia --project=SpeedyWeather/test/differentiability [--check-bounds=yes] <file>
using SpeedyWeather
using Enzyme
using Test
import SpeedyWeather.SpeedyTransforms as ST
import SpeedyWeather.SpeedyTransforms: SpectralTransform, transform!, Spectrum

println("Julia ", VERSION, ", Enzyme ", pkgversion(Enzyme))

@testset "chunked transform Enzyme rules" begin
    trunc = 5
    spectrum = Spectrum(trunc, one_degree_more = true)
    grid = FullGaussianGrid(ST.get_nlat_half(trunc, 2))
    NL = 4
    S_chunked = SpectralTransform(spectrum, grid; NF = Float32, nlayers = NL, transform_batch = [1])
    S_batched = SpectralTransform(spectrum, grid; NF = Float32, nlayers = NL, transform_batch = [1, NL])
    @test ST._needs_chunking(NL, S_chunked)
    @test !ST._needs_chunking(NL, S_batched)

    coeffs0 = rand(ComplexF32, spectrum, NL)
    dfield0 = rand(Float32, grid, NL)
    dcoeffs = map((S_chunked, S_batched)) do S
        coeffs = deepcopy(coeffs0)
        field = zeros(Float32, grid, NL)
        dfield = deepcopy(dfield0)
        dc = make_zero(coeffs)
        autodiff(
            set_runtime_activity(Reverse), transform!, Const,
            Duplicated(field, dfield), Duplicated(coeffs, dc),
            Duplicated(deepcopy(S.scratch_memory), make_zero(deepcopy(S.scratch_memory))), Const(S),
        )
        dc
    end
    @test all(isfinite, dcoeffs[1].data)
    @test any(!iszero, dcoeffs[1].data)
    @test isapprox(dcoeffs[1].data, dcoeffs[2].data, rtol = 1.0e-4)

    field0 = rand(Float32, grid, NL)
    dcoeffs0 = rand(ComplexF32, spectrum, NL)
    dfields = map((S_chunked, S_batched)) do S
        field = deepcopy(field0)
        coeffs = zeros(ComplexF32, spectrum, NL)
        dcoeffs_seed = deepcopy(dcoeffs0)
        df = make_zero(field)
        autodiff(
            set_runtime_activity(Reverse), transform!, Const,
            Duplicated(coeffs, dcoeffs_seed), Duplicated(field, df),
            Duplicated(deepcopy(S.scratch_memory), make_zero(deepcopy(S.scratch_memory))), Const(S),
        )
        df
    end
    @test all(isfinite, dfields[1].data)
    @test any(!iszero, dfields[1].data)
    @test isapprox(dfields[1].data, dfields[2].data, rtol = 1.0e-4)
end
