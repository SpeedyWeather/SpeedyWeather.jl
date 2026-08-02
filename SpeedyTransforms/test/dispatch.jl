using JET

# Regression guard for runtime dynamic dispatch in the spectral transform.

@testset "transform! free of runtime dispatch (JET)" begin
    arch = SpeedyTransforms.Architectures.CPU()
    spectrum = Spectrum(trunc = 31, architecture = arch)
    grid = FullGaussianGrid(RingGrids.get_nlat_half(32), arch)

    # nlayers=1 exercises the serial (K=1) plans; nlayers=8 the batched (K>1) plans. nlayers=1 also
    # leaves the batched Dict empty — it must still be concretely typed, otherwise the (compiled but
    # never-taken) batched branch of `_fourier!` reintroduces the abstract-plan cascade.
    for nlayers in (1, 8)
        S = SpectralTransform(spectrum, grid; NF = Float32, nlayers)
        specs = rand(ComplexF32, spectrum, nlayers)
        field = rand(Float32, grid, nlayers)
        scratch = S.scratch_memory

        # the FFT plans must be concretely typed (the root cause of the dispatch cascade)
        @test isconcretetype(eltype(S.rfft_plan_serial))
        @test isconcretetype(eltype(S.brfft_plan_serial))
        @test isconcretetype(eltype(valtype(S.rfft_plans_batched)))
        @test isconcretetype(eltype(valtype(S.brfft_plans_batched)))

        # and the transforms themselves must be dispatch-free
        @test_opt target_modules = (SpeedyTransforms,) transform!(specs, field, scratch, S)  # grid → spectral
        @test_opt target_modules = (SpeedyTransforms,) transform!(field, specs, scratch, S)  # spectral → grid

        # `curl!`/`divergence!` take `flipsign`/`add` as runtime `Bool`s but encode them as `KernelOP`
        # type parameters; without the union split in `_divergence!` the `KernelOP` widens to an
        # abstract type and `_divergence!` is reached by runtime dispatch.
        u = rand(ComplexF32, spectrum, nlayers)
        v = rand(ComplexF32, spectrum, nlayers)
        out = zeros(ComplexF32, spectrum, nlayers)
        for add in (false, true), flipsign in (false, true)
            @test_opt target_modules = (SpeedyTransforms,) curl!(out, u, v, S; add, flipsign)
            @test_opt target_modules = (SpeedyTransforms,) divergence!(out, u, v, S; add, flipsign)
        end
    end
end
