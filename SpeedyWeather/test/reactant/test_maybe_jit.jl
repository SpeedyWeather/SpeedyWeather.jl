# Tests for @maybe_jit nesting behaviour.
# Nested @maybe_jit must detect within_compile() and skip the inner @jit.

@testset "@maybe_jit nesting" begin

    # -----------------------------------------------------------------------
    # 1. Simple function
    # -----------------------------------------------------------------------

    function _nested_inner!(arch, x)
        x .+= one(eltype(x))   # non-idempotent: double-call gives 2, single gives 1
        return nothing
    end

    function _nested_outer!(arch, x)
        @maybe_jit arch _nested_inner!(arch, x)   # nested @maybe_jit
        return nothing
    end

    @testset "simple - CPU" begin
        x = zeros(Float32, 4)
        @maybe_jit CPU() _nested_outer!(CPU(), x)
        @test all(x .== 1.0f0)
    end

    @testset "simple - Reactant" begin
        arch = ReactantDevice()
        x = Reactant.to_rarray(zeros(Float32, 4))
        @maybe_jit arch _nested_outer!(arch, x)
        @test all(Array(x) .== 1.0f0)   # must be 1, not 2
    end

    # -----------------------------------------------------------------------
    # 2. transform! with MatrixSpectralTransform
    # -----------------------------------------------------------------------

    NF = Float32
    trunc = TRUNC
    spectrum = Spectrum(trunc)
    grid = OctahedralGaussianGrid(SpeedyTransforms.get_nlat_half(trunc))
    M_cpu = MatrixSpectralTransform(spectrum, grid; NF)

    # random spectral coefficients as input
    spec_cpu = randn(Complex{NF}, spectrum)
    field_cpu = zeros(NF, grid)
    scratch_cpu = copy(M_cpu.scratch_memory)

    function _transform_nested!(arch, field, coeffs, scratch, M)
        @maybe_jit arch transform!(field, coeffs, scratch, M)
        return nothing
    end

    @testset "transform! - CPU" begin
        _transform_nested!(CPU(), field_cpu, spec_cpu, scratch_cpu, M_cpu)
        @test all(isfinite, Array(field_cpu))
    end

    @testset "transform! - Reactant" begin
        arch = ReactantDevice()
        spectrum_r = on_architecture(arch, spectrum)
        grid_r = on_architecture(arch, grid)
        M_rea = MatrixSpectralTransform(spectrum_r, grid_r; NF)

        # move CPU arrays to Reactant via on_architecture
        spec_rea = on_architecture(arch, copy(spec_cpu))
        field_rea = on_architecture(arch, zeros(NF, grid))
        scratch_rea = on_architecture(arch, copy(M_cpu.scratch_memory))

        @maybe_jit arch _transform_nested!(arch, field_rea, spec_rea, scratch_rea, M_rea)

        @test isapprox(Array(field_cpu), Array(field_rea); rtol = RTOL, atol = ATOL)
    end

end
