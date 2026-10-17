# The MatrixSpectralTransform is also called outside of compiled code, e.g. for the initial
# conditions, so it has to work eagerly on Reactant arrays in both directions

@testset "MatrixSpectralTransform on Reactant outside of compilation" begin
    arch = ReactantDevice()
    nlayers = 2
    spectral_grid = SpectralGrid(; architecture = arch, truncation = TRUNCATION, nlayers)
    spectral_grid_cpu = SpectralGrid(; truncation = TRUNCATION, nlayers)
    M = MatrixSpectralTransform(spectral_grid)
    M_cpu = MatrixSpectralTransform(spectral_grid_cpu)

    # grid to spectral
    field_cpu = rand(Float32, spectral_grid_cpu.grid, nlayers)
    coeffs_cpu = transform(field_cpu, M_cpu)
    coeffs = transform(on_architecture(arch, field_cpu), M)
    @test on_architecture(SpeedyWeather.CPU(), coeffs.data) ≈ coeffs_cpu.data

    # spectral to grid
    field_back_cpu = transform(coeffs_cpu, M_cpu)
    field_back = transform(coeffs, M)
    @test on_architecture(SpeedyWeather.CPU(), field_back.data) ≈ field_back_cpu.data
end

# The FFT-based SpectralTransform on Reactant: no FFT plans, Fourier transform per group of rings of
# equal length, Legendre transform as a batched matrix multiply (SpeedyTransformsReactantExt).
# Compared with the CPU transform eagerly (outside of compilation) and compiled, on views into a
# fused parent array as the model's variables are. The grids cover rings in north/south pairs
# (octahedral), all rings in one FFT (full grid) and an odd number of rings with an equator (HEALPix).
@testset "SpectralTransform on Reactant" begin
    arch = ReactantDevice()
    nlayers = 3
    relative_error(a, b) = maximum(abs, Array(a) - Array(b)) / maximum(abs, Array(b))

    for Grid in (OctahedralGaussianGrid, FullGaussianGrid, HEALPixGrid)
        spectral_grid = SpectralGrid(; architecture = arch, truncation = TRUNCATION, nlayers, Grid)
        spectral_grid_cpu = SpectralGrid(; truncation = TRUNCATION, nlayers, Grid)
        S = SpectralTransform(spectral_grid)
        S_cpu = SpectralTransform(spectral_grid_cpu)
        @test S.gpu_graphs == false             # XLA runs a compiled program as one graph itself
        @test isnothing(S.rfft_plans_batched)   # no FFT plans, XLA plans the FFTs
        # only the dense copy of the polynomials is used, the others are empty (`length`, not `isempty`: on a
        # concrete Reactant array `isempty` checks whether its buffer was deleted, not its length)
        @test size(S.legendre_polynomials.data, 2) == 0
        @test length(S.legendre_polynomials_transposed) == 0

        # no layer dimension (vectors, as single-layer variables in the model), one layer, several layers
        for layers in ((), (1,), (nlayers,))
            field_cpu = rand(Float32, spectral_grid_cpu.grid, layers...)
            coeffs_cpu = transform(field_cpu, S_cpu)
            field_back_cpu = transform(coeffs_cpu, S_cpu; unscale_coslat = true)

            # outside of compilation
            coeffs = transform(on_architecture(arch, field_cpu), S)
            field_back = transform(coeffs, S; unscale_coslat = true)
            @test relative_error(coeffs.data, coeffs_cpu.data) < 1.0e-5
            @test relative_error(field_back.data, field_back_cpu.data) < 1.0e-5

            # compiled, both directions in one program, input a view into a fused parent array
            K = prod(layers; init = 1)
            parent_cpu = zeros(Float32, spectral_grid_cpu.grid, K + 2)
            parent_cpu.data[:, 2:(K + 1)] .= reshape(field_cpu.data, :, K)
            parent = on_architecture(arch, parent_cpu)
            columns = isempty(layers) ? 2 : (2:(K + 1))
            field_view = Field(view(parent.data, :, columns), spectral_grid.grid)
            coeffs2 = zeros(ComplexF32, spectral_grid.spectrum, layers...)
            field2 = zeros(Float32, spectral_grid.grid, layers...)
            function roundtrip!(coeffs, field_out, field_in, S)
                transform!(coeffs, field_in, S)
                transform!(field_out, coeffs, S; unscale_coslat = true)
                return nothing
            end
            compiled_roundtrip! = @compile compile_options = SpeedyWeather.reactant_compile_options() roundtrip!(coeffs2, field2, field_view, S)
            compiled_roundtrip!(coeffs2, field2, field_view, S)
            @test relative_error(coeffs2.data, coeffs_cpu.data) < 1.0e-5
            @test relative_error(field2.data, field_back_cpu.data) < 1.0e-5
        end
    end
end
