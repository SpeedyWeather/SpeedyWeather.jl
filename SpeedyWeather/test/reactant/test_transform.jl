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

# The SpectralTransform on Reactant (SpeedyTransformsReactantExt), with the Fourier transform as truncated
# DFT (`fourier_dft = true`, the default) or with FFTs per group of rings of equal length, Legendre transform
# as a batched matrix multiply. Compared with the CPU transform (FFTW) eagerly (outside of compilation) and
# compiled, on views into a fused parent array as the model's variables are. The grids cover rings in
# north/south pairs (octahedral), all rings of equal length (full grid), an odd number of rings with an
# equator in a run of equal rings (HEALPix) and an equator ring on its own (OctaHEALPix, OctahedralClenshaw).
@testset "SpectralTransform on Reactant, fourier_dft = $fourier_dft" for fourier_dft in (true, false)
    arch = ReactantDevice()
    nlayers = 3
    relative_error(a, b) = maximum(abs, Array(a) - Array(b)) / maximum(abs, Array(b))

    for Grid in (OctahedralGaussianGrid, FullGaussianGrid, HEALPixGrid, OctaHEALPixGrid, OctahedralClenshawGrid)
        spectral_grid = SpectralGrid(; architecture = arch, truncation = TRUNCATION, nlayers, Grid)
        spectral_grid_cpu = SpectralGrid(; truncation = TRUNCATION, nlayers, Grid)
        S = SpectralTransform(spectral_grid; fourier_dft)
        S_cpu = SpectralTransform(spectral_grid_cpu)
        @test S.gpu_graphs == false             # XLA runs a compiled program as one graph itself
        @test isnothing(S.rfft_plans_batched)   # no FFT plans, XLA plans the FFTs
        # only the dense copy of the polynomials is used, the others are empty (`length`, not `isempty`: on a
        # concrete Reactant array `isempty` checks whether its buffer was deleted, not its length)
        @test size(S.legendre.polynomials.data, 2) == 0
        @test length(S.legendre.polynomials_transposed) == 0
        @test S.fourier_dft == fourier_dft
        (; mmax) = S.spectrum
        dft_size = fourier_dft ? 2mmax * maximum(S.nlons) * spectral_grid.grid.nlat_half : 0
        @test length(S.fourier_matrix_forward) == length(S.fourier_matrix_inverse) == dft_size
        # `show` reports the size of the data (not of the Reactant array wrapper) and the DFT only if used
        @test SpeedyWeather.SpeedyTransforms.nbytes(S.fourier_matrix_forward) == sizeof(Float32) * dft_size
        @test occursin("truncated DFT", sprint(show, S)) == fourier_dft
        @test size(S.scratch_memory.north, 2) == 0         # the transform needs no scratch memory on Reactant

        # no layer dimension (vectors, as single-layer variables in the model), one layer, several layers
        # (only several layers on the grids with a single equator ring, which only differ in the Fourier transform)
        layer_cases = Grid in (OctaHEALPixGrid, OctahedralClenshawGrid) ? ((nlayers,),) : ((), (1,), (nlayers,))
        for layers in layer_cases
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

# The Fourier transform as truncated DFT (`fourier_dft = true`) against the standard FFTs, for correctness and
# precision: its matrices and the traced Fourier stages (DFT and FFTs) ring by ring against the Fourier transform of
# the CPU SpectralTransform (FFTW), and the whole transform against the CPU SpectralTransform in Float64 (agrees to
# rounding) and Float32 (its error vs a Float64 reference not larger than that of the FFT transforms).
const SpeedyTransformsReactantExt = Base.get_extension(SpeedyWeather.SpeedyTransforms, :SpeedyTransformsReactantExt)

"""Fourier coefficients (nlat_half, mmax, 2, K) of the northern (1) and southern (2) rings of `data` from the Fourier
transform of the CPU SpectralTransform `S` (FFTW), weighted with the solid angles and rotated back by the longitude
offsets as its Legendre transform does, orders m = 0..mmax-1."""
function fftw_forward_spectra(data::AbstractMatrix, S)
    (; nfreq_max, lon_offsets, solid_angles) = S
    nlat_half, mmax, K = S.grid.nlat_half, S.spectrum.mmax, size(data, 2)
    f_north, f_south = zeros(ComplexF64, nfreq_max, K, nlat_half), zeros(ComplexF64, nfreq_max, K, nlat_half)
    SpeedyWeather.SpeedyTransforms._fourier!(f_north, f_south, Field(Float64.(data), S.grid), S)
    weights = reshape(solid_angles[1:nlat_half]' .* conj.(lon_offsets), mmax, 1, nlat_half)
    weighted(f) = permutedims(weights .* f[1:mmax, :, :], (3, 1, 2))     # (nlat_half, mmax, K)
    return cat(reshape(weighted(f_north), nlat_half, mmax, 1, K), reshape(weighted(f_south), nlat_half, mmax, 1, K); dims = 3)
end

"""Grid points (npoints, K) from the Fourier coefficients `spectra` (mmax, nlat_half, 2, K) of the northern (1) and
southern (2) rings, rotated by the longitude offsets as the Legendre transform of the CPU SpectralTransform `S` does,
with its inverse Fourier transform (FFTW)."""
function fftw_inverse_grid(spectra::AbstractArray{<:Complex, 4}, S)
    (; nfreq_max, lon_offsets) = S
    nlat_half, mmax, K = S.grid.nlat_half, S.spectrum.mmax, size(spectra, 4)
    g_north, g_south = zeros(ComplexF64, nfreq_max, K, nlat_half), zeros(ComplexF64, nfreq_max, K, nlat_half)
    g_north[1:mmax, :, :] .= permutedims(lon_offsets .* spectra[:, :, 1, :], (1, 3, 2))
    g_south[1:mmax, :, :] .= permutedims(lon_offsets .* spectra[:, :, 2, :], (1, 3, 2))
    field = zeros(Float64, S.grid, K)
    SpeedyWeather.SpeedyTransforms._fourier!(field, g_north, g_south, S)
    return field.data
end

@testset "Truncated DFT against FFTW" begin
    arch = ReactantDevice()
    K = 3
    relative_error(a, b) = maximum(abs, Array(a) .- Array(b)) / maximum(abs, Array(b))
    grids = (OctahedralGaussianGrid, FullGaussianGrid, HEALPixGrid, OctaHEALPixGrid, FullClenshawGrid, OctahedralClenshawGrid)

    # default: the DFT up to DFT_MAX_TRUNCATION on Reactant, never on other architectures
    (; DFT_MAX_TRUNCATION, default_fourier_dft) = SpeedyWeather.SpeedyTransforms
    @test default_fourier_dft(arch, Spectrum(truncation = DFT_MAX_TRUNCATION; architecture = arch))
    @test !default_fourier_dft(arch, Spectrum(truncation = DFT_MAX_TRUNCATION + 1; architecture = arch))
    @test !default_fourier_dft(SpeedyWeather.CPU(), Spectrum(truncation = 32))

    # the matrices on the host, ring by ring vs FFTW rfft/brfft (with the weights and rotations of the FFT path)
    @testset "DFT matrices, $Grid" for Grid in grids
        S = SpectralTransform(SpectralGrid(; truncation = TRUNCATION, nlayers = 1, Grid, NF = Float64))
        (; mmax) = S.spectrum
        forward, inverse = SpeedyTransformsReactantExt.fourier_matrices(S.nlons, S.lon_offsets, S.solid_angles, mmax)
        data = randn(sum(length, S.rings), K)
        spectra = randn(ComplexF64, mmax, S.grid.nlat_half, 2, K)
        spectra_fftw = fftw_forward_spectra(data, S)
        data_fftw = fftw_inverse_grid(spectra, S)
        for j in 1:S.grid.nlat_half, hemisphere in 1:2
            ring = hemisphere == 1 ? j : S.nlat - j + 1
            hemisphere == 2 && ring == j && continue
            nlon = S.nlons[j]
            parts = forward[:, 1:nlon, j] * data[S.rings[ring], :]
            @test complex.(parts[1:mmax, :], parts[(mmax + 1):end, :]) ≈ spectra_fftw[j, :, hemisphere, :] rtol = 1.0e-13
            ring_data = inverse[1:nlon, :, j] * vcat(real.(spectra[:, j, hemisphere, :]), imag.(spectra[:, j, hemisphere, :]))
            @test ring_data ≈ data_fftw[S.rings[ring], :] rtol = 1.0e-13
            @test all(iszero, inverse[(nlon + 1):end, :, j])   # nothing beyond the ring
            @test all(iszero, forward[:, (nlon + 1):end, j])
        end
    end

    # the traced Fourier stages on Reactant (the truncated DFT and the FFTs, which share their interface) vs FFTW,
    # ring by ring
    (; dft_fourier_forward, dft_fourier_inverse, fft_fourier_forward, fft_fourier_inverse) = SpeedyTransformsReactantExt
    @testset "Traced Fourier stages, $fourier_forward, $Grid, $NF" for (fourier_forward, fourier_inverse) in
            ((dft_fourier_forward, dft_fourier_inverse), (fft_fourier_forward, fft_fourier_inverse)),
            Grid in (OctahedralGaussianGrid, HEALPixGrid, OctahedralClenshawGrid), NF in (Float64, Float32)
        spectral_grid = SpectralGrid(; truncation = TRUNCATION, nlayers = K, Grid, NF, architecture = arch)
        S = SpectralTransform(spectral_grid)        # the DFT matrices are only built for the DFT
        S64 = SpectralTransform(SpectralGrid(; truncation = TRUNCATION, nlayers = K, Grid, NF = Float64))
        (; mmax) = S.spectrum
        nlat_half = spectral_grid.grid.nlat_half
        tolerance = NF == Float64 ? 1.0e-13 : 2.0e-6

        data = randn(NF, spectral_grid.npoints, K)
        parts = Array(@jit compile_options = SpeedyWeather.reactant_compile_options() fourier_forward(Reactant.to_rarray(data), S))
        @test size(parts) == (nlat_half, mmax, 2, 2, K)
        @test relative_error(complex.(parts[:, :, 1, :, :], parts[:, :, 2, :, :]), fftw_forward_spectra(data, S64)) < tolerance

        spectra = Complex{NF}.(randn(ComplexF64, mmax, nlat_half, 2, K))
        parts = cat(reshape(real.(spectra), mmax, nlat_half, 1, 2, K), reshape(imag.(spectra), mmax, nlat_half, 1, 2, K); dims = 3)
        data = Array(@jit compile_options = SpeedyWeather.reactant_compile_options() fourier_inverse(Reactant.to_rarray(parts), S))
        @test relative_error(data, fftw_inverse_grid(ComplexF64.(spectra), S64)) < tolerance
    end

    # the whole transform: DFT vs FFTW (CPU) and vs the FFT path on Reactant. Float64: equal up to rounding.
    # Float32: the error vs the Float64 FFTW transform not larger (up to a factor for the randomness) than FFTW's
    @testset "Transform with DFT vs FFTs, $Grid" for Grid in grids
        spectral_grid64 = SpectralGrid(; truncation = TRUNCATION, nlayers = K, Grid, NF = Float64)
        S64 = SpectralTransform(spectral_grid64)
        field64 = rand(Float64, spectral_grid64.grid, K)
        coeffs64 = transform(field64, S64)
        field_back64 = transform(coeffs64, S64; unscale_coslat = true)

        for NF in (Float64, Float32)
            spectral_grid_cpu = SpectralGrid(; truncation = TRUNCATION, nlayers = K, Grid, NF)
            spectral_grid = SpectralGrid(; truncation = TRUNCATION, nlayers = K, Grid, NF, architecture = arch)
            field_cpu = Field(NF.(field64.data), spectral_grid_cpu.grid)
            coeffs_cpu = LowerTriangularArray(Complex{NF}.(coeffs64.data), spectral_grid_cpu.spectrum)

            results = map(
                (
                    dft = SpectralTransform(spectral_grid; fourier_dft = true),
                    fft = SpectralTransform(spectral_grid; fourier_dft = false),
                    fftw = SpectralTransform(spectral_grid_cpu),
                )
            ) do S
                coeffs = transform(on_architecture(S.architecture, field_cpu), S)
                field = transform(on_architecture(S.architecture, coeffs_cpu), S; unscale_coslat = true)
                (on_architecture(SpeedyWeather.CPU(), coeffs.data), on_architecture(SpeedyWeather.CPU(), field.data))
            end

            if NF == Float64
                for direction in 1:2
                    @test relative_error(results.dft[direction], results.fftw[direction]) < 1.0e-13
                    @test relative_error(results.dft[direction], results.fft[direction]) < 1.0e-13
                end
            else
                for (direction, reference) in ((1, coeffs64.data), (2, field_back64.data))
                    error_dft = relative_error(results.dft[direction], reference)
                    error_fftw = relative_error(results.fftw[direction], reference)
                    error_fft = relative_error(results.fft[direction], reference)
                    @test error_dft < 1.0e-6
                    @test error_dft < 4 * max(error_fftw, error_fft)
                end
            end
        end
    end
end
