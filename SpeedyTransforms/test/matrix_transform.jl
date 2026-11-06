@testset "MatrixSpectralTransform: Initialization and roundtrip" begin
    @testset for truncation in (32, 64)
        @testset for NF in (Float32, Float64)
            @testset for Grid in (
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                )

                spectrum = Spectrum(truncation)
                grid = Grid(SpeedyTransforms.get_nlat_half(truncation))
                M = MatrixSpectralTransform(spectrum, grid; NF)

                # initialization checks
                @test M isa MatrixSpectralTransform
                @test eltype(M) == NF
                @test M.nlayers > 0
                nharmonics = LowerTriangularArrays.nonzeros(spectrum)
                npoints = RingGrids.get_npoints(grid)
                @test eltype(M.forward_stacked) == NF            # real matrices, stacked Re/Im parts
                @test size(M.forward_stacked) == (2nharmonics, npoints)
                @test size(M.backward_stacked) == (npoints, 2nharmonics)
                @test size(M.scratch_memory) == (2nharmonics, M.nlayers)

                rtol = NF == Float32 ? 1.0e-3 : 1.0e-7
                atol = NF == Float32 ? 1.0e-3 : 1.0e-7

                # 2D roundtrip: start in spectral space to avoid aliasing issues
                spec = randn(Complex{NF}, spectrum)
                field = transform(spec, M)
                spec_roundtrip = transform(field, M)
                field_roundtrip = transform(spec_roundtrip, M)
                @test field_roundtrip ≈ field atol = atol rtol = rtol

                # 2D roundtrip: start in grid space
                field2 = randn(NF, grid)
                spec2 = transform(field2, M)
                field2_roundtrip = transform(spec2, M)
                spec2_roundtrip = transform(field2_roundtrip, M)
                @test spec2_roundtrip ≈ spec2 atol = atol rtol = rtol

                # 3D roundtrip: start in spectral space
                nlayers = 8
                M3D = MatrixSpectralTransform(spectrum, grid; NF, nlayers)

                spec3D = randn(Complex{NF}, spectrum, nlayers)
                field3D = transform(spec3D, M3D)
                spec3D_roundtrip = transform(field3D, M3D)
                field3D_roundtrip = transform(spec3D_roundtrip, M3D)
                @test field3D_roundtrip ≈ field3D atol = atol rtol = rtol

                # 3D roundtrip: start in grid space
                field3D_2 = randn(NF, grid, nlayers)
                spec3D_2 = transform(field3D_2, M3D)
                field3D_2_roundtrip = transform(spec3D_2, M3D)
                spec3D_2_roundtrip = transform(field3D_2_roundtrip, M3D)
                @test spec3D_2_roundtrip ≈ spec3D_2 atol = atol rtol = rtol
            end
        end
    end
end

@testset "MatrixSpectralTransform: Agreement with SpectralTransform" begin
    @testset for truncation in (32, 64)
        @testset for NF in (Float32, Float64)
            @testset for Grid in (
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                    FullClenshawGrid,
                    OctahedralClenshawGrid,
                    HEALPixGrid,
                    OctaminimalGaussianGrid,
                    OctaHEALPixGrid,
                )

                nlayers = 8
                spectrum = Spectrum(truncation)
                grid = Grid(SpeedyTransforms.get_nlat_half(truncation))
                S = SpectralTransform(spectrum, grid; NF, nlayers)
                M = MatrixSpectralTransform(spectrum, grid; NF, nlayers)

                rtol = NF == Float32 ? 1.0e-3 : 1.0e-7
                atol = NF == Float32 ? 1.0e-3 : 1.0e-7

                # 2D spectral -> grid
                spec = randn(Complex{NF}, spectrum)
                field_S = transform(spec, S)
                field_M = transform(spec, M)
                @test field_M ≈ field_S atol = atol rtol = rtol

                # 2D grid -> spectral
                field = randn(NF, grid)
                spec_S = transform(field, S)
                spec_M = transform(field, M)
                @test spec_M ≈ spec_S atol = atol rtol = rtol

                # 3D spectral -> grid
                spec3D = randn(Complex{NF}, spectrum, nlayers)
                field3D_S = transform(spec3D, S)
                field3D_M = transform(spec3D, M)
                @test field3D_M ≈ field3D_S atol = atol rtol = rtol

                # 3D grid -> spectral
                field3D = randn(NF, grid, nlayers)
                spec3D_S = transform(field3D, S)
                spec3D_M = transform(field3D, M)
                @test spec3D_M ≈ spec3D_S atol = atol rtol = rtol
            end
        end
    end
end

@testset "MatrixSpectralTransform: wide batches and views" begin
    NF = Float32
    spectrum = Spectrum(32)
    grid = OctahedralGaussianGrid(SpeedyTransforms.get_nlat_half(32))
    nlayers = 8
    S = SpectralTransform(spectrum, grid; NF, nlayers)
    M = MatrixSpectralTransform(spectrum, grid; NF, nlayers)

    spec = randn(Complex{NF}, spectrum, nlayers)
    field = transform(spec, S)

    # a batch as wide as the model's tendency batch (9nlayers+1 columns) is a single multiply when
    # the transform (its scratch) is constructed for it; a too narrow scratch throws
    K = 9nlayers + 1
    M_wide = MatrixSpectralTransform(spectrum, grid; NF, nlayers = K)
    spec_wide = randn(Complex{NF}, spectrum, K)
    field_wide = transform(spec_wide, S)
    @test transform(field_wide, M_wide) ≈ transform(field_wide, S) rtol = 1.0e-3
    @test transform(spec_wide, M_wide) ≈ field_wide rtol = 1.0e-3
    nharmonics = LowerTriangularArrays.nonzeros(spectrum)
    scratch_narrow = zeros(NF, 2nharmonics, 3)
    @test_throws DimensionMismatch transform!(spec, field, scratch_narrow, M)
    @test_throws DimensionMismatch transform!(field, spec, scratch_narrow, M)

    # transforms into and out of views (2D view of a 3D parent, 1D slot) as in the model's fused
    # variables agree with plain arrays and land in the parent
    nsteps = 2
    coeffs_parent = zeros(Complex{NF}, spectrum, nlayers, nsteps)
    field_parent = zeros(NF, grid, nlayers, nsteps)
    coeffs_view = LowerTriangularArray(view(coeffs_parent.data, :, :, 2), spectrum)
    field_view = Field(view(field_parent.data, :, :, 2), grid)
    coeffs_view .= spec
    transform!(field_view, coeffs_view, M)
    @test field_view ≈ field rtol = 1.0e-3
    @test field_parent[:, :, 2] ≈ field.data rtol = 1.0e-3
    @test all(iszero, field_parent[:, :, 1])

    coeffs_view .= 0
    transform!(coeffs_view, field_view, M)
    @test coeffs_view ≈ transform(field, M) rtol = 1.0e-3
    @test all(iszero, coeffs_parent[:, :, 1])

    # single layer slot (1D view) in both directions
    coeffs_slot = LowerTriangularArray(view(coeffs_parent.data, :, 3, 2), spectrum)
    field_slot = Field(view(field_parent.data, :, 3, 2), grid)
    transform!(field_slot, coeffs_slot, M)
    @test field_slot ≈ field[:, 3] rtol = 1.0e-3
    transform!(coeffs_slot, field_slot, M)
    @test coeffs_slot ≈ transform(field, M)[:, 3] rtol = 1.0e-3
end
