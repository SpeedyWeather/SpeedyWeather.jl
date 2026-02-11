@testset "MatrixSpectralTransform: Initialization and roundtrip" begin
    @testset for trunc in (31, 63)
        @testset for NF in (Float32, Float64)
            @testset for Grid in (
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                    FullClenshawGrid,
                    OctahedralClenshawGrid,
                )

                spectrum = Spectrum(trunc)
                grid = Grid(SpeedyTransforms.get_nlat_half(trunc))
                M = MatrixSpectralTransform(spectrum, grid; NF)

                # initialization checks
                @test M isa MatrixSpectralTransform
                @test eltype(M) == NF
                @test M.nlayers > 0
                @test size(M.forward, 1) == LowerTriangularArrays.nonzeros(spectrum)
                @test size(M.forward, 2) == RingGrids.get_npoints(grid)
                @test size(M.backward, 1) == RingGrids.get_npoints(grid)
                @test size(M.backward, 2) == LowerTriangularArrays.nonzeros(spectrum)
                @test size(M.backward_real) == size(M.backward)
                @test size(M.backward_imag) == size(M.backward)

                tol = NF == Float32 ? 1.0e-3 : 1.0e-7

                # 2D roundtrip: start in spectral space to avoid aliasing issues
                spec = randn(Complex{NF}, spectrum)
                field = transform(spec, M)
                spec_roundtrip = transform(field, M)
                field_roundtrip = transform(spec_roundtrip, M)
                @test field_roundtrip ≈ field atol = tol rtol = tol

                # 2D roundtrip: start in grid space
                field2 = randn(NF, grid)
                spec2 = transform(field2, M)
                field2_roundtrip = transform(spec2, M)
                spec2_roundtrip = transform(field2_roundtrip, M)
                @test spec2_roundtrip ≈ spec2 atol = tol rtol = tol

                # 3D roundtrip: start in spectral space
                nlayers = 8
                M3D = MatrixSpectralTransform(spectrum, grid; NF, nlayers)

                spec3D = randn(Complex{NF}, spectrum, nlayers)
                field3D = transform(spec3D, M3D)
                spec3D_roundtrip = transform(field3D, M3D)
                field3D_roundtrip = transform(spec3D_roundtrip, M3D)
                @test field3D_roundtrip ≈ field3D atol = tol rtol = tol

                # 3D roundtrip: start in grid space
                field3D_2 = randn(NF, grid, nlayers)
                spec3D_2 = transform(field3D_2, M3D)
                field3D_2_roundtrip = transform(spec3D_2, M3D)
                spec3D_2_roundtrip = transform(field3D_2_roundtrip, M3D)
                @test spec3D_2_roundtrip ≈ spec3D_2 atol = tol rtol = tol
            end
        end
    end
end

@testset "MatrixSpectralTransform: Agreement with SpectralTransform" begin
    @testset for trunc in (31, 63)
        @testset for NF in (Float32, Float64)
            @testset for Grid in (
                    FullGaussianGrid,
                    OctahedralGaussianGrid,
                )

                nlayers = 8
                spectrum = Spectrum(trunc)
                grid = Grid(SpeedyTransforms.get_nlat_half(trunc))
                S = SpectralTransform(spectrum, grid; NF, nlayers)
                M = MatrixSpectralTransform(spectrum, grid; NF, nlayers)

                tol = NF == Float32 ? 1.0e-3 : 1.0e-7

                # 2D spectral -> grid
                spec = randn(Complex{NF}, spectrum)
                field_S = transform(spec, S)
                field_M = transform(spec, M)
                @test field_M ≈ field_S atol = tol rtol = tol

                # 2D grid -> spectral
                field = randn(NF, grid)
                spec_S = transform(field, S)
                spec_M = transform(field, M)
                @test spec_M ≈ spec_S atol = tol rtol = tol

                # 3D spectral -> grid
                spec3D = randn(Complex{NF}, spectrum, nlayers)
                field3D_S = transform(spec3D, S)
                field3D_M = transform(spec3D, M)
                @test field3D_M ≈ field3D_S atol = tol rtol = tol

                # 3D grid -> spectral
                field3D = randn(NF, grid, nlayers)
                spec3D_S = transform(field3D, S)
                spec3D_M = transform(field3D, M)
                @test spec3D_M ≈ spec3D_S atol = tol rtol = tol
            end
        end
    end
end
