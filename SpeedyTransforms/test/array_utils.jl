@testset "similar for LowerTriangularArray and Field" begin
    @testset for NF in (Float32, Float64)
        @testset for Grid in (FullGaussianGrid, OctahedralGaussianGrid)
            trunc = 31
            nlat_half = 9
            nlayers = 3

            # Create spectrum and grid
            spectrum = Spectrum(trunc)
            grid = Grid(nlat_half)

            # Create test LowerTriangularArray and Field
            spec = randn(Complex{NF}, spectrum, nlayers)
            field = randn(NF, grid, nlayers)

            @testset "similar(LowerTriangularArray, grid)" begin
                result = similar(spec, grid, NF)

                # Check type
                @test result isa Field
                @test eltype(result) == NF

                # Check dimensions
                @test size(result, 1) == RingGrids.get_npoints(grid)
                @test size(result, 2) == nlayers

                # Check grid is preserved
                @test result.grid === grid
            end

            @testset "similar(Field, spectrum)" begin
                result = similar(field, spectrum, Complex{NF})

                # Check type
                @test result isa LowerTriangularArray
                @test eltype(result) == Complex{NF}

                # Check dimensions
                @test size(result, 1) == LowerTriangularArrays.nonzeros(spectrum)
                @test size(result, 2) == nlayers

                # Check spectrum is preserved
                @test result.spectrum === spectrum
            end
        end
    end
end

@testset "wrapped_view for LowerTriangularArray and Field" begin
    @testset for NF in (Float32, Float64)
        @testset for Grid in (FullGaussianGrid, OctahedralGaussianGrid)
            trunc = 31
            nlat_half = 9
            nlayers = 3

            # Create spectrum and grid
            spectrum = Spectrum(trunc)
            grid = Grid(nlat_half)

            @testset "wrapped_view on Field with `:, i` indexing" begin
                # Create a 3D field
                field = randn(NF, grid, nlayers)

                # wrapped_view should return a Field with SubArray data
                view_result = wrapped_view(field, :, 1)

                @test view_result isa Field
                @test view_result.data isa SubArray
                @test eltype(view_result) == NF
                @test size(view_result, 1) == RingGrids.get_npoints(grid)
                @test view_result.grid === grid

                # Check that view shares data with original
                @test view_result[1] ≈ field[1, 1]
                field[1, 1] = NF(42.0)
                @test view_result[1] ≈ NF(42.0)
            end

            @testset "wrapped_view on 2D Field with `:` indexing" begin
                # Create a 2D field
                field = randn(NF, grid)

                # wrapped_view should return a Field with SubArray data
                view_result = wrapped_view(field, :)

                @test view_result isa Field
                @test view_result.data isa SubArray
                @test eltype(view_result) == NF
                @test size(view_result) == size(field)
                @test view_result.grid === grid
            end

            @testset "wrapped_view on Field with regular indexing (fallback)" begin
                field = randn(NF, grid, nlayers)

                # Regular indexing should fall back to normal view (SubArray)
                view_result = wrapped_view(field, 1)

                @test view_result isa SubArray
            end

            @testset "wrapped_view on LowerTriangularArray with `:, i` indexing" begin
                # Create a 3D LowerTriangularArray
                lta = randn(Complex{NF}, spectrum, nlayers)

                # wrapped_view should return a LowerTriangularArray with SubArray data
                view_result = wrapped_view(lta, :, 1)

                @test view_result isa LowerTriangularArray
                @test view_result.data isa SubArray
                @test eltype(view_result) == Complex{NF}
                @test view_result.spectrum === lta.spectrum

                # Check that view shares data with original
                original_val = lta[1]
                view_result[1] = original_val * Complex{NF}(2.0)
                @test lta[1] ≈ original_val * Complex{NF}(2.0)
            end

            @testset "wrapped_view on LowerTriangularMatrix with `:` indexing" begin
                # Create a 2D LowerTriangularArray (matrix)
                lta = randn(Complex{NF}, spectrum)

                # wrapped_view should return a LowerTriangularArray with SubArray data
                view_result = wrapped_view(lta, :)

                @test view_result isa LowerTriangularArray
                @test view_result.data isa SubArray
                @test eltype(view_result) == Complex{NF}
                @test view_result.spectrum === lta.spectrum
            end

            @testset "wrapped_view on LowerTriangularArray with regular indexing (fallback)" begin
                lta = randn(Complex{NF}, spectrum, nlayers)

                # Regular indexing should fall back to normal view (SubArray)
                view_result = wrapped_view(lta, 1)

                @test view_result isa SubArray
            end

            @testset "wrapped_view preserves data sharing" begin
                # Test that wrapped_view creates actual views, not copies
                field = randn(NF, grid, nlayers)
                view_result = wrapped_view(field, :, 2)

                # Modify through view
                view_result[1] = NF(999.0)

                # Check original is modified
                @test field[1, 2] ≈ NF(999.0)

                # And vice versa
                field[2, 2] = NF(-999.0)
                @test view_result[2] ≈ NF(-999.0)
            end
        end
    end
end
