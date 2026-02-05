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
