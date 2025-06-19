using KernelAbstractions
import SpeedyWeather: on_architecture, CPU, launch!
@testset "KernelAbstractions tests" begin 

    # To-Do write tests for each type of dims_type in the kernel launching util, 
    # the tests currently below will be removed when the KA becomes the only one

    # standard LTA + running index kernel 

    # Test the kernel with LowerTriangularArrays
    @testset "LowerTriangularArrays kernel test" begin

        @kernel function test_lta_kernel!(A, B, C)
            I = @index(Global, Linear)
            A[I] = B[I] * C[I]
        end

        # Test parameters
        L = 5   # degree
        M = 4   # order
        N = 2   # number of matrices
        NF = Float32

        arch = SpeedyWeather.CPU()

        # Create test arrays
        B = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}}, L, M, N))
        C = on_architecture(arch, rand(LowerTriangularArray{Complex{NF}}, L, M, N))
        A = similar(B)

        expected = B .* C

        # Run the kernel
        launch!(arch, :lmk, size(A), test_lta_kernel!, A, B, C)
        synchronize(arch)

        # Verify results
        @test A ≈ expected
    end

    @testset "RingGrids kernel test" begin

        Grid = HEALPixGrid
        @kernel function test_ringgrids_kernel!(A, B, C)
            I = @index(Global, Linear)
            A[I] = B[I] * C[I]
        end

        # Test parameters
        nlat_half = 5   # degree
        nlayers = 5
        NF = Float32

        arch = SpeedyWeather.CPU()

        # Create test arrays
        B = on_architecture(arch, rand(Grid{NF}, nlat_half, nlayers))
        C = on_architecture(arch, rand(Grid{NF}, nlat_half, nlayers))
        A = similar(B)

        expected = B .* C

        # Run the kernel
        launch!(arch, :ijk, size(A), test_ringgrids_kernel!, A, B, C)
        synchronize(arch)

        # Verify results
        @test A ≈ expected

    end     

end 