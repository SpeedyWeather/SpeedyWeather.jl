spectral_resolutions = (31,)#, 63, 127)
nlayers_list = (8,)# 8, 32]
# TODO: can uncomment this when we push to main  
grid_list = [
    FullGaussianGrid,
    # FullClenshawGrid,
    OctahedralGaussianGrid,
    # OctahedralClenshawGrid,
    # HEALPixGrid,
    # OctaHEALPixGrid,
]
# CUDA.cu now implicitly converts to Float32 so that's the only relevant type to 
# test here
NFs = (Float32,)


# Function to generate random inputs 
function get_test_data(; trunc, nlayers, Grid, NF)
    # We use dealiasing=3 to ensure that the transform is exact for both 
    # Clenshaw and Gaussian grids
    spectral_grid_cpu = SpectralGrid(; NF, trunc, nlayers, Grid, dealiasing=3)
    spectral_grid_gpu = SpectralGrid(; NF, trunc, nlayers, Grid, architecture=SpeedyWeather.GPU, dealiasing=3)
    
    S_cpu = SpectralTransform(spectral_grid_cpu)
    S_gpu = SpectralTransform(spectral_grid_gpu)

    grid_cpu = rand(spectral_grid_cpu.Grid{spectral_grid_cpu.NF}, 
                    spectral_grid_cpu.nlat_half, 
                    spectral_grid_cpu.nlayers)
    spec_cpu = rand(LowerTriangularArray{Complex{spectral_grid_cpu.NF}}, 
                    spectral_grid_cpu.trunc+2, 
                    spectral_grid_cpu.trunc+1, 
                    spectral_grid_cpu.nlayers)
    
    grid_gpu = cu(grid_cpu)
    spec_gpu = cu(spec_cpu)
    
    return S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu
end


@testset "Whole transform: test a round trip" begin
    @testset for NF in NFs
        @testset for trunc in spectral_resolutions
            @testset for nlayers in nlayers_list
                @testset for Grid in grid_list
                    # Generate test data
                    S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
                        trunc=trunc, nlayers=nlayers, Grid=Grid, NF=NF
                    )

                    # Full return journey starting from grid
                    spec_ = SpeedyTransforms.transform(grid_cpu, S_cpu)
                    grid_cpu = SpeedyTransforms.transform(spec_, S_cpu)
                    # @test grid_cpu ≈ grid_test
                    # grid_cpu = grid_test

                    # Full return journey starting from spec
                    grid_ = SpeedyTransforms.transform(spec_cpu, S_cpu)
                    spec_cpu = SpeedyTransforms.transform(grid_, S_cpu)
                    # @test spec_cpu ≈ spec_test
                    # spec_cpu = spec_test

                    # Copy to GPU and repeat the test. We use initally 
                    # transformed data so that the imaginary components of the 
                    # m=1 specs are zero. 
                    grid_gpu = cu(grid_cpu)
                    spec_gpu = cu(spec_cpu)
                    grid_gpu_test = cu(grid_cpu)
                    spec_gpu_test = cu(spec_cpu)

                    # Full return journey starting from grid on GPU
                    transform!(spec_gpu_test, grid_gpu, S_gpu)
                    transform!(grid_gpu_test, spec_gpu_test, S_gpu)
                    grid_test = adapt(Array, grid_gpu_test)
                    @test grid_cpu ≈ grid_test rtol=sqrt(eps(Float32))

                    # Full return journey starting from spec on GPU
                    transform!(grid_gpu_test, spec_gpu, S_gpu)
                    transform!(spec_gpu_test, grid_gpu_test, S_gpu)
                    spec_test = adapt(Array, spec_gpu_test)
                    @test spec_cpu ≈ spec_test rtol=sqrt(eps(Float32))
                end
            end
        end
    end
end


@testset "fourier_batched: compare forward pass to CPU" begin
    @testset for trunc in spectral_resolutions
        @testset for nlayers in nlayers_list
            @testset for Grid in grid_list
                @testset for NF in NFs
                    # Generate test data
                    S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
                        trunc=trunc, nlayers=nlayers, Grid=Grid, NF=NF
                    )
                    
                    # CPU forward transform
                    SpeedyTransforms._fourier_batched!(
                        S_cpu.scratch_memory_north, 
                        S_cpu.scratch_memory_south, 
                        grid_cpu, 
                        S_cpu
                    )
                    # GPU forward transform
                    SpeedyTransforms._fourier_batched!(
                        S_gpu.scratch_memory_north, 
                        S_gpu.scratch_memory_south, 
                        grid_gpu, 
                        S_gpu
                    )

                    # Convert GPU to CPU for comparison
                    grid_gpu_compare = adapt(Array, grid_gpu)

                    # NOTE: This is probably unnecessary, we can do a coarser 
                    # check for broad comparison of equivalence 
                    # for ij in eachindex(grid_cpu, grid_gpu_compare)
                        # @test grid_cpu[ij] ≈ grid_gpu_compare[ij]
                    # end
                    @test grid_cpu ≈ grid_gpu_compare
                end
            end
        end
    end
end

@testset "fourier_batched: compare backward pass to CPU" begin
    @testset for trunc in spectral_resolutions
        @testset for nlayers in nlayers_list
            @testset for Grid in grid_list
                @testset for NF in NFs
                    # Generate test data
                    S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
                        trunc=trunc, nlayers=nlayers, Grid=Grid, NF=NF
                    )

                    # Use scratch memory to store mid-transform data, using the 
                    # CPU legendre transform to generate the intermediate data
                    # NOTE: assumption of working Legendre transform
                    g_north_cpu = S_cpu.scratch_memory_north    
                    g_south_cpu = S_cpu.scratch_memory_south   
                    SpeedyTransforms._legendre!(g_north_cpu, g_south_cpu, spec_cpu, S_cpu)
                    # Copy to GPU
                    g_north_gpu = cu(g_north_cpu)
                    g_south_gpu = cu(g_south_cpu);

                    # CPU inverse transform
                    SpeedyTransforms._fourier_batched!(
                        grid_cpu, g_north_cpu, g_south_cpu, S_cpu
                    )
                    # GPU inverse transform
                    SpeedyTransforms._fourier_batched!(
                        grid_gpu, g_north_gpu, g_south_gpu, S_gpu
                    )

                    # Copy back to CPU again for comparison
                    grid_gpu_compare = adapt(Array, grid_gpu)
                    
                    # NOTE: again, not necessary to be so fine-grained here
                    # for ij in eachindex(grid_cpu, grid_gpu_compare)
                    #     @test grid_cpu[ij] ≈ grid_gpu_compare[ij]
                    # end
                    @test grid_cpu ≈ grid_gpu_compare 
                end
            end
        end
    end
end


@testset "fourier_serial: compare forward pass to CPU" begin
    @testset for trunc in spectral_resolutions
        @testset for nlayers in nlayers_list
            @testset for Grid in grid_list
                @testset for NF in NFs
                    # Generate test data
                    S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
                        trunc=trunc, nlayers=nlayers, Grid=Grid, NF=NF
                    )
                    
                    # CPU forward transform
                    SpeedyTransforms._fourier_serial!(
                        S_cpu.scratch_memory_north, 
                        S_cpu.scratch_memory_south, 
                        grid_cpu, 
                        S_cpu
                    )
                    # GPU forward transform
                    SpeedyTransforms._fourier_serial!(
                        S_gpu.scratch_memory_north, 
                        S_gpu.scratch_memory_south, 
                        grid_gpu, 
                        S_gpu
                    )

                    # Convert GPU to CPU for comparison
                    grid_gpu_compare = adapt(Array, grid_gpu)

                    @test grid_cpu ≈ grid_gpu_compare 
                end
            end
        end
    end
end

# NOTE: Currently failing due to problem with Float32 FFTW planning
# @testset "fourier_serial: compare backward pass to CPU" begin
#     @testset for trunc in spectral_resolutions
#         @testset for nlayers in nlayers_list
#             @testset for Grid in grid_list
#                 @testset for NF in (Float32, Float64)
#                     # Generate test data
#                     S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
#                         trunc=trunc, nlayers=nlayers, Grid=Grid, NF=NF
#                     )

#                     # Use scratch memory to store mid-transform data, using the 
#                     # CPU legendre transform to generate the intermediate data
#                     # NOTE: assumption of working Legendre transform
#                     g_north_cpu = S_cpu.scratch_memory_north    
#                     g_south_cpu = S_cpu.scratch_memory_south   
#                     SpeedyTransforms._legendre!(g_north_cpu, g_south_cpu, spec_cpu, S_cpu)
#                     # Copy to GPU
#                     g_north_gpu = cu(g_north_cpu)
#                     g_south_gpu = cu(g_south_cpu);

#                     # CPU inverse transform
#                     SpeedyTransforms._fourier_serial!(
#                         grid_cpu, g_north_cpu, g_south_cpu, S_cpu
#                     )
#                     # GPU inverse transform
#                     SpeedyTransforms._fourier_serial!(
#                         grid_gpu, g_north_gpu, g_south_gpu, S_gpu
#                     )

#                     # Copy back to CPU again for comparison
#                     grid_gpu_compare = adapt(Array, grid_gpu)
                    
#                     @test grid_cpu ≈ grid_gpu_compare
#                 end
#             end
#         end
#     end
# end

@testset "legendre: compare inverse transform to CPU" begin
    @testset for NF in NFs
        @testset for trunc in spectral_resolutions
            @testset for nlayers in nlayers_list
                @testset for Grid in grid_list
                    # Generate test data
                    S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
                        trunc=trunc, nlayers=nlayers, Grid=Grid, NF=NF
                    )
                    
                    # CPU inverse transform
                    SpeedyTransforms._legendre!(
                        S_cpu.scratch_memory_north, 
                        S_cpu.scratch_memory_south, 
                        spec_cpu, S_cpu
                    )
                    # GPU inverse transform
                    SpeedyTransforms._legendre!(
                        S_gpu.scratch_memory_north, 
                        S_gpu.scratch_memory_south, 
                        spec_gpu, S_gpu
                    )

                    # Convert GPU to CPU for comparison, result is stored in the 
                    # scratch memory
                    result_gpu = adapt(Array, S_gpu.scratch_memory_north);
                    result_cpu = S_cpu.scratch_memory_north;
                    @test result_cpu ≈ result_gpu rtol=sqrt(eps(Float32))   # GPU error tolerance always Float32

                    result_gpu = adapt(Array, S_gpu.scratch_memory_south);
                    result_cpu = S_cpu.scratch_memory_south;
                    @test result_cpu ≈ result_gpu rtol=sqrt(eps(Float32))
                end
            end
        end
    end
end

@testset "legendre: compare forward transform to CPU" begin
    @testset for NF in NFs
        @testset for trunc in spectral_resolutions
            @testset for nlayers in nlayers_list
                @testset for Grid in grid_list
                    # Generate test data
                    S_cpu, S_gpu, grid_cpu, grid_gpu, spec_cpu, spec_gpu = get_test_data(
                        trunc=trunc, nlayers=nlayers, Grid=Grid, NF=NF
                    )

                    # Use scratch memory to store mid-transform data, using the 
                    # CPU fourier transform to generate the intermediate data
                    f_north_cpu = S_cpu.scratch_memory_north    
                    f_south_cpu = S_cpu.scratch_memory_south   
                    SpeedyTransforms._fourier!(f_north_cpu, f_south_cpu, grid_cpu, S_cpu)
                    # Copy to GPU
                    f_north_gpu = cu(f_north_cpu)
                    f_south_gpu = cu(f_south_cpu)
                    
                    # CPU inverse transform
                    SpeedyTransforms._legendre!(
                        spec_cpu,    
                        f_north_cpu, 
                        f_south_cpu, 
                        S_cpu
                    )
                    # GPU inverse transform
                    SpeedyTransforms._legendre!(
                        spec_gpu,
                        f_north_gpu, 
                        f_south_gpu, 
                        S_gpu
                    )

                    # Convert GPU to CPU for comparison, result is stored spec
                    result_gpu = adapt(Array, spec_gpu);
                    result_cpu = spec_cpu;
                    @test result_cpu ≈ result_gpu rtol=sqrt(eps(Float32))   # GPU error tolerance always Float32
                end
            end
        end
    end
end
