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

    field_cpu  = rand(NF, spectral_grid_cpu.grid,
                        spectral_grid_cpu.nlayers)
    coeffs_cpu = rand(LowerTriangularArray{Complex{NF}}, 
                        spectral_grid_cpu.spectrum, 
                        spectral_grid_cpu.nlayers)
    
    grid_gpu = on_architecture(S_gpu.architecture, field_cpu)
    spec_gpu = on_architecture(S_gpu.architecture, coeffs_cpu)
    
    return S_cpu, S_gpu, field_cpu, grid_gpu, coeffs_cpu, spec_gpu
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

                    gpu_arch = S_gpu.architecture
                    cpu_arch = S_cpu.architecture
                    
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
                    grid_gpu = on_architecture(gpu_arch, grid_cpu)
                    spec_gpu = on_architecture(gpu_arch, spec_cpu)
                    grid_gpu_test = on_architecture(gpu_arch, grid_cpu)
                    spec_gpu_test = on_architecture(gpu_arch, spec_cpu)

                    # Full return journey starting from grid on GPU
                    transform!(spec_gpu_test, grid_gpu, S_gpu)
                    transform!(grid_gpu_test, spec_gpu_test, S_gpu)
                    grid_test = on_architecture(cpu_arch, grid_gpu_test)
                    @test grid_cpu ≈ grid_test rtol=sqrt(eps(Float32))

                    # Full return journey starting from spec on GPU
                    transform!(grid_gpu_test, spec_gpu, S_gpu)
                    transform!(spec_gpu_test, grid_gpu_test, S_gpu)
                    spec_test = on_architecture(cpu_arch, spec_gpu_test)
                    @test spec_cpu ≈ spec_test rtol=sqrt(eps(Float32))
                end
            end
        end
    end
end

@testset "Whole transform: works in 2D" begin
    S_cpu, S_gpu, field_cpu, field_gpu, spec_cpu, spec_gpu = get_test_data(
                        trunc=spectral_resolutions[1], nlayers=nlayers_list[1], Grid=grid_list[1], NF=NFs[1]
                    )

    cpu_arch = S_cpu.architecture
                    
    field_2d_cpu = field_cpu[:,1]
    spec_2d_cpu = spec_cpu[:,1]

    spec_2d_cpu_res = transform(field_2d_cpu, S_cpu)
    field_2d_cpu_res = transform(spec_2d_cpu, S_cpu)

    field_2d_gpu = field_gpu[:,1]
    spec_2d_gpu = spec_gpu[:,1]

    spec_2d_gpu_res = transform(field_2d_gpu, S_gpu)
    field_2d_gpu_res = transform(spec_2d_gpu, S_gpu)

    @test field_2d_cpu_res ≈ on_architecture(cpu_arch, field_2d_gpu_res)
    @test spec_2d_cpu_res ≈ on_architecture(cpu_arch, spec_2d_gpu_res)
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

                    cpu_arch = S_cpu.architecture
                    gpu_arch = S_gpu.architecture
                    
                    # CPU forward transform
                    SpeedyTransforms._fourier_batched!(
                        S_cpu.scratch_memory.north, 
                        S_cpu.scratch_memory.south, 
                        grid_cpu, 
                        S_cpu.scratch_memory.column,
                        S_cpu
                    )
                    # GPU forward transform
                    SpeedyTransforms._fourier_batched!(
                        S_gpu.scratch_memory.north, 
                        S_gpu.scratch_memory.south, 
                        grid_gpu, 
                        S_gpu.scratch_memory.column,
                        S_gpu
                    )

                    # Convert GPU to CPU for comparison
                    grid_gpu_compare = on_architecture(cpu_arch, grid_gpu)

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

                    cpu_arch = S_cpu.architecture
                    gpu_arch = S_gpu.architecture

                    # Use scratch memory to store mid-transform data, using the 
                    # CPU legendre transform to generate the intermediate data
                    # NOTE: assumption of working Legendre transform
                    g_north_cpu = S_cpu.scratch_memory.north    
                    g_south_cpu = S_cpu.scratch_memory.south   
                    SpeedyTransforms._legendre!(g_north_cpu, g_south_cpu, spec_cpu, S_cpu)
                    # Copy to GPU
                    g_north_gpu = cu(g_north_cpu)
                    g_south_gpu = cu(g_south_cpu);

                    # CPU inverse transform
                    SpeedyTransforms._fourier_batched!(
                        grid_cpu, g_north_cpu, g_south_cpu, S_cpu.scratch_memory.column, S_cpu
                    )
                    # GPU inverse transform
                    SpeedyTransforms._fourier_batched!(
                        grid_gpu, g_north_gpu, g_south_gpu, S_gpu.scratch_memory.column, S_gpu
                    )

                    # Copy back to CPU again for comparison
                    grid_gpu_compare = on_architecture(cpu_arch, grid_gpu)
                    
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

                    gpu_arch = S_gpu.architecture
                    cpu_arch = S_cpu.architecture
                    
                    # CPU forward transform
                    SpeedyTransforms._fourier_serial!(
                        S_cpu.scratch_memory.north, 
                        S_cpu.scratch_memory.south, 
                        grid_cpu, 
                        S_cpu.scratch_memory.column,
                        S_cpu
                    )
                    # GPU forward transform
                    SpeedyTransforms._fourier_serial!(
                        S_gpu.scratch_memory.north, 
                        S_gpu.scratch_memory.south, 
                        grid_gpu, 
                        S_gpu.scratch_memory.column,
                        S_gpu
                    )

                    # Convert GPU to CPU for comparison
                    grid_gpu_compare = on_architecture(cpu_arch, grid_gpu)

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

                    gpu_arch = S_gpu.architecture
                    cpu_arch = S_cpu.architecture   
                    
                    # CPU inverse transform
                    SpeedyTransforms._legendre!(
                        S_cpu.scratch_memory.north, 
                        S_cpu.scratch_memory.south, 
                        spec_cpu, S_cpu.scratch_memory.column, S_cpu
                    )
                    # GPU inverse transform
                    SpeedyTransforms._legendre!(
                        S_gpu.scratch_memory.north, 
                        S_gpu.scratch_memory.south, 
                        spec_gpu, S_gpu.scratch_memory.column, S_gpu
                    )

                    # Convert GPU to CPU for comparison, result is stored in the 
                    # scratch memory
                    result_gpu = on_architecture(cpu_arch, S_gpu.scratch_memory.north);
                    result_cpu = S_cpu.scratch_memory.north;
                    @test result_cpu ≈ result_gpu rtol=sqrt(eps(Float32))   # GPU error tolerance always Float32

                    result_gpu = on_architecture(cpu_arch, S_gpu.scratch_memory.south);
                    result_cpu = S_cpu.scratch_memory.south;
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

                    gpu_arch = S_gpu.architecture
                    cpu_arch = S_cpu.architecture

                    # Use scratch memory to store mid-transform data, using the 
                    # CPU fourier transform to generate the intermediate data
                    f_north_cpu = S_cpu.scratch_memory.north    
                    f_south_cpu = S_cpu.scratch_memory.south   
                    SpeedyTransforms._fourier!(f_north_cpu, f_south_cpu, grid_cpu, S_cpu)
                    # Copy to GPU
                    f_north_gpu = cu(f_north_cpu)
                    f_south_gpu = cu(f_south_cpu)
                    
                    # CPU inverse transform
                    SpeedyTransforms._legendre!(
                        spec_cpu,    
                        f_north_cpu, 
                        f_south_cpu,
                        S_cpu.scratch_memory.column,
                        S_cpu
                    )
                    # GPU inverse transform
                    SpeedyTransforms._legendre!(
                        spec_gpu,
                        f_north_gpu, 
                        f_south_gpu, 
                        S_gpu.scratch_memory.column,
                        S_gpu
                    )

                    # Convert GPU to CPU for comparison, result is stored spec
                    result_gpu = on_architecture(cpu_arch, spec_gpu);
                    result_cpu = spec_cpu;
                    @test result_cpu ≈ result_gpu rtol=sqrt(eps(Float32))   # GPU error tolerance always Float32
                end
            end
        end
    end
end
