using CUDA
using SpeedyWeather
using Adapt
using Test

spectral_resolutions = (31, 63)
nlayers_list = [8]
grid_list = [
    FullGaussianGrid,
]
function get_test_data(; trunc, nlayers, Grid, NF)
    spectral_grid_cpu = SpectralGrid(; NF, trunc, nlayers, Grid)
    spectral_grid_gpu = SpectralGrid(; NF, trunc, nlayers, Grid, device=SpeedyWeather.GPU())
    
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

@testset "legendre: compare forward transform to CPU" begin
    @testset for NF in (Float32,)
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