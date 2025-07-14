import SpeedyWeather: array_type

@testset "Architecture / Device handling" begin
    # test basic transfer of data and grids between CPU and GPU
    NF = Float32
    nlayers = 2
    trunc = 32
    nlat_half = 6
    Grid = FullGaussianGrid

    arch_gpu = SpeedyWeather.GPU()
    arch_cpu = SpeedyWeather.CPU()

    grid_cpu = Grid(nlat_half, arch_cpu)

    # transfering Grid 
    grid_gpu = on_architecture(arch_gpu, grid_cpu)
    @test typeof(grid_gpu.architecture) <: typeof(arch_gpu)
    
    # transfering Field
    field_cpu = rand(NF, grid_cpu, nlayers)
    field_gpu = on_architecture(arch_gpu, field_cpu)
    
    @test array_type(field_gpu) <: array_type(arch_gpu) <: CUDA.CuArray
    @test typeof(field_gpu.grid.architecture) <: typeof(arch_gpu) <: SpeedyWeather.GPU

    spectrum_cpu = Spectrum(; trunc, architecture=arch_cpu)

    # transfering Spectrum
    spectrum_gpu = on_architecture(arch_gpu, spectrum_cpu)
    @test typeof(spectrum_gpu.architecture) <: typeof(arch_gpu) <: SpeedyWeather.GPU
    
    spec_cpu = rand(LowerTriangularArray{Complex{NF}}, spectrum_cpu, nlayers)
    spec_gpu = on_architecture(arch_gpu, spec_cpu)
    
    @test array_type(spec_gpu) <: array_type(arch_gpu) <: CUDA.CuArray
    @test typeof(spec_gpu.spectrum.architecture) <: typeof(arch_gpu) <: SpeedyWeather.GPU
end 