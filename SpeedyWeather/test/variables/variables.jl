@testset "Variables on_architecture CPU → JLArray GPU → CPU" begin
    using JLArrays

    NF = Float32
    spectral_grid = SpectralGrid(; NF, nlayers = 2, Grid = FullGaussianGrid)
    model = PrimitiveWetModel(spectral_grid)
    vars_cpu = Variables(model)
    vars_cpu.prognostic.vorticity .= zeros(Complex{NF}, spectral_grid.spectrum, 2, 2)
    jl_arch = SpeedyWeather.GPU(JLArrays.JLBackend())
    vars_gpu = on_architecture(jl_arch, vars_cpu)

    # vorticity should now live on JLArray
    @test parent(vars_gpu.prognostic.vorticity) isa JLArray

    # round-trip back to CPU should reproduce original values
    vars_back = on_architecture(SpeedyWeather.CPU(), vars_gpu)
    @test Array(parent(vars_back.prognostic.vorticity)) == Array(parent(vars_cpu.prognostic.vorticity))
    @test Array(parent(vars_back.prognostic.temperature)) == Array(parent(vars_cpu.prognostic.temperature))
end
