@testset "Variables on_architecture CPU → JLArray GPU → CPU" begin
    using JLArrays

    NF = Float32
    spectral_grid = SpectralGrid(; NF, nlayers = 2, Grid = FullGaussianGrid)
    model = PrimitiveWetModel(spectral_grid)
    vars_cpu = Variables(model)
    vars_cpu.prognostic.vor .= zeros(Complex{NF}, spectral_grid.spectrum, 2, 2)
    jl_arch = SpeedyWeather.GPU(JLArrays.JLBackend())
    vars_gpu = on_architecture(jl_arch, vars_cpu)

    # vor should now live on JLArray
    @test parent(vars_gpu.prognostic.vor) isa JLArray

    # round-trip back to CPU should reproduce original values
    vars_back = on_architecture(SpeedyWeather.CPU(), vars_gpu)
    @test Array(parent(vars_back.prognostic.vor)) == Array(parent(vars_cpu.prognostic.vor))
    @test Array(parent(vars_back.prognostic.temp)) == Array(parent(vars_cpu.prognostic.temp))
end
