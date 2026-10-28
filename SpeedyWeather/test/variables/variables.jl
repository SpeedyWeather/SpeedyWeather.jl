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

@testset "Output variables group" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 2)
    model = PrimitiveWetModel(spectral_grid)
    vars = Variables(model)

    # accumulated precipitation and runoff live in the output group, not the parameterizations
    for name in (:rain_convection, :snow_convection, :rain_large_scale, :snow_large_scale)
        @test haskey(vars.output, name)
        @test !haskey(vars.parameterizations, name)
    end
    @test haskey(vars.output.land, :river_runoff)
    @test !haskey(vars.parameterizations.land, :river_runoff)

    # OutputVariable sorts into the output group with namespaces
    v = OutputVariable(:my_accumulation, SpeedyWeather.Grid2D(), namespace = :land)
    @test SpeedyWeather.nonparametric_type(v) == OutputVariable
end
