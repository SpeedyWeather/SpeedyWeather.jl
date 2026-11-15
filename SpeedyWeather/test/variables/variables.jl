@testset "Variables on_architecture CPU → JLArray GPU → CPU" begin
    using JLArrays

    NF = Float32
    spectral_grid = SpectralGrid(; NF, nlayers = 2, Grid = FullGaussianGrid)
    model = PrimitiveWetModel(spectral_grid)
    vars_cpu = Variables(model)
    vars_cpu.prognostic.vorticity .= zeros(Complex{NF}, spectral_grid.spectrum, 2, 2)
    jl_arch = SpeedyWeather.GPU(JLArrays.JLBackend())
    vars_gpu = on_architecture(jl_arch, vars_cpu)

    # vorticity should now live on JLArray, as a view of its fused parent
    @test parent(vars_gpu.fused.prognostic).data isa JLArray
    @test parent(vars_gpu.prognostic.vorticity.data) === parent(vars_gpu.fused.prognostic).data

    # round-trip back to CPU should reproduce original values
    vars_back = on_architecture(SpeedyWeather.CPU(), vars_gpu)
    @test Array(parent(vars_back.prognostic.vorticity)) == Array(parent(vars_cpu.prognostic.vorticity))
    @test Array(parent(vars_back.prognostic.temperature)) == Array(parent(vars_cpu.prognostic.temperature))
end

@testset "Implicit solver variables: TensorDim, ValDim, traceable ScalarDim" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 5)
    model = PrimitiveDryModel(spectral_grid)
    vars = Variables(model)
    implicit = vars.dynamics.implicit

    # the implicit solver holds only options, its arrays are variables allocated for the model
    @test implicit.S⁻¹ isa Array{spectral_grid.NF, 3}
    @test size(implicit.S⁻¹) == (spectral_grid.truncation + 1, 5, 5)
    @test size(implicit.R) == (5, 5)
    @test size(implicit.temp_profile) == (5,)
    @test implicit.nlayers === SpeedyWeather.StaticValue{5}()
    @test Val(implicit.nlayers) === Val(5)
    @test implicit.Δt isa Base.RefValue{spectral_grid.NF}
    @test on_architecture(SpeedyWeather.CPU(), implicit.nlayers) === implicit.nlayers

    # the dimensions on their own
    @test SpeedyWeather.ScalarDim(1) == SpeedyWeather.ScalarDim(value = 1, traceable = false)
    @test SpeedyWeather.ScalarDim(traceable = true).traceable

    # no implicit variables for models without the primitive equation solver
    for Model in (BarotropicModel, ShallowWaterModel)
        vars_2D = Variables(Model(SpectralGrid(truncation = 21, nlayers = 1)))
        @test !haskey(vars_2D.dynamics, :implicit)
    end
end
