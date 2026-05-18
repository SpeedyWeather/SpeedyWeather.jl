using JLD2

@testset "GPU Barotropic" begin
    spectral_grid = SpectralGrid(trunc = 32, nlayers = 1, architecture = SpeedyWeather.GPU())
    model = BarotropicModel(spectral_grid = spectral_grid)
    simulation = initialize!(model)
    run!(simulation, steps = 4)

    @test simulation.model.feedback.nans_detected == false
end

@testset "GPU Barotropic with JLD2Output" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_gpu_jld2_")
    spectral_grid = SpectralGrid(trunc = 32, nlayers = 1, architecture = SpeedyWeather.GPU())
    output = JLD2Output(path = tmp_output_path, id = "gpu-jld2", write_restart = false)
    model = BarotropicModel(; spectral_grid, output)
    simulation = initialize!(model)
    run!(simulation, steps = 3, output = true)

    @test simulation.model.feedback.nans_detected == false

    # JLD2Output should have transferred data to CPU before writing
    f = jldopen(joinpath(output.run_path, output.filename), "r")
    @test length(f["output_vector"]) > 0

    snapshot = f["output_vector"][1]
    @test parent(snapshot.prognostic.vorticity) isa Array
    @test all(isfinite, snapshot.prognostic.vorticity)
    close(f)
end

@testset "GPU Barotropic with ArrayOutput" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_gpu_jld2_")
    spectral_grid = SpectralGrid(trunc = 32, nlayers = 1, architecture = SpeedyWeather.GPU())
    output = ArrayOutput()
    model = BarotropicModel(; spectral_grid, output)
    simulation = initialize!(model)
    run!(simulation, steps = 3, output = true)

    @test simulation.model.feedback.nans_detected == false

    # check IC snapshot (transfer to CPU first to avoid scalar indexing on GPU)
    ic = output.output[1].prognostic
    @test all(isfinite, on_architecture(SpeedyWeather.CPU(), ic.vorticity))

    # check final snapshot
    final_snapshot = output.output[end].prognostic
    @test all(isfinite, on_architecture(SpeedyWeather.CPU(), final_snapshot.vorticity))
end
