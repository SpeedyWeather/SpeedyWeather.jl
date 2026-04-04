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
    @test parent(snapshot.prognostic.vor) isa Array
    @test all(isfinite, snapshot.prognostic.vor)
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

    # check IC snapshot
    ic = output.output[1].prognostic
    @test all(isfinite, ic.vor)
    @test all(isfinite, ic.div)
    @test all(isfinite, ic.humid)
    @test all(isfinite, ic.pres)
    @test all(isfinite, ic.temp)

    # check final snapshot
    final_snapshot = output.output[end].prognostic
    @test all(isfinite, final_snapshot.vor)
    @test all(isfinite, final_snapshot.div)
    @test all(isfinite, final_snapshot.humid)
    @test all(isfinite, final_snapshot.pres)
    @test all(isfinite, final_snapshot.temp)
end
