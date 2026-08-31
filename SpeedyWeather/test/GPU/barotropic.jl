using JLD2
using Zarr

@testset "GPU Barotropic (with MatrixSpectralTransform)" begin
    spectral_grid = SpectralGrid(truncation = 33, nlayers = 1, architecture = SpeedyWeather.GPU())
    spectral_transform = MatrixSpectralTransform(spectral_grid)
    model = BarotropicModel(spectral_grid; spectral_transform)
    simulation = initialize!(model)
    run!(simulation, steps = 4)

    @test simulation.model.feedback.nans_detected == false
end

@testset "GPU Barotropic with JLD2Output" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_gpu_jld2_")
    spectral_grid = SpectralGrid(truncation = 33, nlayers = 1, architecture = SpeedyWeather.GPU())
    output = JLD2Output(path = tmp_output_path, id = "gpu-jld2", write_restart = false)
    model = BarotropicModel(spectral_grid; output)
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
    spectral_grid = SpectralGrid(truncation = 33, nlayers = 1, architecture = SpeedyWeather.GPU())
    output = ArrayOutput()
    model = BarotropicModel(spectral_grid; output)
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

@testset "GPU Barotropic with HEALPixOutput, interpolation skipped" begin
    # model already on a HEALPixGrid *on GPU*, output on the same grid on CPU. `grids_match`
    # compares grid type and nlat_half only, so no interpolator is built and the output path
    # takes the plain `copyto!` branch — after the GPU -> CPU transfer of the model field.
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_gpu_healpix_skip_")

    spectral_grid = SpectralGrid(
        truncation = 33, nlayers = 1, Grid = HEALPixGrid,
        architecture = SpeedyWeather.GPU(),
    )
    output = HEALPixOutput(
        spectral_grid, Barotropic;
        path = tmp_output_path, id = "gpu-healpix-skip", write_restart = false,
    )
    @test isnothing(output.interpolator)                # same grid in and out
    @test output.field2D.grid.nlat_half == spectral_grid.grid.nlat_half

    model = BarotropicModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, steps = 3, output = true)

    @test simulation.model.feedback.nans_detected == false

    g = Zarr.zopen(joinpath(output.run_path, output.filename))
    npix = RingGrids.get_npoints(spectral_grid.grid)
    @test size(g["vor"], 1) == npix
    @test all(isfinite, g["vor"][:, :, :])

    # coordinates are the model grid's own, computed on CPU
    londs, latds = RingGrids.get_londlatds(on_architecture(SpeedyWeather.CPU(), spectral_grid.grid))
    @test g["lat"][:] ≈ latds
    @test g["lon"][:] ≈ londs
end
