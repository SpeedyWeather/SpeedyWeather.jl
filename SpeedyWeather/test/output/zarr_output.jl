using Zarr, Dates

@testset "ZarrOutput type and defaults" begin
    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
    output = ZarrOutput(spectral_grid)
    @test output isa SpeedyWeather.ZarrOutput
    @test output.active == false
    @test output.filename == "output.zarr"
end

@testset "ZarrOutput for ShallowWaterModel" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_zarrtests_sw_")
    period = Day(1)

    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
    output = ZarrOutput(
        spectral_grid, ShallowWater;
        path = tmp_output_path, write_restart = false,
    )
    model = ShallowWaterModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true; period)
    @test simulation.model.feedback.nans_detected == false

    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))

    # all declared output variables made it into the store
    for var in values(output.variables)
        @test haskey(g.arrays, var.name)
    end

    # coordinate arrays were written
    for c in ("lon", "lat", "layer", "time")
        @test haskey(g.arrays, c)
    end

    # time axis length matches Number of expected outputs (IC + period/output_dt)
    expected_times = Int(period / output.output_dt) + 1
    @test length(g["time"][:]) == expected_times
    @test g["time"][1] == 0.0

    # spatial dims and layer count for a 3D variable
    nx, ny, nz, nt = size(g["vor"])
    @test (nx, ny) == RingGrids.matrix_size(output.field2D)
    @test nz == spectral_grid.nlayers
    @test nt == expected_times

    # values are finite (not stuck on the fill value)
    @test all(isfinite, g["vor"][:, :, 1, :])
end

@testset "ZarrOutput for PrimitiveWetModel with soil layers" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_zarrtests_pw_")
    period = Day(1)

    spectral_grid = SpectralGrid(trunc = 5, nlayers = 4)
    output = ZarrOutput(
        spectral_grid, PrimitiveWet;
        path = tmp_output_path, write_restart = false,
    )
    model = PrimitiveWetModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true; period)
    @test simulation.model.feedback.nans_detected == false

    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))
    @test haskey(g.arrays, "temp")
    @test haskey(g.arrays, "humid")
    @test haskey(g.arrays, "soil_layer")    # soil dim coordinate

    # 3D atmosphere variable has correct vertical dim
    @test size(g["temp"], 3) == spectral_grid.nlayers
end

@testset "ZarrOutput add!/delete! variables" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_zarrtests_add_")
    period = Day(1)

    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
    output = ZarrOutput(
        spectral_grid, ShallowWater;
        path = tmp_output_path, write_restart = false,
    )

    # add an extra variable, run, check it shows up
    div_output = SpeedyWeather.DivergenceOutput()
    add!(output, div_output)
    @test haskey(output.variables, Symbol(div_output.name))

    model = ShallowWaterModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true; period)
    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))
    @test haskey(g.arrays, div_output.name)

    # delete!: variable is removed from the dict (Zarr arrays from the previous
    # run remain on disk; we only test the dict-level behaviour here)
    delete!(output, div_output.name)
    @test !haskey(output.variables, Symbol(div_output.name))
end

@testset "ZarrOutput xarray-compatible dimension order" begin
    # _ARRAY_DIMENSIONS must be in row-major (C) order to match the on-disk
    # shape that Zarr stores; otherwise xarray will assign dim names to the
    # wrong axes. Regression test for the original implementation that wrote
    # them in Julia (column-major) order.
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_zarrtests_dims_")
    period = Day(1)

    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
    output = ZarrOutput(
        spectral_grid, ShallowWater;
        path = tmp_output_path, write_restart = false,
    )
    model = ShallowWaterModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true; period)

    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))

    # 4D variable: shape (time, layer, lat, lon) on disk and dims
    # attribute in matching order
    z_vor = g["vor"]
    nlon = length(g["lon"][:])
    nlat = length(g["lat"][:])
    nlayer = length(g["layer"][:])
    @test size(z_vor) == (nlon, nlat, nlayer, Int(period / output.output_dt) + 1)
    @test z_vor.attrs["_ARRAY_DIMENSIONS"] == ["time", "layer", "lat", "lon"]

    # 3D variable (no layer dim): shape (time, lat, lon)
    z_eta = g["eta"]
    @test size(z_eta) == (nlon, nlat, Int(period / output.output_dt) + 1)
    @test z_eta.attrs["_ARRAY_DIMENSIONS"] == ["time", "lat", "lon"]

    # 1D coordinates: dim attribute matches name
    @test g["lon"].attrs["_ARRAY_DIMENSIONS"] == ["lon"]
    @test g["time"].attrs["_ARRAY_DIMENSIONS"] == ["time"]
end

@testset "ZarrOutput compressor and time_chunk options" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_zarrtests_opt_")
    period = Day(1)

    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
    output = ZarrOutput(
        spectral_grid, ShallowWater;
        path = tmp_output_path,
        write_restart = false,
        time_chunk = 3,
        compressor = Zarr.BloscCompressor(clevel = 5),
    )
    model = ShallowWaterModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true; period)

    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))
    z_vor = g["vor"]
    # chunk along the time axis matches `time_chunk`
    @test z_vor.metadata.chunks[end] == 3
    # compressor was applied
    @test z_vor.metadata.compressor isa Zarr.BloscCompressor
end
