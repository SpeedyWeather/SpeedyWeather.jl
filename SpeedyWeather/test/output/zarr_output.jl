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

    # time axis length matches Number of expected outputs (IC + period/interval)
    expected_times = Int(period / output.interval) + 1
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

    # MSLP is a 2D variable derived on-the-fly from the 3D temperature/log(p_s)
    # in a PrimitiveWet default output set. Regression test: an earlier
    # ZarrOutput tried to interpolate the 3D temperature into a 2D output and
    # threw a DimensionMismatch. Skip the t=0 snapshot — physics hasn't run yet
    # so the surface-air temperature buffer is still uninitialized there.
    @test haskey(g.arrays, "mslp")
    nx, ny = RingGrids.matrix_size(output.field2D)
    @test size(g["mslp"]) == (nx, ny, Int(period / output.interval) + 1)
    @test all(isfinite, g["mslp"][:, :, 2:end])
end

@testset "ZarrOutput for variables with custom output! (2D-from-3D)" begin
    # tsurf, u10/v10 each have specialized output! methods that derive a 2D
    # field from 3D model state on the fly. Make sure they all interpolate
    # cleanly through the ZarrOutput backend.
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_zarrtests_2dfrom3d_")
    period = Day(1)

    spectral_grid = SpectralGrid(trunc = 5, nlayers = 4)
    output = ZarrOutput(
        spectral_grid, PrimitiveWet;
        path = tmp_output_path, write_restart = false,
    )
    add!(
        output,
        SpeedyWeather.SurfaceTemperatureOutput(),
        SpeedyWeather.ZonalVelocity10mOutput(),
        SpeedyWeather.MeridionalVelocity10mOutput(),
    )
    model = PrimitiveWetModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true; period)
    @test simulation.model.feedback.nans_detected == false

    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))
    nx, ny = RingGrids.matrix_size(output.field2D)
    nt = Int(period / output.interval) + 1

    for name in ("tsurf", "u10", "v10")
        @test haskey(g.arrays, name)
        @test size(g[name]) == (nx, ny, nt)
        @test all(isfinite, g[name][:, :, 2:end])
    end
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
    @test size(z_vor) == (nlon, nlat, nlayer, Int(period / output.interval) + 1)
    @test z_vor.attrs["_ARRAY_DIMENSIONS"] == ["time", "layer", "lat", "lon"]

    # 3D variable (no layer dim): shape (time, lat, lon)
    z_eta = g["eta"]
    @test size(z_eta) == (nlon, nlat, Int(period / output.interval) + 1)
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

@testset "ZarrOutput spatial chunking" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_zarrtests_spatial_")
    period = Day(1)

    spectral_grid = SpectralGrid(trunc = 5, nlayers = 4)
    output = ZarrOutput(
        spectral_grid, PrimitiveDry;
        path = tmp_output_path, write_restart = false,
        lon_chunk = 4, lat_chunk = 4, vertical_chunk = 2,
    )
    model = PrimitiveDryModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true; period)

    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))

    # 4D variable: chunks come back in Julia (column-major) order: (lon, lat, layer, time)
    z_temp = g["temp"]
    nlon, nlat = RingGrids.matrix_size(output.field2D)
    @test z_temp.metadata.chunks[1] == 4   # lon_chunk
    @test z_temp.metadata.chunks[2] == 4   # lat_chunk
    @test z_temp.metadata.chunks[3] == 2   # vertical_chunk
    @test z_temp.metadata.chunks[4] == max(output.time_chunk, 1)

    # 3D variable (no layer dim): (lon, lat, time)
    z_mslp = g["mslp"]
    @test z_mslp.metadata.chunks[1] == 4
    @test z_mslp.metadata.chunks[2] == 4
    @test z_mslp.metadata.chunks[3] == max(output.time_chunk, 1)

    # values still finite (skip t=0; physics buffers not populated at IC for mslp)
    @test all(isfinite, g["temp"][:, :, :, :])
end

@testset "ZarrOutput chunk clamping to dim size" begin
    # Requesting chunk sizes larger than the array must be clamped down rather
    # than rejected by Zarr (chunk ≤ shape is required by the format).
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_zarrtests_clamp_")
    period = Day(1)

    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
    output = ZarrOutput(
        spectral_grid, ShallowWater;
        path = tmp_output_path, write_restart = false,
        lon_chunk = 100_000, lat_chunk = 100_000, vertical_chunk = 100_000,
    )
    model = ShallowWaterModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true; period)

    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))
    z_vor = g["vor"]
    nlon, nlat = RingGrids.matrix_size(output.field2D)
    @test z_vor.metadata.chunks[1] == nlon
    @test z_vor.metadata.chunks[2] == nlat
    @test z_vor.metadata.chunks[3] == 1     # nlayers=1 ⇒ z clamped to 1
end

@testset "ZarrOutput ensemble members into one store" begin
    # Emulate parallel ensemble members within a single process by running members
    # 1..N sequentially: member 1 (creator) builds the shared store + readiness marker,
    # members 2..N find the marker, open the existing store and write their own slice.
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_zarrtests_ensemble_")
    period = Day(1)
    ensemble_size = 3

    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
    initial_conditions = ZonalJet(spectral_grid)    # deterministic IC, identical for all members
    store_path = ""
    for member in 1:ensemble_size
        output = ZarrOutput(
            spectral_grid, ShallowWater;
            path = tmp_output_path,
            ensemble_index = member, ensemble_size = ensemble_size,
        )
        model = ShallowWaterModel(spectral_grid; output, initial_conditions)
        simulation = initialize!(model)
        run!(simulation, output = true; period)
        @test simulation.model.feedback.nans_detected == false
        # all members must resolve to the same shared run folder / store
        member == 1 && (store_path = joinpath(model.output.run_path, model.output.filename))
        @test joinpath(model.output.run_path, model.output.filename) == store_path
    end

    # side files: members run as parallel processes but share one run folder, so the
    # creator keeps the default filenames and writer members get a _member suffix
    run_path = dirname(store_path)
    for filename in ("parameters", "progress", "restart")
        extension = filename == "restart" ? ".jld2" : ".txt"
        @test isfile(joinpath(run_path, filename * extension))              # creator (member 1)
        for member in 2:ensemble_size                                       # writer members
            @test isfile(joinpath(run_path, filename * "_member$member" * extension))
        end
    end

    # a writer member configured inconsistently with the creator's store errors early
    output = ZarrOutput(
        spectral_grid, ShallowWater;
        path = tmp_output_path,
        ensemble_index = 2, ensemble_size = ensemble_size + 1,  # ensemble size mismatch
    )
    model = ShallowWaterModel(spectral_grid; output, initial_conditions)
    simulation = initialize!(model)
    @test_throws ErrorException run!(simulation, output = true; period)

    # metadata was consolidated by the creator for faster opening (xarray etc.)
    @test isfile(joinpath(store_path, ".zmetadata"))

    g = Zarr.zopen(store_path)

    # ensemble coordinate exists and has length ensemble_size
    @test haskey(g.arrays, "ensemble")
    @test g["ensemble"][:] == collect(1:ensemble_size)

    # a single shared time axis was written (only the creator writes it)
    nlon = length(g["lon"][:])
    nlat = length(g["lat"][:])
    nt = length(g["time"][:])

    # 3D variable gains a trailing ensemble axis; ensemble is first in row-major dims
    z_vor = g["vor"]
    @test ndims(z_vor) == 5
    @test size(z_vor) == (nlon, nlat, spectral_grid.nlayers, nt, ensemble_size)
    @test z_vor.attrs["_ARRAY_DIMENSIONS"] == ["ensemble", "time", "layer", "lat", "lon"]
    # ensemble axis is chunked with size 1 so members write disjoint chunk files
    @test z_vor.metadata.chunks[end] == 1

    # every member wrote finite data into its own ensemble slice; since all members share
    # the same deterministic IC and there's no stochastic physics in ShallowWater, every
    # ensemble slice should be (approximately) identical
    vor_1 = g["vor"][:, :, :, :, 1]
    @test all(isfinite, vor_1)
    for e in 2:ensemble_size
        vor_e = g["vor"][:, :, :, :, e]
        @test all(isfinite, vor_e)
        @test vor_e ≈ vor_1
    end
end

@testset "ZarrOutput second run! creates a new store" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_zarrtests_rerun_")
    period = Day(1)

    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
    output = ZarrOutput(
        spectral_grid, ShallowWater;
        path = tmp_output_path, write_restart = false, id = "rerun",
    )
    model = ShallowWaterModel(spectral_grid; output)
    simulation = initialize!(model)

    run!(simulation, output = true; period)
    first_path = model.output.run_path
    n_first = length(Zarr.zopen(joinpath(first_path, output.filename))["time"][:])

    run!(simulation, output = true; period)
    second_path = model.output.run_path

    # paths differ (new run folder) and the old store wasn't grown by the
    # second run
    @test first_path != second_path
    @test isdir(joinpath(first_path, output.filename))
    @test isdir(joinpath(second_path, output.filename))
    @test length(Zarr.zopen(joinpath(first_path, output.filename))["time"][:]) == n_first
end
