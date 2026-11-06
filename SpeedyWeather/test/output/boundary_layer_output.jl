using NCDatasets

@testset "Surface roughness output writers" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
    period = Day(1)

    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    output = NetCDFOutput(spectral_grid, path = tmp_output_path)
    model = PrimitiveWetModel(spectral_grid; output)
    model.feedback.verbose = false

    # the surface roughness and boundary layer drag output groups under test
    add!(
        model,
        OutputVariables.SurfaceRoughness(),
        OutputVariables.BoundaryLayerDrag(),
    )

    simulation = initialize!(model)
    run!(simulation, period = Day(1), output = true)
    @test simulation.model.feedback.nans_detected == false

    ds = NCDataset(joinpath(model.output.run_path, model.output.filename))
    for key in ("sr_h", "sr_m", "bld_h", "bld_m", "bld_q")
        @test haskey(ds, key)

        # 2D + time variable (no vertical dimension)
        nx, ny, nt = size(ds[key])
        @test (nx, ny) == RingGrids.matrix_size(output.field2D)
        @test nt == Int(period / output.interval) + 1

        # values must be finite, i.e. the `path` actually resolved (a broken path would
        # be caught and written as the missing_value NaN by the generic output writer)
        @test all(x -> ismissing(x) || isfinite(x), ds[key][:, :, end])
    end
    close(ds)
end
