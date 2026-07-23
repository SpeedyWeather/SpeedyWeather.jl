using Zarr, Dates

@testset "ZarrOutput ensemble across OS processes" begin
    # A genuine multi-processing test for the ensemble ZarrOutput: every member runs as
    # its own OS process (not sequentially in one process like the other ensemble test),
    # so this exercises the real parallel path — the shared run-folder race, the creator's
    # readiness marker, the writer members blocking on it, and the concurrent chunk writes
    # into one shared store. Each member is a fresh `julia` launched with the same project
    # as this test; it needs SpeedyWeather + Zarr, both of which are in the test project.
    ensemble_size = 2

    tmp_output_path = mktempdir(pwd(), prefix = "tmp_zarr_mp_")   # the shared output parent
    work_dir = mktempdir()                                        # worker script + logs

    # worker script: run a single ensemble member. `ensemble_timeout` is kept short so a
    # writer fails fast (rather than hanging for the 600s default) if the creator crashed.
    worker_script = joinpath(work_dir, "run_member.jl")
    write(
        worker_script, """
        using SpeedyWeather, Zarr, Dates
        member        = parse(Int, ARGS[1])
        ensemble_size = parse(Int, ARGS[2])
        path          = ARGS[3]

        spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
        initial_conditions = ZonalJet(spectral_grid)   # deterministic, identical for all members
        output = ZarrOutput(
            spectral_grid, ShallowWater;
            path,
            ensemble_index   = member,
            ensemble_size    = ensemble_size,
            ensemble_timeout = 300,
        )
        model = ShallowWaterModel(spectral_grid; output, initial_conditions)
        simulation = initialize!(model)
        run!(simulation, output = true, period = Day(1))
        exit(simulation.model.feedback.nans_detected ? 1 : 0)
        """
    )

    # inherit this process' project (the test sandbox under Pkg.test, so Zarr is available)
    # and its julia flags (e.g. --check-bounds). Launch ALL members at once with wait=false.
    project = dirname(Base.active_project())
    logs = [joinpath(work_dir, "member$m.log") for m in 1:ensemble_size]
    procs = map(1:ensemble_size) do member
        cmd = `$(Base.julia_cmd()) --project=$project $worker_script $member $ensemble_size $tmp_output_path`
        run(pipeline(cmd; stdout = logs[member], stderr = logs[member]); wait = false)
    end

    # wait for every member and assert a clean exit; surface the log on failure
    for (member, p) in enumerate(procs)
        wait(p)
        success(p) || @error "ensemble member $member failed:\n$(read(logs[member], String))"
        @test success(p)
    end

    # all members must have resolved to exactly one shared run folder
    run_folders = filter(f -> isdir(joinpath(tmp_output_path, f)) && startswith(f, "run"), readdir(tmp_output_path))
    @test length(run_folders) == 1
    run_path = joinpath(tmp_output_path, only(run_folders))
    store_path = joinpath(run_path, "output.zarr")

    g = Zarr.zopen(store_path)

    # the store was built with the ensemble axis and a single shared time axis
    @test haskey(g.arrays, "ensemble")
    @test g["ensemble"][:] == collect(1:ensemble_size)

    # every member wrote finite data into its own ensemble slice (deterministic ShallowWater
    # IC + no stochastic physics ⇒ all slices are approximately identical)
    z_vor = g["vor"]
    @test size(z_vor)[end] == ensemble_size
    vor_1 = g["vor"][:, :, :, :, 1]
    @test all(isfinite, vor_1)
    for e in 2:ensemble_size
        vor_e = g["vor"][:, :, :, :, e]
        @test all(isfinite, vor_e)
        @test vor_e ≈ vor_1
    end

    # each member (creator included) wrote its side files under the _member suffix
    for member in 1:ensemble_size
        @test isfile(joinpath(run_path, "progress_member$member.txt"))
    end
end
