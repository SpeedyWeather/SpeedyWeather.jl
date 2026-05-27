using NCDatasets, Dates

@testset "Output for BarotropicModel" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
    period = Day(1)

    for Grid in (FullGaussianGrid, FullClenshawGrid, OctahedralGaussianGrid, OctahedralClenshawGrid, HEALPixGrid, OctaHEALPixGrid)
        spectral_grid = SpectralGrid(nlayers = 1)
        model = BarotropicModel(spectral_grid)
        output = NetCDFOutput(model, path = tmp_output_path)
        simulation = initialize!(model; output)
        run!(simulation, output = true; period)
        @test simulation.feedback.nans_detected == false

        # read netcdf file and check that all variables exist
        ds = NCDataset(joinpath(output.run_path, output.filename))
        for key in keys(output.variables)
            @test haskey(ds, key)

            # test dimensions
            nx, ny, nz, nt = size(ds[key])
            @test (nx, ny) == RingGrids.matrix_size(output.field2D)
            @test nz == spectral_grid.nlayers
            @test nt == Int(period / output.interval) + 1
        end
    end
end

@testset "Output for ShallowWaterModel" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
    period = Day(1)

    for output_NF in (Float32, Float64)
        spectral_grid = SpectralGrid(nlayers = 1)
        model = ShallowWaterModel(spectral_grid)
        output = NetCDFOutput(model, path = tmp_output_path; output_NF)
        simulation = initialize!(model; output)
        run!(simulation, output = true; period)
        @test simulation.feedback.nans_detected == false

        # read netcdf file and check that all variables exist
        ds = NCDataset(joinpath(output.run_path, output.filename))
        for key in keys(output.variables)
            @test haskey(ds, key)

            # test output number format, NCDatasets returns Union{Missing, Float32} for Float32 so do <: to check
            @test output_NF <: eltype(ds[key][:])
        end

        # add divergence output
        div_output = SpeedyWeather.DivergenceOutput()
        add!(output.variables, div_output)
        run!(simulation, output = true; period)
        ds = NCDataset(joinpath(output.run_path, output.filename))
        @test haskey(ds, div_output.name)

        # add orography output
        orog_output = SpeedyWeather.OrographyOutput()
        add!(output.variables, orog_output)
        run!(simulation, output = true; period)
        ds = NCDataset(joinpath(output.run_path, output.filename))
        @test haskey(ds, orog_output.name)

        nx, ny = size(ds[orog_output.name])
        @test (nx, ny) == RingGrids.matrix_size(output.field2D)

        # delete divergence output
        delete!(output, div_output.name)
        run!(simulation, output = true; period)
        ds = NCDataset(joinpath(output.run_path, output.filename))
        @test ~haskey(ds, div_output.name)
    end
end

@testset "Output for PrimitiveDryModel" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
    period = Day(1)

    # test also output at various resolutions
    for nlat_half in (24, 32, 48, 64)
        spectral_grid = SpectralGrid(nlayers = 8)
        output_grid = RingGrids.full_grid_type(typeof(spectral_grid.grid))(nlat_half)
        model = PrimitiveDryModel(spectral_grid)
        output = NetCDFOutput(model, path = tmp_output_path; output_grid)
        simulation = initialize!(model; output)
        run!(simulation, output = true; period)
        @test simulation.feedback.nans_detected == false

        # read netcdf file and check that all variables exist
        ds = NCDataset(joinpath(output.run_path, output.filename))
        for key in keys(output.variables)
            @test haskey(ds, key)
        end

        nx, ny, nz, nt = size(ds[:vor])
        @test nx == 2ny
        @test ny == 2nlat_half
    end
end

@testset "Output for PrimitiveWetModel" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
    period = Day(1)

    # test also output at various resolutions
    spectral_grid = SpectralGrid(nlayers = 8)
    for interval in (Hour(1), Minute(120), Hour(3), Hour(6), Day(1))
        model = PrimitiveWetModel(spectral_grid)
        output = NetCDFOutput(model, path = tmp_output_path; interval)
        simulation = initialize!(model; output)
        run!(simulation, output = true; period)
        @test simulation.feedback.nans_detected == false

        # read netcdf file and check that all variables exist
        ds = NCDataset(joinpath(output.run_path, output.filename))
        for key in keys(output.variables)
            @test haskey(ds, key)

            # test time
            nt = size(ds[key])[end]
            @test nt == Int(period / output.interval) + 1
        end
    end

    # test outputting other model defaults
    model = PrimitiveWetModel(spectral_grid)
    output = NetCDFOutput(model, path = tmp_output_path)

    # Add surface variables output for testing them
    add!(
        output,
        SpeedyWeather.ZonalVelocity10mOutput(),
        SpeedyWeather.MeridionalVelocity10mOutput(),
        SpeedyWeather.SurfaceTemperatureOutput(),
        #         SpeedyWeather.MeanSeaLevelPressureOutput(),   # this should be default now
        #         SpeedyWeather.SurfacePressureOutput(),        # don't output surface pressure too
    )

    simulation = initialize!(model; output)
    run!(simulation, output = true; period)

    @test simulation.feedback.nans_detected == false
    ds = NCDataset(joinpath(output.run_path, output.filename))

    # test
    @test haskey(ds, "temp")
    @test haskey(ds, "humid")
    @test ~haskey(ds, "pres")   # with MSLP as default this should not be contained in the nc file
    @test haskey(ds, "mslp")    # but this variable

    # Test reasonable scale for mean
    p₀ = model.atmosphere.reference_pressure / 100      # Pa -> hPa
    mslp = ds["mslp"].var[:, :, end]    # variable at last time step `.var` to read the raw data ignoring any mask

    # should be within ~800 to ~1200hPa
    @test all(0.8 .< mslp ./ p₀ .< 1.2)

    ## test u10, v10 existence
    @test haskey(ds, "u")
    @test haskey(ds, "u10")
    @test haskey(ds, "v")
    @test haskey(ds, "v10")

    # and 10m values cannot exceed lowermost layer of u,v value
    @test maximum(abs.(ds["u10"].var[:, :, end])) < maximum(abs.(ds["u"].var[:, :, end, end]))
    @test maximum(abs.(ds["v10"].var[:, :, end])) < maximum(abs.(ds["u"].var[:, :, end, end]))

    ## surface temperature should be within 60-130% of
    T₀ = model.atmosphere.reference_temperature     # in K
    Tsurf = ds["tsurf"].var[:, :, end] .+ 273.15    # last timestep from ˚C to K
    @test all(0.6 .< (Tsurf ./ T₀) .< 1.3)
end

@testset "Restart from restart file" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits

    spectral_grid = SpectralGrid()
    model = PrimitiveDryModel(spectral_grid)
    output = NetCDFOutput(model, path = tmp_output_path, id = "restart-test")
    simulation = initialize!(model; output)
    run!(simulation, output = true; period = Day(1))

    initial_conditions = StartFromFile(path = tmp_output_path, id = "restart-test")
    model_new = PrimitiveDryModel(spectral_grid; initial_conditions)
    simulation_new = initialize!(model_new)

    progn_old = simulation.variables.prognostic
    progn_new = simulation_new.variables.prognostic

    for varname in (:vorticity, :divergence, :temperature, :pressure)
        var_old = getfield(progn_old, varname)
        var_new = getfield(progn_new, varname)
        @test all(var_old .== var_new)
    end
end

@testset "Restart from restart file without output" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_restart_")  # Cleaned up when the process exits

    spectral_grid = SpectralGrid()
    model = PrimitiveDryModel(spectral_grid)
    callbacks = SpeedyWeather.CallbackDict(
        :restart_file => WriteVariablesRestartFile(path = tmp_output_path, filename = "myrestart.jld2"),
    )
    simulation = initialize!(model; callbacks)
    run!(simulation, period = Day(1))

    initial_conditions = StartFromFile(run_folder = tmp_output_path, filename = "myrestart.jld2")
    model_new = PrimitiveDryModel(spectral_grid; initial_conditions)
    simulation_new = initialize!(model_new)

    progn_old = simulation.variables.prognostic
    progn_new = simulation_new.variables.prognostic

    for varname in (:vorticity, :divergence, :temperature, :pressure)
        var_old = getfield(progn_old, varname)
        var_new = getfield(progn_new, varname)
        @test all(var_old .== var_new)
    end
end

@testset "Time axis" begin

    function manual_time_axis(startdate, dt, n_timesteps)
        time_axis = zeros(typeof(startdate), n_timesteps + 1)
        for i in 0:n_timesteps
            time_axis[i + 1] = startdate + dt * i
        end
        time_axis
    end

    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits

    spectral_grid = SpectralGrid()
    model = PrimitiveDryModel(spectral_grid)
    output = NetCDFOutput(model, path = tmp_output_path, id = "dense-output-test", interval = Hour(0))
    simulation = initialize!(model; output)
    run!(simulation, output = true; period = Day(1))

    progn = simulation.variables.prognostic
    tmp_read_path = joinpath(output.run_path, output.filename)
    t = NCDataset(tmp_read_path)["time"][:]
    @test t == manual_time_axis(output.startdate, model.time_stepping.Δt_millisec, progn.clock.n_timesteps)

    # do a simulation with the adjust_Δt_with_output turned on
    time_stepping = Leapfrog(spectral_grid, adjust_with_output = true)
    model = PrimitiveDryModel(spectral_grid; time_stepping)
    output = NetCDFOutput(model, path = tmp_output_path, id = "adjust_dt_with_output-test", interval = Minute(70))
    simulation = initialize!(model; output)
    run!(simulation, output = true; period = Day(1))
    t = SpeedyWeather.load_trajectory("time", simulation)
    @test all(y -> y == diff(t)[1], diff(t)) # all elements equal
    @test diff(t)[1] == Millisecond(Minute(70))

    # this is a nonsense simulation with way too large timesteps, but it's here to test the time axis output
    # for future tests: This simulation blows up because of too large time steps but only a warning is thrown
    # at the moment, no error. With the refactored simulation/output split, `later_timestep!`
    # short-circuits once NaNs are detected so output also stops at that point — so the written
    # time axis is a prefix of the full requested axis, not the full one.
    # 1kyrs simulation
    spectral_grid = SpectralGrid()
    time_stepping = Leapfrog(spectral_grid, Δt_at_T31 = Day(3650))
    model = PrimitiveDryModel(spectral_grid; time_stepping)
    output = NetCDFOutput(model, path = tmp_output_path, id = "long-output-test", interval = Day(3650))
    simulation = initialize!(model; output)
    model.implicit.initialized = true     # bypass implicit initialization with this insanely large timestep (threw SingularException)
    model.implicit.reinitialize = false
    run!(simulation, output = true, period = Day(365000))

    progn = simulation.variables.prognostic
    tmp_read_path = joinpath(output.run_path, output.filename)
    t = NCDataset(tmp_read_path)["time"][:]
    expected_full = manual_time_axis(output.startdate, model.time_stepping.Δt_millisec, progn.clock.n_timesteps)
    # what we got is a prefix of the full requested time axis
    @test t == expected_full[1:length(t)]
    @test t == SpeedyWeather.load_trajectory("time", simulation)
end

@testset "get_output_path" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")

    # output inactive: should throw an error
    spectral_grid = SpectralGrid(nlayers = 1)
    model = BarotropicModel(spectral_grid)
    simulation = initialize!(model)
    @test_throws ErrorException SpeedyWeather.get_output_path(simulation)

    # output active: should return the correct path
    model = BarotropicModel(spectral_grid)
    output = NetCDFOutput(model, path = tmp_output_path)
    simulation = initialize!(model; output)
    run!(simulation, output = true, period = Day(1))
    expected_path = joinpath(output.run_path, output.filename)
    @test SpeedyWeather.get_output_path(simulation) == expected_path
    @test isfile(SpeedyWeather.get_output_path(simulation))
end

@testset "set!(output, model; interval)" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")

    spectral_grid = SpectralGrid(nlayers = 1)
    model = BarotropicModel(spectral_grid)
    output = NetCDFOutput(model, path = tmp_output_path, interval = Hour(6))
    simulation = initialize!(model; output)
    run!(simulation, output = true, period = Day(1))

    Δt_ms = model.time_stepping.Δt_millisec.value

    # change interval to a clean multiple of the model time step (no rounding)
    new_interval_ms = 4 * Δt_ms
    new_interval = Millisecond(new_interval_ms)
    @test_logs set!(output, model, simulation.callbacks; interval = new_interval)   # no @info log
    @test output.core.output_every_n_steps == 4
    @test Millisecond(output.interval).value == new_interval_ms

    # counters reset by re-initialization
    @test output.core.timestep_counter == 0
    @test output.core.output_counter == 0

    # change to a different multiple
    set!(output, model, simulation.callbacks; interval = Millisecond(2 * Δt_ms))
    @test output.core.output_every_n_steps == 2
    @test Millisecond(output.interval).value == 2 * Δt_ms

    # request an interval that is NOT a multiple of the model time step:
    # use Δt + 1ms to force rounding back to a single time step (and an @info)
    non_multiple = Millisecond(Δt_ms + 1)
    @test_logs (:info,) set!(output, model, simulation.callbacks; interval = non_multiple)
    @test output.core.output_every_n_steps == 1
    @test Millisecond(output.interval).value == Δt_ms

    # check that run! after set! uses the new output frequency
    period = Day(1)
    set!(output, model, simulation.callbacks; interval = Hour(2))   # = 3 * Δt at default settings
    expected_n = round(Int, Millisecond(Hour(2)).value / Δt_ms)
    @test output.core.output_every_n_steps == expected_n
    run!(simulation, output = true; period)
    ds = NCDataset(joinpath(output.run_path, output.filename))
    nt = size(ds["vor"])[end]
    @test nt == Int(Millisecond(period).value ÷ Millisecond(output.interval).value) + 1
end
