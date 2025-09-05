using NCDatasets, Dates

@testset "Output for BarotropicModel" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
    period = Day(1)

    for Grid in (FullGaussianGrid, FullClenshawGrid, OctahedralGaussianGrid, OctahedralClenshawGrid, HEALPixGrid, OctaHEALPixGrid)
        spectral_grid = SpectralGrid(nlayers=1)
        output = NetCDFOutput(spectral_grid, path=tmp_output_path)
        model = BarotropicModel(spectral_grid; output)
        simulation = initialize!(model)
        run!(simulation, output=true; period)
        @test simulation.model.feedback.nans_detected == false

        # read netcdf file and check that all variables exist
        ds = NCDataset(joinpath(model.output.run_path, model.output.filename))
        for key in keys(output.variables)
            @test haskey(ds, key)

            # test dimensions
            nx, ny, nz, nt = size(ds[key])
            @test (nx, ny) == RingGrids.matrix_size(output.field2D)
            @test nz == spectral_grid.nlayers
            @test nt == Int(period / output.output_dt) + 1
        end
    end
end

@testset "Output for ShallowWaterModel" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
    period = Day(1)

    for output_NF in (Float32, Float64)
        spectral_grid = SpectralGrid(nlayers=1)
        output = NetCDFOutput(spectral_grid, ShallowWater, path=tmp_output_path; output_NF)
        model = ShallowWaterModel(spectral_grid; output)
        simulation = initialize!(model)
        run!(simulation, output=true; period)
        @test simulation.model.feedback.nans_detected == false

        # read netcdf file and check that all variables exist
        ds = NCDataset(joinpath(model.output.run_path, model.output.filename))
        for key in keys(output.variables)
            @test haskey(ds, key)

            # test output number format, NCDatasets returns Union{Missing, Float32} for Float32 so do <: to check
            @test output_NF <: eltype(ds[key][:])
        end

        # add divergence output
        div_output = SpeedyWeather.DivergenceOutput()
        add!(output.variables, div_output)
        run!(simulation, output=true; period)
        ds = NCDataset(joinpath(model.output.run_path, model.output.filename))
        @test haskey(ds, div_output.name)

        # add orography output
        orog_output = SpeedyWeather.OrographyOutput()
        add!(output.variables, orog_output)
        run!(simulation, output=true; period)
        ds = NCDataset(joinpath(model.output.run_path, model.output.filename))
        @test haskey(ds, orog_output.name)
        
        nx, ny = size(ds[orog_output.name])
        @test (nx, ny) == RingGrids.matrix_size(output.field2D)

        # delete divergence output
        delete!(output, div_output.name)
        run!(simulation, output=true; period)
        ds = NCDataset(joinpath(model.output.run_path, model.output.filename))
        @test ~haskey(ds, div_output.name)
    end
end

@testset "Output for PrimitiveDryModel" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
    period = Day(1)

    # test also output at various resolutions
    for nlat_half in (24, 32, 48, 64)
        spectral_grid = SpectralGrid(nlayers=8)
        output_grid = RingGrids.full_grid_type(typeof(spectral_grid.grid))(nlat_half)
        output = NetCDFOutput(spectral_grid, ShallowWater, path=tmp_output_path; output_grid)
        model = PrimitiveDryModel(spectral_grid; output)
        simulation = initialize!(model)
        run!(simulation, output=true; period)
        @test simulation.model.feedback.nans_detected == false

        # read netcdf file and check that all variables exist
        ds = NCDataset(joinpath(model.output.run_path, model.output.filename))
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
    spectral_grid = SpectralGrid(nlayers=8)
    for output_dt in (Hour(1), Minute(120), Hour(3), Hour(6), Day(1))
        output = NetCDFOutput(spectral_grid, PrimitiveWet, path=tmp_output_path; output_dt)
        model = PrimitiveWetModel(spectral_grid; output)
        simulation = initialize!(model)
        run!(simulation, output=true; period)
        @test simulation.model.feedback.nans_detected == false

        # read netcdf file and check that all variables exist
        ds = NCDataset(joinpath(model.output.run_path, model.output.filename))
        for key in keys(output.variables)
            @test haskey(ds, key)
           
            # test time
            nt = size(ds[key])[end]
            @test nt == Int(period / output.output_dt) + 1
        end
    end

    # test outputting other model defaults
    output = NetCDFOutput(spectral_grid, Barotropic, path=tmp_output_path)
    model = PrimitiveWetModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output=true; period)
    @test simulation.model.feedback.nans_detected == false
    ds = NCDataset(joinpath(model.output.run_path, model.output.filename))
    @test ~haskey(ds, "temp")
    @test ~haskey(ds, "humid")
    @test ~haskey(ds, "pres")
end

@testset "Restart from output file" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits

    spectral_grid = SpectralGrid()
    output = NetCDFOutput(spectral_grid, PrimitiveDry, path=tmp_output_path, id="restart-test")
    model = PrimitiveDryModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output=true; period=Day(1))

    initial_conditions = StartFromFile(path=tmp_output_path, id="restart-test")
    model_new = PrimitiveDryModel(spectral_grid; initial_conditions)
    simulation_new = initialize!(model_new)

    progn_old = simulation.prognostic_variables
    progn_new = simulation_new.prognostic_variables

    for varname in (:vor, :div, :temp, :pres)
        var_old = getfield(progn_old, varname)[1]
        var_new = getfield(progn_new, varname)[1]
        @test all(var_old .== var_new)
    end
end 

@testset "Time axis" begin 

    function manual_time_axis(startdate, dt, n_timesteps)
        time_axis = zeros(typeof(startdate), n_timesteps+1)
        for i=0:n_timesteps
            time_axis[i+1] = startdate + dt*i
        end 
        time_axis 
    end 

    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
    
    spectral_grid = SpectralGrid()
    output = NetCDFOutput(spectral_grid, PrimitiveDry, path=tmp_output_path, id="dense-output-test", output_dt=Hour(0))
    model = PrimitiveDryModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output=true; period=Day(1))
    
    progn = simulation.prognostic_variables
    tmp_read_path = joinpath(model.output.run_path, model.output.filename)
    t = NCDataset(tmp_read_path)["time"][:]
    @test t == manual_time_axis(model.output.startdate, model.time_stepping.Δt_millisec, progn.clock.n_timesteps)
    
    # do a simulation with the adjust_Δt_with_output turned on 
    output = NetCDFOutput(spectral_grid, PrimitiveDry, path=tmp_output_path, id="adjust_dt_with_output-test", output_dt=Minute(70))
    time_stepping = Leapfrog(spectral_grid, adjust_with_output=true)
    model = PrimitiveDryModel(spectral_grid; output, time_stepping)
    simulation = initialize!(model)
    run!(simulation, output=true; period=Day(1))
    t = SpeedyWeather.load_trajectory("time", model)
    @test all(y->y==diff(t)[1], diff(t)) # all elements equal 
    @test diff(t)[1] == Minute(70)

    # this is a nonsense simulation with way too large timesteps, but it's here to test the time axis output
    # for future tests: This simulation blows up because of too large time steps but only a warning is thrown
    # at the moment, no error
    # 1kyrs simulation
    spectral_grid = SpectralGrid()
    time_stepping = Leapfrog(spectral_grid, Δt_at_T31=Day(3650))
    output = NetCDFOutput(spectral_grid, PrimitiveDry, path=tmp_output_path, id="long-output-test", output_dt=Day(3650))
    model = PrimitiveDryModel(spectral_grid; output, time_stepping)
    simulation = initialize!(model)
    run!(simulation, output=true, period=Day(365000))

    progn = simulation.prognostic_variables
    tmp_read_path = joinpath(model.output.run_path, model.output.filename)
    t = NCDataset(tmp_read_path)["time"][:]
    @test t == manual_time_axis(model.output.startdate, model.time_stepping.Δt_millisec, progn.clock.n_timesteps)
    @test t == SpeedyWeather.load_trajectory("time", model)
end
