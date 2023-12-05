using NCDatasets, Dates

@testset "Output on various grids" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
    period = Day(1)

    # default grid, Float64, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float64,nlev=1)
    output = OutputWriter(spectral_grid,ShallowWater,path=tmp_output_path)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;period)
    @test simulation.model.feedback.nars_detected == false

    # default grid, Float32, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float32,nlev=1)
    output = OutputWriter(spectral_grid,ShallowWater,path=tmp_output_path)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;period)
    @test simulation.model.feedback.nars_detected == false

    # FullClenshawGrid, Float32, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float32,Grid=FullClenshawGrid,nlev=1)
    output = OutputWriter(spectral_grid,ShallowWater,path=tmp_output_path)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;period)
    @test simulation.model.feedback.nars_detected == false

    # OctahedralClenshawGrid, Float32, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float32,Grid=OctahedralClenshawGrid,nlev=1)
    output = OutputWriter(spectral_grid,ShallowWater,path=tmp_output_path)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;period)
    @test simulation.model.feedback.nars_detected == false

    # HEALPixGrid, Float32, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float32,Grid=HEALPixGrid,nlev=1)
    output = OutputWriter(spectral_grid,ShallowWater,path=tmp_output_path)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;period)
    @test simulation.model.feedback.nars_detected == false

    # OctaHEALPixGrid, Float32, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float32,Grid=OctaHEALPixGrid,nlev=1)
    output = OutputWriter(spectral_grid,ShallowWater,path=tmp_output_path)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;period)
    @test simulation.model.feedback.nars_detected == false

    # OctahedralClenshawGrid, as matrix, Float32, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float32,Grid=OctahedralClenshawGrid,nlev=1)
    output = OutputWriter(spectral_grid,ShallowWater,path=tmp_output_path,as_matrix=true)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;period)
    @test simulation.model.feedback.nars_detected == false

    # OctaHEALPixGrid, as matrix, Float32, PrimitiveDry
    spectral_grid = SpectralGrid(;NF=Float32,Grid=OctaHEALPixGrid)
    output = OutputWriter(spectral_grid,PrimitiveDry,path=tmp_output_path,as_matrix=true)
    model = PrimitiveDryModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;period)
    @test simulation.model.feedback.nars_detected == false

    # OctaHEALPixGrid, as matrix, Float32, but output Float64 PrimitiveDry
    spectral_grid = SpectralGrid(;NF=Float32,Grid=OctaHEALPixGrid)
    output = OutputWriter(spectral_grid,PrimitiveDry,path=tmp_output_path,as_matrix=true,NF=Float64)
    model = PrimitiveDryModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;period)
    @test simulation.model.feedback.nars_detected == false
end

@testset "Restart from output file" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits

    spectral_grid = SpectralGrid()
    output = OutputWriter(spectral_grid,PrimitiveDry,path=tmp_output_path,id="restart-test")
    model = PrimitiveDryModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;period=Day(1))

    initial_conditions = StartFromFile(path=tmp_output_path,id="restart-test")
    model2 = PrimitiveDryModel(;spectral_grid,initial_conditions)
    simulation2 = initialize!(model2)

    p1 = simulation.prognostic_variables
    p2 = simulation2.prognostic_variables

    for varname in propertynames(p1.layers[1].timesteps[1])
        if SpeedyWeather.has(p1, varname)
            for (var_new, var_old) in zip(SpeedyWeather.get_var(p1, varname), SpeedyWeather.get_var(p2, varname))
                @test all(var_new .== var_old)
            end
        end
    end
    @test all(SpeedyWeather.get_pressure(p1) .== SpeedyWeather.get_pressure(p2))
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
    output = OutputWriter(spectral_grid,PrimitiveDry,path=tmp_output_path,id="dense-output-test",output_dt=Hour(0))
    model = PrimitiveDryModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;period=Day(1))
    
    progn = simulation.prognostic_variables
    tmp_read_path = joinpath(model.output.run_path,model.output.filename)
    t = NCDataset(tmp_read_path)["time"][:]
    @test t == manual_time_axis(model.output.startdate,model.time_stepping.Δt_millisec,progn.clock.n_timesteps)
    
    # do a simulation with the adjust_Δt_with_output turned on 
    output = OutputWriter(spectral_grid,PrimitiveDry,path=tmp_output_path,id="adjust_dt_with_output-test",output_dt=Minute(70))
    time_stepping = Leapfrog(spectral_grid, adjust_with_output=true)
    model = PrimitiveDryModel(;spectral_grid,output,time_stepping)
    simulation = initialize!(model)
    run!(simulation,output=true;period=Day(1))
    t = SpeedyWeather.load_trajectory("time", model)
    @test all(y->y==diff(t)[1], diff(t)) # all elements equal 
    @test diff(t)[1] == Minute(70)

    # this is a nonsense simulation with way too large timesteps, but it's here to test the time axis output
    # for future tests: This simulation blows up because of too large time steps but only a warning is thrown
    # at the moment, no error
    # 1kyrs simulation
    spectral_grid = SpectralGrid()
    time_stepping = Leapfrog(spectral_grid,Δt_at_T31=Day(3650))
    output = OutputWriter(spectral_grid,PrimitiveDry,path=tmp_output_path,id="long-output-test",output_dt=Day(3650))
    model = PrimitiveDryModel(;spectral_grid,output,time_stepping)
    simulation = initialize!(model)
    run!(simulation,output=true,period=Day(365000))

    progn = simulation.prognostic_variables
    tmp_read_path = joinpath(model.output.run_path,model.output.filename)
    t = NCDataset(tmp_read_path)["time"][:]
    @test t == manual_time_axis(model.output.startdate,model.time_stepping.Δt_millisec,progn.clock.n_timesteps)
    @test t == SpeedyWeather.load_trajectory("time", model)
end
