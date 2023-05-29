import NetCDF

@testset "Output on various grids" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
    n_days = 1

    # default grid, Float64, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float64)
    output = OutputWriter(spectral_grid,path=tmp_output_path)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;n_days)
    @test simulation.model.feedback.nars_detected == false

    # default grid, Float32, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float32)
    output = OutputWriter(spectral_grid,path=tmp_output_path)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;n_days)
    @test simulation.model.feedback.nars_detected == false

    # FullClenshawGrid, Float32, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float32,Grid=FullClenshawGrid)
    output = OutputWriter(spectral_grid,path=tmp_output_path)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;n_days)
    @test simulation.model.feedback.nars_detected == false

    # OctahedralClenshawGrid, Float32, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float32,Grid=OctahedralClenshawGrid)
    output = OutputWriter(spectral_grid,path=tmp_output_path)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;n_days)
    @test simulation.model.feedback.nars_detected == false

    # HEALPixGrid, Float32, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float32,Grid=HEALPixGrid)
    output = OutputWriter(spectral_grid,path=tmp_output_path)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;n_days)
    @test simulation.model.feedback.nars_detected == false

    # OctaHEALPixGrid, Float32, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float32,Grid=OctaHEALPixGrid)
    output = OutputWriter(spectral_grid,path=tmp_output_path)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;n_days)
    @test simulation.model.feedback.nars_detected == false

    # OctahedralClenshawGrid, as matrix, Float32, ShallowWater
    spectral_grid = SpectralGrid(;NF=Float32,Grid=OctahedralClenshawGrid)
    output = OutputWriter(spectral_grid,path=tmp_output_path,as_matrix=true)
    model = ShallowWaterModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;n_days)
    @test simulation.model.feedback.nars_detected == false

    # OctaHEALPixGrid, as matrix, Float32, PrimitiveDry
    spectral_grid = SpectralGrid(;NF=Float32,Grid=OctaHEALPixGrid)
    output = OutputWriter(spectral_grid,path=tmp_output_path,as_matrix=true)
    model = PrimitiveDryModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;n_days)
    @test simulation.model.feedback.nars_detected == false

    # OctaHEALPixGrid, as matrix, Float32, but output Float64 PrimitiveDry
    spectral_grid = SpectralGrid(;NF=Float32,Grid=OctaHEALPixGrid)
    output = OutputWriter(spectral_grid,path=tmp_output_path,as_matrix=true,NF=Float64)
    model = PrimitiveDryModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true;n_days)
    @test simulation.model.feedback.nars_detected == false
end

@testset "Restart from output file" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
    
    spectral_grid = SpectralGrid()
    output = OutputWriter(spectral_grid,path=tmp_output_path,id="restart-test")
    model = PrimitiveDryModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true,n_days=1)

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

    function manual_time_axis(dt, n_timesteps)
        
        time_axis = zeros(Float64, n_timesteps+1)
        for i=0:n_timesteps
            time_axis[i+1] = Float64(dt) * i
        end 
        time_axis 
    end 

    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits
    
    spectral_grid = SpectralGrid()
    output = OutputWriter(spectral_grid,path=tmp_output_path,id="dense-output-test",output_dt=0)
    model = PrimitiveDryModel(;spectral_grid,output)
    simulation = initialize!(model)
    run!(simulation,output=true,n_days=1)
    
    tmp_read_path = joinpath(model.output.run_path,model.output.filename)
    t = NetCDF.ncread(tmp_read_path, "time")
    @test t ≈ Int64.(manual_time_axis(model.time_stepping.Δt_sec, model.clock.n_timesteps))
    
    # this is a nonsense simulation with way too large timesteps, but it's here to test the time axis output
    # for future tests: This simulation blows up because of too large time steps but only a warning is thrown
    # at the moment, no error
    # 1kyrs simulation
    spectral_grid = SpectralGrid()
    time_stepping = Leapfrog(spectral_grid,Δt_at_T31=60*24*365*10)
    output = OutputWriter(spectral_grid,path=tmp_output_path,id="long-output-test",output_dt=24*365*10)
    model = PrimitiveDryModel(;spectral_grid,output,time_stepping)
    simulation = initialize!(model)
    run!(simulation,output=true,n_days=365000)

    tmp_read_path = joinpath(model.output.run_path,model.output.filename)
    t = NetCDF.ncread(tmp_read_path, "time")
    @test t ≈ Int64.(manual_time_axis(model.time_stepping.Δt_sec, model.clock.n_timesteps))
    @test t ≈ SpeedyWeather.load_trajectory("time", model)
end
