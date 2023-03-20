@testset "NetCDF output" begin 
    import NetCDF

    @testset "Time axis" begin 

        function manual_time_axis(dt, n_timesteps)
            
            time_axis = zeros(Float64, n_timesteps+1)
            for i=0:n_timesteps
                time_axis[i+1] = Float64(dt) * i
            end 
            time_axis 
        end 

        p, d, m = initialize_speedy(Float32, output=true, n_days=1, output_dt=0, run_id="dense-output-test")
        SpeedyWeather.time_stepping!(p, d, m)

        t = NetCDF.ncread("run-dense-output-test/output.nc", "time")
        @test t ≈ Int64.(manual_time_axis(m.constants.Δt_sec, m.constants.n_timesteps))
        
        # this is a nonsense simulation with way too large timesteps, but it's here to test the time axis output
        # 1kyrs simulation         
        p, d, m = initialize_speedy(Float32, output=true, n_days=365000, Δt_at_T31=60*24*365*10, output_dt=24*365*10, run_id="long-output-test")
        SpeedyWeather.time_stepping!(p, d, m)

        t = NetCDF.ncread("run-long-output-test/output.nc", "time")
        @test t ≈ Int64.(manual_time_axis(m.constants.Δt_sec, m.constants.n_timesteps))

        @test t ≈ SpeedyWeather.load_trajectory("time", m)
    end 
end 


