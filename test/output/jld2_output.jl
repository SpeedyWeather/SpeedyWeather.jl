using JLD2 

@testset "JLD2 Output" begin 
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits

    spectral_grid = SpectralGrid()   
    
    # write-restart false is important to not mutate the final state in the simulation object
    output = JLD2Output(path=tmp_output_path, id="jld2-test", output_diagnostic=true, write_restart=false)
    model = PrimitiveWetModel(; spectral_grid, output) 
    simulation = initialize!(model)  
    initialize!(simulation)

    progn_ic = deepcopy(simulation.prognostic_variables)
    diagn_ic = deepcopy(simulation.diagnostic_variables)

    run!(simulation, period=Day(2), output=true) 

    f = jldopen(joinpath(output.path,"run_jld2-test/output.jld2"), "r")
    
    @test length(f["output_vector"]) == simulation.prognostic_variables.clock.n_timesteps + 1 # + 1 for the IC

    for i in eachindex(progn_ic.vor)
        @test f["output_vector"][1][1].vor[i] == progn_ic.vor[i]
        @test f["output_vector"][1][1].div[i] == progn_ic.div[i]
        @test f["output_vector"][1][1].humid[i] == progn_ic.humid[i]
        @test f["output_vector"][1][1].pres[i] == progn_ic.pres[i]
        @test f["output_vector"][1][1].temp[i] == progn_ic.temp[i]
    end 

    # the last output in the simuluation is unscaled, shifted and compressed by `write_restart_file!` 
    # for comparision we need to do that as well
    final_output = f["output_vector"][end][1] 
    SpeedyWeather.unscale!(final_output)
    
    for i in eachindex(progn_ic.vor)
        @test final_output.vor[i] == simulation.prognostic_variables.vor[i]
        @test final_output.div[i] == simulation.prognostic_variables.div[i]
        @test final_output.humid[i] == simulation.prognostic_variables.humid[i]
        @test final_output.pres[i] == simulation.prognostic_variables.pres[i]
        @test final_output.temp[i] == simulation.prognostic_variables.temp[i]
    end
end 