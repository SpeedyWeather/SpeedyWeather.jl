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

    @test f["output_vector"][1][1].vor == progn_ic.vor
    @test f["output_vector"][1][1].div == progn_ic.div
    @test f["output_vector"][1][1].humid == progn_ic.humid
    @test f["output_vector"][1][1].pres == progn_ic.pres
    @test f["output_vector"][1][1].temp == progn_ic.temp

    # the last output in the simuluation is unscaled, shifted and compressed by `write_restart_file!` 
    # for comparision we need to do that as well
    final_output = f["output_vector"][end][1] 
    SpeedyWeather.unscale!(final_output)
    
    @test final_output.vor == simulation.prognostic_variables.vor
    @test final_output.div == simulation.prognostic_variables.div
    @test final_output.humid == simulation.prognostic_variables.humid
    @test final_output.pres == simulation.prognostic_variables.pres
    @test final_output.temp == simulation.prognostic_variables.temp
end 