using JLD2

@testset "JLD2 Output" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_jld2tests_")  # Cleaned up when the process exits

    spectral_grid = SpectralGrid()

    # write-restart false is important to not mutate the final state in the simulation object
    output = JLD2Output(output_dt = Hour(1), path = tmp_output_path, id = "jld2-test", write_restart = false)
    model = PrimitiveWetModel(; spectral_grid, output)
    simulation = initialize!(model)
    run!(simulation, period = Day(2), output = true)

    f = jldopen(joinpath(output.run_path, output.filename), "r")

    @test length(f["output_vector"]) == 2 * 24 + 1 # 2 Days with 24 Hours + 1 for the IC

    # check that IC output contains finite values for all prognostic variables
    ic = f["output_vector"][1].prognostic
    @test all(isfinite, ic.vor)
    @test all(isfinite, ic.div)
    @test all(isfinite, ic.humid)
    @test all(isfinite, ic.pres)
    @test all(isfinite, ic.temp)

    # check last output also contains finite values
    final_output = f["output_vector"][end].prognostic
    @test all(isfinite, final_output.vor)
    @test all(isfinite, final_output.div)
    @test all(isfinite, final_output.humid)
    @test all(isfinite, final_output.pres)
    @test all(isfinite, final_output.temp)
end
