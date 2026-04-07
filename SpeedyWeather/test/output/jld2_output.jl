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
    @test all(isfinite, ic.vorticity)
    @test all(isfinite, ic.divergence)
    @test all(isfinite, ic.humidity)
    @test all(isfinite, ic.pressure)
    @test all(isfinite, ic.temperature)

    # check last output also contains finite values
    final_output = f["output_vector"][end].prognostic
    @test all(isfinite, final_output.vorticity)
    @test all(isfinite, final_output.divergence)
    @test all(isfinite, final_output.humidity)
    @test all(isfinite, final_output.pressure)
    @test all(isfinite, final_output.temperature)
end
