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

    # default groups = (:all,) outputs the full Variables struct
    snapshot = f["output_vector"][1]
    @test snapshot isa Variables
    @test hasproperty(snapshot, :prognostic)
    @test hasproperty(snapshot, :grid)
    @test hasproperty(snapshot, :tendencies)
    @test hasproperty(snapshot, :dynamics)
    @test hasproperty(snapshot, :parameterizations)
    @test hasproperty(snapshot, :particles)
    @test hasproperty(snapshot, :scratch)
    close(f)
end

@testset "JLD2 Output with group selection" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_jld2groups_")

    spectral_grid = SpectralGrid()

    # only save prognostic and grid groups
    output = JLD2Output(
        output_dt = Hour(6),
        path = tmp_output_path,
        id = "jld2-groups",
        write_restart = false,
        groups = (:prognostic, :grid),
    )
    model = PrimitiveWetModel(; spectral_grid, output)
    simulation = initialize!(model)
    run!(simulation, period = Day(1), output = true)

    f = jldopen(joinpath(output.run_path, output.filename), "r")
    @test length(f["output_vector"]) == 4 + 1  # 4 outputs at 6h intervals + IC

    snapshot = f["output_vector"][1]
    @test haskey(snapshot, :prognostic)
    @test haskey(snapshot, :grid)
    @test !haskey(snapshot, :tendencies)
    @test !haskey(snapshot, :dynamics)
    @test !haskey(snapshot, :parameterizations)

    # prognostic data is still valid
    @test all(isfinite, snapshot.prognostic.vor)
    close(f)
end

@testset "ArrayOutput" begin
    spectral_grid = SpectralGrid()

    output = ArrayOutput(output_dt = Hour(1))
    model = PrimitiveWetModel(; spectral_grid, output)
    simulation = initialize!(model)
    run!(simulation, period = Day(2), output = true)

    @test length(output.output) == 2 * 24 + 1  # 2 Days with 24 Hours + 1 for the IC

    # check IC snapshot
    ic = output.output[1].prognostic
    @test all(isfinite, ic.vor)
    @test all(isfinite, ic.div)
    @test all(isfinite, ic.humid)
    @test all(isfinite, ic.pres)
    @test all(isfinite, ic.temp)

    # check final snapshot
    final_snapshot = output.output[end].prognostic
    @test all(isfinite, final_snapshot.vor)
    @test all(isfinite, final_snapshot.div)
    @test all(isfinite, final_snapshot.humid)
    @test all(isfinite, final_snapshot.pres)
    @test all(isfinite, final_snapshot.temp)

    # default groups = (:all,) outputs the full Variables struct
    @test output.output[1] isa Variables
    @test hasproperty(output.output[1], :prognostic)
    @test hasproperty(output.output[1], :grid)
    @test hasproperty(output.output[1], :tendencies)
    @test hasproperty(output.output[1], :dynamics)
    @test hasproperty(output.output[1], :parameterizations)
    @test hasproperty(output.output[1], :particles)
    @test hasproperty(output.output[1], :scratch)

    # snapshots should be independent copies (deepcopy), not aliased
    @test output.output[1].prognostic.vor !== output.output[end].prognostic.vor
end

@testset "ArrayOutput with group selection" begin
    spectral_grid = SpectralGrid()

    # only save prognostic group
    output = ArrayOutput(output_dt = Hour(6), groups = (:prognostic,))
    model = PrimitiveWetModel(; spectral_grid, output)
    simulation = initialize!(model)
    run!(simulation, period = Day(1), output = true)

    @test length(output.output) == 4 + 1  # 4 outputs at 6h intervals + IC

    snapshot = output.output[1]
    @test haskey(snapshot, :prognostic)
    @test !haskey(snapshot, :grid)
    @test !haskey(snapshot, :tendencies)
    @test !haskey(snapshot, :dynamics)

    # prognostic data is valid
    @test all(isfinite, snapshot.prognostic.vor)
end

