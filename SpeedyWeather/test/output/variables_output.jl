using JLD2

@testset "JLD2 Output" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_jld2tests_")  # Cleaned up when the process exits

    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)

    # write-restart false is important to not mutate the final state in the simulation object
    output = JLD2Output(output_dt = Hour(6), path = tmp_output_path, id = "jld2-test", write_restart = false)
    model = BarotropicModel(; spectral_grid, output)
    simulation = initialize!(model)
    run!(simulation, period = Day(1), output = true)

    f = jldopen(joinpath(output.run_path, output.filename), "r")

    @test length(f["output_vector"]) == 4 + 1 # 1 Day with 4 6-hourly outputs + 1 for the IC

    # check that IC output contains finite values for the prognostic vorticity
    ic = f["output_vector"][1].prognostic
    @test all(isfinite, ic.vorticity)

    # check last output also contains finite values
    final_output = f["output_vector"][end].prognostic
    @test all(isfinite, final_output.vorticity)

    # the final snapshot was saved before unscale!, so vorticity is scaled by radius;
    # after run! the simulation's vorticity has been unscaled
    radius = simulation.model.planet.radius
    @test Array(final_output.vorticity) ≈ Array(simulation.variables.prognostic.vorticity) * radius

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

    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)

    # only save prognostic and grid groups
    output = JLD2Output(
        output_dt = Hour(6),
        path = tmp_output_path,
        id = "jld2-groups",
        write_restart = false,
        groups = (:prognostic, :grid),
    )
    model = BarotropicModel(; spectral_grid, output)
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
    @test all(isfinite, snapshot.prognostic.vorticity)
    close(f)
end

@testset "ArrayOutput" begin
    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)

    output = ArrayOutput(output_dt = Hour(6))
    model = BarotropicModel(; spectral_grid, output)
    simulation = initialize!(model)
    run!(simulation, period = Day(1), output = true)

    @test length(output.output) == 4 + 1  # 1 Day with 4 6-hourly outputs + 1 for the IC

    # check IC snapshot
    ic = output.output[1].prognostic
    @test all(isfinite, ic.vorticity)

    # check final snapshot
    final_snapshot = output.output[end].prognostic
    @test all(isfinite, final_snapshot.vorticity)

    # the final snapshot was saved before unscale!, so vorticity is scaled by radius;
    # after run! the simulation's vorticity has been unscaled
    radius = simulation.model.planet.radius
    @test Array(final_snapshot.vorticity) ≈ Array(simulation.variables.prognostic.vorticity) * radius

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
    @test output.output[1].prognostic.vorticity !== output.output[end].prognostic.vorticity
end

@testset "ArrayOutput with group selection" begin
    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)

    # only save prognostic group
    output = ArrayOutput(output_dt = Hour(6), groups = (:prognostic,))
    model = BarotropicModel(; spectral_grid, output)
    simulation = initialize!(model)
    run!(simulation, period = Day(1), output = true)

    @test length(output.output) == 4 + 1  # 4 outputs at 6h intervals + IC

    snapshot = output.output[1]
    @test haskey(snapshot, :prognostic)
    @test !haskey(snapshot, :grid)
    @test !haskey(snapshot, :tendencies)
    @test !haskey(snapshot, :dynamics)

    # prognostic data is valid
    @test all(isfinite, snapshot.prognostic.vorticity)
end
