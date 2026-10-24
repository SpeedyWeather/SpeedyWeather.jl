@testset "Progress and parameter without output" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_feedback_")  # Cleaned up when the process exits

    spectral_grid = SpectralGrid(nlayers = 1)
    model = BarotropicModel(spectral_grid)
    model.feedback.verbose = false
    simulation = initialize!(model)

    add!(model, ProgressTxt(path = tmp_output_path, write_only_with_output = true))
    add!(model, ParametersTxt(path = tmp_output_path, write_only_with_output = true))

    run!(simulation, period = Day(1))

    # test that files are not created because output=false
    @test ~isfile(joinpath(tmp_output_path, "parameters.txt"))
    @test ~isfile(joinpath(tmp_output_path, "progress.txt"))

    add!(model, ProgressTxt(path = tmp_output_path, write_only_with_output = false))
    add!(model, ParametersTxt(path = tmp_output_path, write_only_with_output = false))

    run!(simulation, period = Day(1))

    # test that files are created even if output=false
    @test isfile(joinpath(tmp_output_path, "parameters.txt"))
    @test isfile(joinpath(tmp_output_path, "progress.txt"))
end

@testset "Feedback diagnostics opt-in for Reactant" begin
    show_diagnostic = SpeedyWeather.show_diagnostic
    reactant = SpeedyWeather.ReactantDevice(nothing)   # for dispatch only, no Reactant needed
    for arch in (SpeedyWeather.CPU(), reactant), option in (true, false)
        @test show_diagnostic(option, arch) == option   # explicit options apply everywhere
    end
    @test show_diagnostic(nothing, SpeedyWeather.CPU())   # default: shown
    @test !show_diagnostic(nothing, reactant)              # default for Reactant: not shown, opt-in

    # the default feedback leaves it to the architecture
    feedback = Feedback()
    @test isnothing(feedback.show_umax)
    @test isnothing(feedback.show_temperature_range)
end
