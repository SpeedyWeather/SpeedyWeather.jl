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

@testset "Feedback progress elements" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 4)

    @testset "no speedstring piracy" begin
        @test !any(m -> m.module === SpeedyWeather, methods(SpeedyWeather.ProgressMeter.speedstring))
        @test !isdefined(SpeedyWeather, :FEEDBACK_UMAX)
    end

    @testset "default elements" begin
        feedback = Feedback(verbose = false)
        elements = default_elements(feedback)
        @test any(e -> e isa SimulationSpeed, elements)
        @test !any(e -> e isa VerticalCourantNumber, elements)
        @test length(default_elements(Feedback(showspeed = false))) == 4
    end

    @testset "elements on redraw, VerticalCourantNumber" begin
        feedback = Feedback(verbose = false, elements = (default_elements(Feedback())..., VerticalCourantNumber()))
        model = PrimitiveDryModel(spectral_grid; feedback)
        simulation = initialize!(model)
        run!(simulation, period = Day(1))
        p = model.feedback.progress_meter
        line = join(map(e -> SpeedyWeather.ProgressMeter.Elements.print_element(e, p), p.elements[5:end]))
        @test occursin("Cᵥ = ", line)
        @test occursin("m/s", line)
        @test occursin("˚C", line)
        courant = SpeedyWeather.vertical_courant_number(
            simulation.variables.dynamics.w, model.geometry.σ_levels_thick, model.time_stepping.Δt
        )
        @test courant > 0 && isfinite(courant)

        # hand-computed for a single column with σ̇ = 1/s at the interfaces, Δσ = 0.25
        w = zeros(Float32, spectral_grid.grid, 4)
        w.data[:, 1:3] .= 1
        @test SpeedyWeather.vertical_courant_number(w, fill(0.25f0, 4), 10) ≈ 40

        # the interface above layer 2 (σ̇ = 1/s) with the thin Δσ = 0.1 of layer 2 dominates
        w.data .= 0
        w.data[:, 1] .= 1
        @test SpeedyWeather.vertical_courant_number(w, Float32[0.5, 0.1, 0.2, 0.2], 10) ≈ 100

        # bound elements hold the Variables but print compactly
        @test repr(p.elements[end]) == "VerticalCourantNumber()"
    end

    @testset "barotropic has no temperature/Courant" begin
        model = BarotropicModel(SpectralGrid(truncation = 15, nlayers = 1); feedback = Feedback(verbose = false, elements = (VerticalCourantNumber(), TemperatureRange())))
        simulation = initialize!(model)
        run!(simulation, period = Day(1))
        p = model.feedback.progress_meter
        @test all(e -> SpeedyWeather.ProgressMeter.Elements.print_element(e, p) == "", p.elements)
    end
end
