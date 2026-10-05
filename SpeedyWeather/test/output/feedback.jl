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
        status = SpeedyWeather.ProgressMeter.ProgressStatus(time(), 1.0, false)
        line = join(map(e -> SpeedyWeather.ProgressMeter.print_element(e, p, status), p.elements[5:end]))
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
    end

    @testset "barotropic has no temperature/Courant" begin
        model = BarotropicModel(SpectralGrid(truncation = 15, nlayers = 1); feedback = Feedback(verbose = false, elements = (VerticalCourantNumber(), TemperatureRange())))
        simulation = initialize!(model)
        run!(simulation, period = Day(1))
        p = model.feedback.progress_meter
        status = SpeedyWeather.ProgressMeter.ProgressStatus(time(), 1.0, false)
        @test all(e -> SpeedyWeather.ProgressMeter.print_element(e, p, status) == "", p.elements)
    end
end
