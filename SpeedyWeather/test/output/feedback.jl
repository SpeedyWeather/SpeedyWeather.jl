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
    spectral_grid = SpectralGrid(truncation = 16, nlayers = 4)
    (; render_lines) = SpeedyWeather.ProgressMeters

    @testset "no ProgressMeter.jl, no globals" begin
        @test !isdefined(SpeedyWeather, :ProgressMeter)
        @test !isdefined(SpeedyWeather, :FEEDBACK_UMAX)
    end

    @testset "default elements" begin
        feedback = Feedback(verbose = false)
        elements = ProgressElements.default_elements()
        @test feedback.elements == elements
        @test any(e -> e isa ProgressElements.SimulationSpeed, elements)
        @test any(e -> e isa ProgressElements.MaximumWindSpeed, elements)
        @test any(e -> e isa ProgressElements.VerticalCourantNumber, elements)
        @test any(e -> e isa ProgressElements.EnergyImbalance, elements)
        @test count(e -> e isa ProgressElements.NewLine, elements) == 1      # two lines
        @test any(e -> e isa ProgressElements.Bar, elements)
        # the generic elements are available under ProgressElements too
        @test ProgressElements.Bar === SpeedyWeather.ProgressMeters.Bar
        @test !hasfield(Feedback, :showspeed) && !hasfield(Feedback, :show_umax) && !hasfield(Feedback, :progress_bar_length)
        @test !hasfield(Feedback, :debug) && !hasfield(Feedback, :interval)
        @test feedback.nan_detection isa NaNDetector
    end

    @testset "default lines" begin
        model = PrimitiveWetModel(spectral_grid; feedback = Feedback(verbose = false))
        simulation = initialize!(model)
        run!(simulation, period = Day(1))
        (; variables) = simulation
        p = model.feedback.progress_meter
        @test p.counter == p.n
        lines = render_lines(p, variables, model)
        @test length(lines) == 2
        @test occursin(r"^ *\d+% [━ ]+ ETA: [\d:]+, \d{4}-\d\d-\d\d, .*/day$", lines[1])
        # (Cᵥ is wrong after the run, as `w` is still radius-scaled, but `scale` isn't anymore)
        @test occursin(r"^ {5}\d+ m/s, \[-?\d+, -?\d+\] ˚C, Cᵥ = [\d.]+, precip \d+\.\d\d mm/day, TOA [+-]\d+\.\d W/m²$", lines[2])
        # without the simulation the simulation elements print nothing, so only the first line
        @test length(render_lines(p)) == 1

        # the vertical Courant number
        @test size(variables.scratch.vertical_velocity_maximum) == (spectral_grid.nlayers,)
        w_max = reshape(variables.scratch.vertical_velocity_maximum, 1, :)
        (; w) = variables.dynamics
        Δσ = model.geometry.σ_levels_thick
        courant = ProgressElements.vertical_courant_number(w_max, w, Δσ, Float64(model.time_stepping.Δt))
        @test courant > 0 && isfinite(courant)

        # hand-computed for a single column with σ̇ = 1/s at the interfaces, Δσ = 0.25
        w = zeros(Float32, spectral_grid.grid, 4)
        w_max = zeros(Float32, 1, 4)
        w.data[:, 1:3] .= 1
        @test ProgressElements.vertical_courant_number(w_max, w, fill(0.25f0, 4), 10) ≈ 40

        # only the interface below each layer is used: σ̇ = -1/s between layers 1 and 2 counts for
        # layer 1 (Δσ = 0.5) but not for the thinner layer 2 (Δσ = 0.1) below it
        w.data .= 0
        w.data[:, 1] .= -1
        @test ProgressElements.vertical_courant_number(w_max, w, Float32[0.5, 0.1, 0.2, 0.2], 10) ≈ 20

        # precipitation and energy imbalance, hand-computed
        @test sum(model.geometry.area_weights) ≈ 1
        (; parameterizations) = variables
        parameterizations.rain_rate .= 3 / 1000 / 86400     # 3 mm/day in m/s
        parameterizations.snow_rate .= 1 / 1000 / 86400
        @test ProgressElements.print_element(ProgressElements.Precipitation(), p, variables, model) == "precip 4.00 mm/day"
        parameterizations.outgoing_shortwave .= 101.25     # incoming 1365 / 4 = 341.25 W/m²
        parameterizations.outgoing_longwave .= 200
        @test ProgressElements.print_element(ProgressElements.EnergyImbalance(), p, variables, model) == "TOA +40.0 W/m²"

        # diagnostics are not shown with Reactant unless wrapped in OnReactant
        reactant_model = (; architecture = SpeedyWeather.ReactantDevice(nothing), geometry = model.geometry, planet = model.planet)
        @test ProgressElements.print_element(ProgressElements.EnergyImbalance(), p, variables, reactant_model) == ""
        @test ProgressElements.print_element(ProgressElements.OnReactant(ProgressElements.EnergyImbalance()), p, variables, reactant_model) == "TOA +40.0 W/m²"
    end

    @testset "barotropic leaves out temperature/Courant, separator, bar" begin
        (; Percentage, Bar, SimulationTime, MaximumWindSpeed, TemperatureRange, VerticalCourantNumber, NewLine, EnergyImbalance) = ProgressElements
        elements = (Percentage(), Bar(length = 4), VerticalCourantNumber(), SimulationTime(), TemperatureRange(), MaximumWindSpeed(), NewLine(), EnergyImbalance())
        feedback = Feedback(; verbose = false, elements, separator = " | ", description = "Test")
        model = BarotropicModel(SpectralGrid(truncation = 16, nlayers = 1); feedback)
        simulation = initialize!(model)
        run!(simulation, period = Day(1))
        # single spaces around the bar, the separator between other elements, no second line
        lines = render_lines(model.feedback.progress_meter, simulation.variables, model)
        @test length(lines) == 1
        @test occursin(r"^Test +\d+% [━ ]{4} \d{4}-\d\d-\d\d \| \d+ m/s$", lines[1])
    end
end

# a user-defined element extends `ProgressElements.print_element`
struct CustomElement <: ProgressElements.AbstractSimulationElement end
ProgressElements.print_element(::CustomElement, p, vars, model) = "custom $(vars.prognostic.clock.step_counter)"

@testset "custom progress element" begin
    elements = (ProgressElements.Counter(), CustomElement())
    model = BarotropicModel(SpectralGrid(truncation = 16, nlayers = 1); feedback = Feedback(; verbose = false, elements))
    simulation = initialize!(model)
    run!(simulation, steps = 2)
    line = SpeedyWeather.ProgressMeters.render_line(model.feedback.progress_meter, simulation.variables, model)
    @test occursin(r"^2, custom \d+$", line)
end

@testset "progress! with several steps" begin
    model = BarotropicModel(SpectralGrid(truncation = 16, nlayers = 1); feedback = Feedback(verbose = false))
    simulation = initialize!(model)
    run!(simulation, steps = 2)
    p = model.feedback.progress_meter
    p.counter = 0
    SpeedyWeather.progress!(model.feedback, simulation.variables, model, 2)
    @test p.counter == 2
end

@testset "NaN detection reports the time step of the blow up" begin
    nan_detection = NaNDetector(interval = 10)
    model = BarotropicModel(SpectralGrid(truncation = 16, nlayers = 1); feedback = Feedback(; verbose = false, nan_detection))
    simulation = initialize!(model)
    (; variables) = simulation
    @test length(variables.dynamics.nan_detection.nan_record) == 10
    run!(simulation, steps = 2)
    initialize!(model.feedback, variables, model)       # reset the progress meter and the NaN detection
    model.feedback.progress_meter.n = 100
    vor = SpeedyWeather.get_prognostic_step(variables.prognostic.vorticity, model.time_stepping, nan_detection)

    # recording a time step does not allocate (measured inside a function, see above)
    record(d, vars, model) = @allocated SpeedyWeather.detect_nans!(d, vars, model, 1)
    record(nan_detection, variables, model)
    @test record(nan_detection, variables, model) == 0
    initialize!(nan_detection, variables, model)

    # blow up at time step 4, only checked (and reported) at step 10, the interval
    for step in 1:10
        step == 4 && (vor.data .= NaN)
        if step < 10
            SpeedyWeather.progress!(model.feedback, variables, model)
            @test !SpeedyWeather.nans_detected(model.feedback)
        else
            @test_logs (:warn, r"time step 4 ") SpeedyWeather.progress!(model.feedback, variables, model)
        end
    end
    @test SpeedyWeather.nans_detected(model.feedback)
    @test nan_detection.step == 4
    @test nan_detection.time == variables.prognostic.clock.time

    # without NaN detection
    model = BarotropicModel(SpectralGrid(truncation = 16, nlayers = 1); feedback = Feedback(verbose = false, nan_detection = nothing))
    simulation = initialize!(model)
    @test !hasproperty(simulation.variables.dynamics, :nan_detection)
    run!(simulation, steps = 2)
    @test !SpeedyWeather.nans_detected(model.feedback)
    @test !SpeedyWeather.nans_detected(nothing)
end
