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
    render_line = SpeedyWeather.ProgressMeters.render_line

    @testset "no ProgressMeter.jl, no globals" begin
        @test !isdefined(SpeedyWeather, :ProgressMeter)
        @test !isdefined(SpeedyWeather, :FEEDBACK_UMAX)
    end

    @testset "default elements" begin
        feedback = Feedback(verbose = false)
        @test isnothing(feedback.elements)      # default_elements(architecture) in initialize!
        elements = ProgressElements.default_elements()
        @test elements == ProgressElements.default_elements(SpeedyWeather.CPU())
        @test any(e -> e isa ProgressElements.SimulationSpeed, elements)
        @test any(e -> e isa ProgressElements.MaximumWindSpeed, elements)
        @test any(e -> e isa ProgressElements.VerticalCourantNumber, elements)
        @test any(e -> e isa ProgressElements.EnergyImbalance, elements)
        @test count(e -> e isa ProgressElements.NewLine, elements) == 1      # two lines
        @test any(e -> e isa ProgressElements.Bar, elements)
        # the diagnostics are opt-in for Reactant
        reactant = SpeedyWeather.ReactantDevice(nothing)   # for dispatch only, no Reactant needed
        elements = ProgressElements.default_elements(reactant)
        @test !any(e -> e isa ProgressElements.AbstractSimulationElement && !(e isa Union{ProgressElements.SimulationTime, ProgressElements.SimulationSpeed}), elements)
        @test !any(e -> e isa ProgressElements.NewLine, elements)
        @test any(e -> e isa ProgressElements.SimulationTime, elements)
        # the generic elements are available under ProgressElements too
        @test ProgressElements.Bar === SpeedyWeather.ProgressMeters.Bar
        @test !hasfield(Feedback, :showspeed) && !hasfield(Feedback, :show_umax) && !hasfield(Feedback, :progress_bar_length)
    end

    @testset "default lines" begin
        model = PrimitiveWetModel(spectral_grid; feedback = Feedback(verbose = false))
        simulation = initialize!(model)
        run!(simulation, period = Day(1))
        p = model.feedback.progress_meter
        @test p.counter == p.n
        lines = SpeedyWeather.ProgressMeters.render_lines(p)
        @test length(lines) == 2
        @test occursin(r"^ *\d+% [━ ]+ ETA: [\d:]+, \d{4}-\d\d-\d\d, .*/day$", lines[1])
        # (Cᵥ is wrong after the run, as `w` is still radius-scaled, but `scale` isn't anymore)
        @test occursin(r"^ {5}\d+ m/s, \[-?\d+, -?\d+\] ˚C, Cᵥ = [\d.]+, precip \d+\.\d\d mm/day, TOA [+-]\d+\.\d W/m²$", lines[2])

        # the vertical Courant number
        element = p.elements[findfirst(e -> e isa ProgressElements.VerticalCourantNumber, p.elements)]
        # the vertical scratch vector is allocated for the element, the redraw doesn't allocate
        @test element.w_max isa AbstractMatrix && size(element.w_max) == (1, spectral_grid.nlayers)
        (; w) = simulation.variables.dynamics
        Δt = Float64(model.time_stepping.Δt)
        courant = ProgressElements.vertical_courant_number!(element.w_max, element.w_max_cpu, w, element.Δσ, Δt)
        @test courant > 0 && isfinite(courant)
        # measure inside a function, a call from the testset scope is a dynamic dispatch that allocates
        allocations(e, w, Δt) = @allocated ProgressElements.vertical_courant_number!(e.w_max, e.w_max_cpu, w, e.Δσ, Δt)
        allocations(element, w, Δt)
        @test allocations(element, w, Δt) == 0

        # hand-computed for a single column with σ̇ = 1/s at the interfaces, Δσ = 0.25
        w = zeros(Float32, spectral_grid.grid, 4)
        w_max, w_max_cpu = zeros(Float32, 1, 4), zeros(Float32, 4)
        w.data[:, 1:3] .= 1
        @test ProgressElements.vertical_courant_number!(w_max, w_max_cpu, w, fill(0.25f0, 4), 10) ≈ 40

        # only the interface below each layer is used: σ̇ = -1/s between layers 1 and 2 counts for
        # layer 1 (Δσ = 0.5) but not for the thinner layer 2 (Δσ = 0.1) below it
        w.data .= 0
        w.data[:, 1] .= -1
        @test ProgressElements.vertical_courant_number!(w_max, w_max_cpu, w, Float32[0.5, 0.1, 0.2, 0.2], 10) ≈ 20

        # bound elements hold the Variables but print compactly
        @test repr(element) == "VerticalCourantNumber()"

        # precipitation and energy imbalance, hand-computed
        @test sum(ProgressElements.global_mean_weights(model)) ≈ 1
        (; parameterizations) = simulation.variables
        parameterizations.rain_rate .= 3 / 1000 / 86400     # 3 mm/day in m/s
        parameterizations.snow_rate .= 1 / 1000 / 86400
        precipitation = p.elements[findfirst(e -> e isa ProgressElements.Precipitation, p.elements)]
        @test ProgressElements.print_element(precipitation, p) == "precip 4.00 mm/day"
        parameterizations.outgoing_shortwave .= 101.25     # incoming 1365 / 4 = 341.25 W/m²
        parameterizations.outgoing_longwave .= 200
        energy = p.elements[findfirst(e -> e isa ProgressElements.EnergyImbalance, p.elements)]
        @test ProgressElements.print_element(energy, p) == "TOA +40.0 W/m²"
    end

    @testset "barotropic leaves out temperature/Courant, separator, bar" begin
        (; Percentage, Bar, SimulationTime, MaximumWindSpeed, TemperatureRange, VerticalCourantNumber) = ProgressElements
        elements = (Percentage(), Bar(length = 4), VerticalCourantNumber(), SimulationTime(), TemperatureRange(), MaximumWindSpeed())
        feedback = Feedback(; verbose = false, elements, separator = " | ", description = "Test")
        model = BarotropicModel(SpectralGrid(truncation = 16, nlayers = 1); feedback)
        simulation = initialize!(model)
        run!(simulation, period = Day(1))
        p = model.feedback.progress_meter
        @test p.elements[1] isa Percentage && p.elements[2] isa Bar
        @test p.elements[3] isa SimulationTime && p.elements[4] isa MaximumWindSpeed
        @test length(p.elements) == 4
        # single spaces around the bar, the separator between other elements
        @test occursin(r"^Test +\d+% [━ ]{4} \d{4}-\d\d-\d\d \| +\d+ m/s$", render_line(p))
    end
end

# a user-defined element extends `ProgressElements.print_element`
struct CustomElement <: ProgressElements.AbstractProgressElement end
ProgressElements.print_element(::CustomElement, p) = "custom"

@testset "custom progress element" begin
    elements = (ProgressElements.Counter(), CustomElement())
    model = BarotropicModel(SpectralGrid(truncation = 16, nlayers = 1); feedback = Feedback(; verbose = false, elements))
    simulation = initialize!(model)
    run!(simulation, steps = 2)
    @test SpeedyWeather.ProgressMeters.render_line(model.feedback.progress_meter) == "2, custom"
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
    feedback = Feedback(verbose = false, interval = 10)
    model = BarotropicModel(SpectralGrid(truncation = 16, nlayers = 1); feedback)
    simulation = initialize!(model)
    (; variables) = simulation
    run!(simulation, steps = 2)
    initialize!(feedback, variables, model)     # reset the progress meter and the NaN records
    feedback.progress_meter.n = 100
    vor = SpeedyWeather.get_prognostic_step(variables.prognostic.vorticity, model.time_stepping, feedback)

    # recording a time step does not allocate (measured inside a function, see above)
    allocations(feedback, variables, model) = @allocated SpeedyWeather.record_nans!(feedback, variables, model)
    allocations(feedback, variables, model)
    @test allocations(feedback, variables, model) == 0
    feedback.nan_record_count = 0

    # blow up at time step 4, only checked (and reported) at step 10, the interval
    for step in 1:10
        step == 4 && (vor.data .= NaN)
        if step < 10
            SpeedyWeather.progress!(feedback, variables, model)
            @test !feedback.nans_detected
        else
            @test_logs (:warn, r"time step 4 ") SpeedyWeather.progress!(feedback, variables, model)
        end
    end
    @test feedback.nans_detected
    @test feedback.nans_detected_step == 4
    @test feedback.nans_detected_time == variables.prognostic.clock.time
end
