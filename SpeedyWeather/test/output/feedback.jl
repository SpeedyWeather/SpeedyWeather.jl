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
        @test feedback.elements == ProgressElements.default_elements()
        @test any(e -> e isa ProgressElements.SimulationSpeed, feedback.elements)
        @test !any(e -> e isa ProgressElements.VerticalCourantNumber, feedback.elements)
        # ProgressMeter's generic elements are available under ProgressElements too
        @test ProgressElements.Bar === SpeedyWeather.ProgressMeter.Elements.Bar
        @test !hasfield(Feedback, :showspeed) && !hasfield(Feedback, :show_umax)
    end

    @testset "elements on redraw, VerticalCourantNumber" begin
        feedback = Feedback(verbose = false, elements = (ProgressElements.default_elements()..., ProgressElements.VerticalCourantNumber()))
        model = PrimitiveDryModel(spectral_grid; feedback)
        simulation = initialize!(model)
        run!(simulation, period = Day(1))
        p = model.feedback.progress_meter
        line = join(map(e -> ProgressElements.print_element(e, p), p.elements))
        # the simulation elements are separated by ", ", also the appended one
        @test occursin(r"ETA: .*, \d{4}-\d\d-\d\d, .*/day, +\d+ m/s, \[ *-?\d+, +-?\d+\] ˚C, Cᵥ = \d", line)
        @test !occursin("(", line)
        # the vertical scratch vector is allocated for the element, the redraw doesn't allocate
        element = p.elements[end]
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
        @test repr(p.elements[end]) == "VerticalCourantNumber()"
    end

    @testset "barotropic leaves out temperature/Courant, separators" begin
        (; Description, Percentage, SimulationTime, MaximumWindSpeed, TemperatureRange, VerticalCourantNumber) = ProgressElements
        elements = (Description(), Percentage(), VerticalCourantNumber(), SimulationTime(), TemperatureRange(), MaximumWindSpeed())
        model = BarotropicModel(SpectralGrid(truncation = 15, nlayers = 1); feedback = Feedback(; verbose = false, elements, separator = " | "))
        simulation = initialize!(model)
        run!(simulation, period = Day(1))
        bound = model.feedback.progress_meter.elements
        @test bound[1] isa Description && bound[2] isa Percentage
        @test bound[3] == " | " && bound[4] isa SimulationTime && bound[5] == " | " && bound[6] isa MaximumWindSpeed
        @test length(bound) == 6

        # no separator after only the description or strings
        bound = ProgressElements.bind_elements((Description(), "[", SimulationTime(), MaximumWindSpeed()), simulation.variables, model)
        @test bound[2] == "[" && bound[3] isa SimulationTime && bound[4] == ", " && bound[5] isa MaximumWindSpeed
    end
end

# a user-defined element extends `ProgressElements.print_element`
struct CustomElement <: ProgressElements.AbstractProgressElement end
ProgressElements.print_element(::CustomElement, p) = " custom"

@testset "custom progress element" begin
    model = BarotropicModel(SpectralGrid(truncation = 15, nlayers = 1); feedback = Feedback(verbose = false, elements = (ProgressElements.Counter(), CustomElement())))
    simulation = initialize!(model)
    run!(simulation, steps = 2)
    @test ProgressElements.print_element(model.feedback.progress_meter.elements[2], model.feedback.progress_meter) == " custom"
end
