# Smoke-test that `output!` actually writes data when the model is run via
# Reactant. The Reactant extension routes output around the compiled
# `time_step!` (output is skipped during tracing and
# `output!` is invoked explicitly from `time_stepping!` between compiled steps), so this
# test guards against regressions in that path for both NetCDF and Zarr.
#
# Notes:
#  - Uses `period=Hour(20)` rather than `steps=20` because the Reactant clock
#    converts its `period` field to integer-seconds; `steps * Δt` can land on a
#    fractional second when `interval` forces Δt adjustment.
#  - Picks `interval = Hour(1)`, a multiple of the default Barotropic Δt so
#    `adjust_with_output` doesn't shift Δt away from a whole second.

using NCDatasets
using Zarr

@testset "Reactant output writers" begin
    nsteps = 20
    interval = Hour(1)
    period = nsteps * interval

    arch = ReactantDevice()
    SG = SpectralGrid(architecture = arch, truncation = TRUNCATION, nlayers = 1)
    M = MatrixSpectralTransform(SG)

    @testset "NetCDFOutput" begin
        tmp = mktempdir(pwd(), prefix = "tmp_reactant_nc_")
        output = NetCDFOutput(SG, BarotropicModel; path = tmp, interval)
        model = BarotropicModel(SG; spectral_transform = M, feedback = nothing, output)
        simulation = initialize!(model)
        initialize!(simulation; period, output = true)

        SpeedyWeather.time_stepping!(simulation)
        SpeedyWeather.finalize!(simulation)

        # IC + one snapshot per output interval
        @test output.output_counter == nsteps + 1

        ds = NCDataset(joinpath(output.run_path, output.filename))
        @test haskey(ds, "vor")
        @test haskey(ds, "u")
        @test haskey(ds, "v")
        @test haskey(ds, "time")

        @test length(ds["time"][:]) == nsteps + 1
        nx, ny, nz, nt = size(ds["vor"])
        @test nt == nsteps + 1
        @test nz == SG.nlayers

        vor = ds["vor"][:, :, :, :]
        # Reactant simulation should produce finite, non-zero output that
        # actually changes from step to step.
        @test all(isfinite, vor)
        @test any(!iszero, vor)
        @test vor[:, :, :, 1] != vor[:, :, :, end]

        close(ds)
    end

    @testset "ZarrOutput" begin
        tmp = mktempdir(pwd(), prefix = "tmp_reactant_zarr_")
        output = ZarrOutput(SG, BarotropicModel; path = tmp, interval)
        model = BarotropicModel(SG; spectral_transform = M, feedback = nothing, output)
        simulation = initialize!(model)
        initialize!(simulation; period, output = true)

        SpeedyWeather.time_stepping!(simulation)
        SpeedyWeather.finalize!(simulation)

        @test output.output_counter == nsteps + 1

        zg = Zarr.zopen(joinpath(output.run_path, output.filename))
        @test haskey(zg.arrays, "vor")
        @test haskey(zg.arrays, "u")
        @test haskey(zg.arrays, "v")
        @test haskey(zg.arrays, "time")

        @test length(zg["time"][:]) == nsteps + 1
        nx, ny, nz, nt = size(zg["vor"])
        @test nt == nsteps + 1
        @test nz == SG.nlayers

        vor = zg["vor"][:, :, :, :]
        @test all(isfinite, vor)
        @test any(!iszero, vor)
        @test vor[:, :, :, 1] != vor[:, :, :, end]
    end

    # Output writers interpolate on the Reactant device (compiled), only the interpolated fields
    # are copied to the CPU for writing
    @testset "Output writers interpolate on the Reactant device" begin
        arch = ReactantDevice()
        cpu = SpeedyWeather.CPU()

        spectral_grid = SpectralGrid(; architecture = arch, truncation = TRUNCATION)
        spectral_grid_cpu = SpectralGrid(; truncation = TRUNCATION)
        output = NetCDFOutput(spectral_grid)
        output_cpu = NetCDFOutput(spectral_grid_cpu)

        @test SpeedyWeather.ismatching(arch, output.field2D.data)
        @test SpeedyWeather.ismatching(arch, output.interpolator.locator.ij_as)
        @test SpeedyWeather.ismatching(cpu, output.host2D.data)     # written from the host

        # a field on the Reactant device interpolates to the same as on the CPU
        field_cpu = rand(Float32, spectral_grid_cpu.grid)
        field = on_architecture(arch, field_cpu)
        SpeedyWeather.interpolate_output!(output, output.field2D, field)
        SpeedyWeather.interpolate_output!(output_cpu, output_cpu.field2D, field_cpu)
        @test Array(output.field2D.data) ≈ output_cpu.field2D.data
    end

    # With Reactant every output step is one compiled program for all output variables
    # (`interpolate_output!` of every variable, everything but the I/O), compiled once with the
    # initial conditions. Same files as on the CPU, and an output step costs about as much as on the CPU
    @testset "Compiled output of all variables" begin
        variables = (
            OutputVariables.Divergence(), OutputVariables.SurfaceTemperature(),
            # the height of the bottom layer is computed in a kernel, in a broadcast Reactant traces the
            # layer index for which `pressure(k::Integer, ...)` doesn't apply
            OutputVariables.ZonalVelocity10m(),
            OutputVariables.PrecipitableWater(),
            OutputVariables.LandSeaMask(),      # no time dimension, only written with the initial conditions
        )
        simulations = map((SpeedyWeather.CPU(), ReactantDevice())) do architecture
            output = NetCDFOutput(
                SpectralGrid(; architecture, nlayers = 8, truncation = TRUNCATION), PrimitiveWetModel;
                path = mktempdir(pwd(), prefix = "tmp_reactant_nc_")
            )
            add!(output, variables...)
            initialize!(create_model(PrimitiveWetModel, architecture; output))
        end
        sync_variables!(simulations...)     # same state on CPU and Reactant

        # write the initial conditions (compiles the output program on Reactant) and time an output step
        output_step_times = map(simulations) do simulation
            initialize!(simulation; steps = 1, output = true)
            (; output) = simulation.model
            minimum(@elapsed(SpeedyWeather.output!(output, output.variables, simulation)) for _ in 1:3)
        end

        output_cpu, output = map(simulation -> simulation.model.output, simulations)
        compiled = output.compiled
        @test compiled isa SpeedyWeatherReactantExt.CompiledOutput
        @test length(compiled.variables) == length(output.variables)
        @test all(field -> SpeedyWeather.ismatching(ReactantDevice(), field.data), compiled.fields)

        SpeedyWeather.output!(output, output.variables, simulations[2])
        @test output.compiled === compiled  # reused for every output step, not recompiled

        # performance: once compiled, an output step costs about as much as on the CPU (the writing,
        # i.e. compression, dominates), not 10-100x as when interpolating eagerly, variable by variable
        @info "Output step CPU vs Reactant" output_step_times
        @test output_step_times[2] < 5 * output_step_times[1]

        foreach(finalize!, simulations)
        datasets = map(output -> NCDataset(joinpath(output.run_path, output.filename)), (output_cpu, output))
        for variable in values(output.variables)
            a, b = map(dataset -> Array(dataset[variable.name].var), datasets)
            @testset "$(variable.name)" begin
                @test size(a) == size(b)
                @test isnan.(a) == isnan.(b)
                finite = .!isnan.(a) .& .!isnan.(b)
                @test isapprox(a[finite], b[finite], rtol = 1.0e-2)
            end
        end
        foreach(close, datasets)
    end
end
