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
  
    # Output writers can't interpolate on Reactant arrays outside of compiled code, so they
    # interpolate on the CPU for a Reactant model, moving the model's fields there first.
    @testset "Output writers interpolate on the CPU for Reactant" begin
        arch = ReactantDevice()
        cpu = SpeedyWeather.CPU()
        @test SpeedyWeather.output_architecture(arch) isa SpeedyWeather.CPU

        spectral_grid = SpectralGrid(; architecture = arch, truncation = TRUNCATION)
        spectral_grid_cpu = SpectralGrid(; truncation = TRUNCATION)
        output = NetCDFOutput(spectral_grid)
        output_cpu = NetCDFOutput(spectral_grid_cpu)

        @test SpeedyWeather.ismatching(cpu, output.field2D.data)
        @test SpeedyWeather.ismatching(cpu, output.interpolator.locator.ij_as)
        @test output.host2D === output.field2D          # written from directly, no extra copy

        # a field on the Reactant device interpolates to the same as on the CPU
        field_cpu = rand(Float32, spectral_grid_cpu.grid)
        field = on_architecture(arch, field_cpu)
        SpeedyWeather.interpolate_output!(output, output.field2D, field)
        SpeedyWeather.interpolate_output!(output_cpu, output_cpu.field2D, field_cpu)
        @test output.field2D.data == output_cpu.field2D.data
    end

    # the height of the bottom layer is computed in a kernel, in a broadcast Reactant traces the
    # layer index for which `pressure(k::Integer, ...)` doesn't apply
    @testset "10m wind output on Reactant" begin
        variable = SpeedyWeather.ZonalVelocity10mOutput()
        simulations = map((SpeedyWeather.CPU(), ReactantDevice())) do architecture
            spectral_grid = SpectralGrid(; architecture, nlayers = 8, truncation = TRUNCATION)
            spectral_transform = MatrixSpectralTransform(spectral_grid)
            path = mktempdir(pwd(), prefix = "tmp_reactant_nc_")
            output = NetCDFOutput(spectral_grid, PrimitiveDryModel; path)
            add!(output, variable)
            model = PrimitiveDryModel(spectral_grid; spectral_transform, feedback = nothing, output)
            initialize!(model)
        end
        sync_variables!(simulations...)     # same state on CPU and Reactant

        z_bottom = map(simulations) do simulation
            initialize!(simulation; steps = 1, output = true)   # writes the initial conditions incl. u10
            SpeedyWeather.output!(simulation.model.output, variable, simulation)
            z_bottom = on_architecture(SpeedyWeather.CPU(), simulation.variables.scratch.grid.a_2D.data)
            finalize!(simulation)
            z_bottom
        end
        @test all(isfinite, z_bottom[2])
        @test z_bottom[2] ≈ z_bottom[1]
    end
end
