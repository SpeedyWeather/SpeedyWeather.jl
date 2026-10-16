using Zarr
import SpeedyWeather.RingGrids
import SpeedyWeather.NCDatasets: NCDataset

@testset "GPU output: interpolation on the model's architecture" begin
    arch = SpeedyWeather.GPU()
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 4, architecture = arch)
    spectral_grid_cpu = SpectralGrid(truncation = 31, nlayers = 4)

    for Writer in (NetCDFOutput, ZarrOutput)
        output = Writer(spectral_grid, PrimitiveWet)
        output_cpu = Writer(spectral_grid_cpu, PrimitiveWet)

        # interpolator and scratch fields live on the GPU, host copies on the CPU
        @test SpeedyWeather.ismatching(arch, output.interpolator.locator.ij_as)
        @test SpeedyWeather.ismatching(arch, output.field3D.data)
        @test SpeedyWeather.ismatching(arch, output.land_fraction.data)
        @test SpeedyWeather.ismatching(SpeedyWeather.CPU(), output.host3D.data)
        @test output.host3D !== output.field3D
        @test size(output.host3D) == size(output.field3D)

        # same interpolation on GPU and CPU, 2D and 3D
        for (field, field_cpu) in ((output.field2D, output_cpu.field2D), (output.field3D, output_cpu.field3D))
            src_cpu = rand(Float32, spectral_grid_cpu.grid, size(field)[2:end]...)
            src = on_architecture(arch, src_cpu)
            SpeedyWeather.interpolate_output!(output, field, src)
            SpeedyWeather.interpolate_output!(output_cpu, field_cpu, src_cpu)
            @test on_architecture(SpeedyWeather.CPU(), field.data) ≈ field_cpu.data
        end
    end

    # HEALPixOutput: with interpolation onto a coarser grid, and without on the model's grid
    output = HEALPixOutput(spectral_grid, PrimitiveWet, nside = 8)
    @test SpeedyWeather.ismatching(arch, output.interpolator.locator.ij_as)
    @test SpeedyWeather.ismatching(arch, output.field2D.data)

    spectral_grid_healpix = SpectralGrid(truncation = 31, nlayers = 4, Grid = HEALPixGrid, architecture = arch)
    output = HEALPixOutput(spectral_grid_healpix, PrimitiveWet)
    @test isnothing(output.interpolator)
    @test SpeedyWeather.ismatching(arch, output.field2D.data)
end

@testset "GPU output: vertical interpolation onto pressure layers" begin
    arch = SpeedyWeather.GPU()
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 8, architecture = arch)
    spectral_grid_cpu = SpectralGrid(truncation = 31, nlayers = 8)
    coordinates = SigmaCoordinates(spectral_grid)
    coordinates_cpu = SigmaCoordinates(spectral_grid_cpu)
    NF = spectral_grid.NF

    pₛ_cpu = fill!(zeros(NF, spectral_grid_cpu.grid), NF(1000.0e2))
    in_field_cpu = 200 .+ rand(NF, spectral_grid_cpu.grid, spectral_grid.nlayers)
    p_cpu = NF[10.0e2, 500.0e2, 990.0e2, 1010.0e2]  # above, inside, below the model layers, below ground
    out_cpu = zeros(NF, spectral_grid_cpu.grid, length(p_cpu))

    pₛ, in_field, p = on_architecture(arch, pₛ_cpu), on_architecture(arch, in_field_cpu), on_architecture(arch, p_cpu)
    out = on_architecture(arch, out_cpu)

    # parameters in the model's number format, as set by sync_extrapolations! at initialize!
    adiabatic = SpeedyWeather.DryAdiabaticExtrapolation(κ = NF(2 / 7))
    for extrapolation in (
            SpeedyWeather.ConstantExtrapolation(),
            adiabatic,
            SpeedyWeather.SubsurfaceMask(above_surface = adiabatic, missing_value = NF(NaN)),
        )
        SpeedyWeather.interpolate_pressure_layers!(out, in_field, pₛ, p, coordinates, SpeedyWeather.LinearInLogPressure(), extrapolation)
        SpeedyWeather.interpolate_pressure_layers!(out_cpu, in_field_cpu, pₛ_cpu, p_cpu, coordinates_cpu, SpeedyWeather.LinearInLogPressure(), extrapolation)
        @test on_architecture(SpeedyWeather.CPU(), out.data) ≈ out_cpu.data nans = true
    end
end

@testset "GPU output: all variables, all writers" begin
    arch = SpeedyWeather.GPU()
    path = mktempdir(pwd(), prefix = "tmp_gpu_output_")
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 8, architecture = arch)

    writers = (
        NetCDFOutput(spectral_grid, PrimitiveWet; path),
        NetCDFOutput(spectral_grid, PrimitiveWet; path, layers = PressureLayers(spectral_grid)),
        ZarrOutput(spectral_grid, PrimitiveWet; path),
        HEALPixOutput(spectral_grid, PrimitiveWet; path),
    )

    for output in writers
        add!(output, SpeedyWeather.AllOutputVariables()...)
        model = PrimitiveWetModel(spectral_grid; output)
        simulation = initialize!(model)

        # write every time step so that the last snapshot is the final state
        set!(output, model, interval = model.time_stepping.Δt_millisec)
        run!(simulation, steps = 3, output = true)
        @test simulation.model.feedback.nans_detected == false

        # the last written temperature agrees with the final state interpolated on the CPU
        temperature = SpeedyWeather.get_prognostic_step(simulation.variables.grid.temperature, model.time_stepping, output)
        temperature = on_architecture(SpeedyWeather.CPU(), temperature)     # the step written out
        output.layers isa PressureLayers && continue    # vertical interpolation on top, not compared
        output isa HEALPixOutput && continue            # different output grid, covered above
        expected = RingGrids.interpolate(on_architecture(SpeedyWeather.CPU(), output.field3D.grid), temperature)
        grid = output.host3D.grid
        expected = reshape(expected.data, length(RingGrids.get_lond(grid)), length(RingGrids.get_latd(grid)), :)   # (lon, lat, layer)

        file = joinpath(output.run_path, output.filename)
        written = output isa NetCDFOutput ?
            NCDataset(file) do ds
                ds["temp"].var[:, :, :, end]
            end : zopen(file)["temp"][:, :, :, end]

        # in ˚C and rounded to keepbits=10, i.e. a relative error of ~1e-3 in K
        @test written .+ 273.15f0 ≈ expected rtol = 2.0e-3
    end
end
