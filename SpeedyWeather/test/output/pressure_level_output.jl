using NCDatasets, Dates

@testset "Output on pressure levels" begin

    @testset "PressureLevels construction" begin
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8)

        levels = SpeedyWeather.PressureLevels()
        @test SpeedyWeather.get_nlayers(levels, spectral_grid) == length(levels.pressure)
        @test SpeedyWeather.vertical_dimension(levels) == "pressure"
        @test SpeedyWeather.vertical_dimension(SpeedyWeather.ModelLevels()) == "layer"

        # model levels is the default and sizes the 3D scratch field with nlayers
        output = NetCDFOutput(spectral_grid, PrimitiveWet)
        @test output.levels isa SpeedyWeather.ModelLevels
        @test size(output.field3D, 2) == spectral_grid.nlayers

        # pressure levels size it with the number of pressure levels instead
        p = [850, 500, 200] .* 100.0
        output = NetCDFOutput(spectral_grid, PrimitiveWet, levels = SpeedyWeather.PressureLevels(p))
        @test size(output.field3D, 2) == length(p)
        @test size(output.field3Dland, 2) == SpeedyWeather.DEFAULT_NLAYERS_SOIL   # unaffected

        # non-monotonic pressure levels are rejected at initialize!
        model = PrimitiveWetModel(spectral_grid)
        levels = SpeedyWeather.PressureLevels([500, 200, 850] .* 100.0)
        @test_throws AssertionError SpeedyWeather.initialize!(levels, model)
        levels = SpeedyWeather.PressureLevels([-500, 200] .* 100.0)
        @test_throws AssertionError SpeedyWeather.initialize!(levels, model)
    end

    @testset "Extrapolation defaults" begin
        # temperature descends dry-adiabatically below the lowest model level,
        # everything else is constant
        @test SpeedyWeather.output_extrapolation(SpeedyWeather.TemperatureOutput()) isa
            SpeedyWeather.DryAdiabaticExtrapolation
        for var in (SpeedyWeather.VorticityOutput(), SpeedyWeather.ZonalVelocityOutput(), SpeedyWeather.HumidityOutput())
            @test SpeedyWeather.output_extrapolation(var) isa SpeedyWeather.ConstantExtrapolation
        end

        # κ is taken from the model's atmosphere at initialize!
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8)
        output = NetCDFOutput(spectral_grid, PrimitiveWet, levels = SpeedyWeather.PressureLevels())
        model = PrimitiveWetModel(spectral_grid; output)
        SpeedyWeather.sync_extrapolations!(output, model)
        @test output.variables[:temp].extrapolation.κ == model.atmosphere.κ
        @test output.variables[:temp].extrapolation.κ isa spectral_grid.NF

        # a mask keeps its inner extrapolation synced too
        output.variables[:temp].extrapolation =
            SpeedyWeather.SubsurfaceMask(above_surface = SpeedyWeather.DryAdiabaticExtrapolation())
        SpeedyWeather.sync_extrapolations!(output, model)
        @test output.variables[:temp].extrapolation.above_surface.κ == model.atmosphere.κ
    end

    @testset "NetCDFOutput on pressure levels" begin
        # full grid model so that the horizontal interpolation onto the (full) output grid
        # is an exact copy, which lets us compare the file against the interpolation directly
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8, Grid = FullGaussianGrid)
        NF = spectral_grid.NF
        p = [850, 500, 200] .* 100.0

        output = NetCDFOutput(
            spectral_grid, PrimitiveWet,
            path = mktempdir(), write_restart = false, interval = Hour(6),
            levels = SpeedyWeather.PressureLevels(p),
        )
        # no bitrounding and soil output so that the comparison below is exact
        add!(output, SpeedyWeather.TemperatureOutput(keepbits = 23), SpeedyWeather.SoilTemperatureOutput())

        model = PrimitiveWetModel(spectral_grid; output)
        simulation = initialize!(model)
        run!(simulation, period = Hour(6), output = true)
        @test model.feedback.nans_detected == false

        NCDataset(SpeedyWeather.get_output_path(model.output)) do ds
            # the pressure dimension replaces the sigma layer dimension
            @test haskey(ds.dim, "pressure")
            @test !haskey(ds.dim, "layer")
            @test ds.dim["pressure"] == length(p)
            @test ds["pressure"][:] ≈ p ./ 100                  # written in hPa
            @test ds["pressure"].attrib["units"] == "hPa"

            # 3D atmospheric variables are on it, soil variables and 2D variables are not
            for name in ("temp", "u", "v", "vor", "humid")
                @test NCDatasets.dimnames(ds[name]) == ("lon", "lat", "pressure", "time")
            end
            @test NCDatasets.dimnames(ds["st"]) == ("lon", "lat", "soil_layer", "time")
            @test NCDatasets.dimnames(ds["mslp"]) == ("lon", "lat", "time")
            @test haskey(ds.dim, "soil_layer")

            # values match interpolating the final state directly
            temp = SpeedyWeather.get_prognostic_step(
                simulation.variables.grid.temperature, model.time_stepping, model.output,
            )
            pₛ = simulation.variables.parameterizations.surface_pressure
            reference = zeros(NF, spectral_grid.grid, length(p))
            SpeedyWeather.interpolate_pressure_levels!(
                reference, temp, pₛ, NF.(p), model.geometry.vertical_coordinates,
                model.output.levels.interpolation,
                model.output.variables[:temp].extrapolation,
            )

            nlon, nlat = length(ds["lon"]), length(ds["lat"])
            from_file = ds["temp"][:, :, :, end] .+ 273.15      # ˚C back to K
            reference_lonlat = reshape(Array(reference.data), nlon, nlat, length(p))
            @test from_file ≈ reference_lonlat rtol = 1.0e-5

            # and are physically sensible: colder aloft, no missing values
            @test all(isfinite, from_file)
            @test sum(from_file[:, :, 1]) > sum(from_file[:, :, 2]) > sum(from_file[:, :, 3])
        end
    end

    @testset "Model levels unchanged by default" begin
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8)
        output = NetCDFOutput(
            spectral_grid, PrimitiveWet,
            path = mktempdir(), write_restart = false, interval = Hour(6),
        )
        model = PrimitiveWetModel(spectral_grid; output)
        simulation = initialize!(model)
        run!(simulation, period = Hour(6), output = true)

        NCDataset(SpeedyWeather.get_output_path(model.output)) do ds
            @test haskey(ds.dim, "layer")
            @test !haskey(ds.dim, "pressure")
            @test ds["layer"][:] ≈ model.geometry.σ_levels_full
            @test ds["layer"].attrib["long_name"] == "sigma layer"
            @test NCDatasets.dimnames(ds["temp"]) == ("lon", "lat", "layer", "time")
        end
    end

end
