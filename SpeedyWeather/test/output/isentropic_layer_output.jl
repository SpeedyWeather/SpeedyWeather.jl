using NCDatasets, Dates
import SpeedyWeather: interpolate_isentropic_layers!, pressure, dry_potential_temperature
import SpeedyWeather: ConstantExtrapolation, DryAdiabaticExtrapolation, SubsurfaceMask

# temperature [K] on the model layers from a profile θ(k, p) of dry potential temperature,
# T = θ (p/p₀)^κ, the same in every column for constant pₛ
function isentropic_test_temperature(spectral_grid, coordinates, pₛ, θ_profile::Function, κ, p₀)
    temperature = zeros(spectral_grid.NF, spectral_grid.grid, spectral_grid.nlayers)
    for k in 1:spectral_grid.nlayers, ij in eachgridpoint(temperature)
        p = pressure(k, pₛ[ij], coordinates)
        temperature[ij, k] = θ_profile(k, p) * (p / p₀)^κ
    end
    return temperature
end

# a field that is an exact function f of pressure on the model layers
function isentropic_test_field(spectral_grid, coordinates, pₛ, f::Function)
    field = zeros(spectral_grid.NF, spectral_grid.grid, spectral_grid.nlayers)
    for k in 1:spectral_grid.nlayers, ij in eachgridpoint(field)
        field[ij, k] = f(pressure(k, pₛ[ij], coordinates))
    end
    return field
end

@testset "Vertical interpolation onto isentropic layers" begin
    κ, p₀ = 2 / 7, 1.0e5

    @testset "θ linear in ln(p) with $VerticalCoordinates" for VerticalCoordinates in
        (SigmaCoordinates, SigmaPressureCoordinates)

        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8, NF = Float64)
        coordinates = VerticalCoordinates(spectral_grid)
        nlayers = spectral_grid.nlayers

        # surface pressure varying in space, as it would be over orography
        pₛ = zeros(Float64, spectral_grid.grid)
        for ij in eachgridpoint(pₛ)
            pₛ[ij] = 1000.0e2 - 200.0e2 * (ij / length(pₛ))
        end

        # θ linear in ln(p) everywhere, so the isentrope's pressure is known exactly,
        # inside the column and extrapolated below it
        θ_of_p(p) = 300 + 40 * log(p₀ / p)
        p_of_θ(θ) = p₀ * exp(-(θ - 300) / 40)
        temperature = isentropic_test_temperature(spectral_grid, coordinates, pₛ, (k, p) -> θ_of_p(p), κ, p₀)

        θ = [302.0, 320.0, 350.0, 400.0]
        out = zeros(Float64, spectral_grid.grid, length(θ))

        # a field linear in ln(p) is reproduced exactly at the isentrope's pressure p*,
        # and θ itself gives the target θ
        f(p) = 2log(p) + 3
        in_field = isentropic_test_field(spectral_grid, coordinates, pₛ, f)
        interpolate_isentropic_layers!(out, in_field, temperature, pₛ, θ, coordinates, κ, p₀)
        θ_field = isentropic_test_field(spectral_grid, coordinates, pₛ, θ_of_p)
        out_θ = zeros(Float64, spectral_grid.grid, length(θ))
        interpolate_isentropic_layers!(out_θ, θ_field, temperature, pₛ, θ, coordinates, κ, p₀)

        n = 0
        for k in eachindex(θ), ij in eachgridpoint(out)
            p_top = pressure(1, pₛ[ij], coordinates)
            p_bottom = pressure(nlayers, pₛ[ij], coordinates)
            p_top <= p_of_θ(θ[k]) <= p_bottom || continue       # inside the column
            @test out[ij, k] ≈ f(p_of_θ(θ[k])) rtol = 1.0e-10
            @test out_θ[ij, k] ≈ θ[k] rtol = 1.0e-10
            n += 1
        end
        @test n > 0                         # don't pass vacuously

        # isentropes equal to the θ of the model layers return the layer values,
        # constant pₛ so that they are the same in every column
        pₛ_const = fill!(zeros(Float64, spectral_grid.grid), 1000.0e2)
        temperature_const = isentropic_test_temperature(spectral_grid, coordinates, pₛ_const, (k, p) -> θ_of_p(p), κ, p₀)
        in_field_const = zeros(Float64, spectral_grid.grid, nlayers)
        for k in 1:nlayers, ij in eachgridpoint(in_field_const)
            in_field_const[ij, k] = k^2 + ij
        end
        θ_layers = [θ_of_p(pressure(k, 1000.0e2, coordinates)) for k in 1:nlayers]
        out_layers = zeros(Float64, spectral_grid.grid, nlayers)
        interpolate_isentropic_layers!(out_layers, in_field_const, temperature_const, pₛ_const, θ_layers, coordinates, κ, p₀)
        @test out_layers.data ≈ in_field_const.data rtol = 1.0e-10
    end

    @testset "Extrapolation as for pressure layers" begin
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8, NF = Float64)
        coordinates = SigmaCoordinates(spectral_grid)
        nlayers = spectral_grid.nlayers

        pₛ = fill!(zeros(Float64, spectral_grid.grid), 1000.0e2)
        θ_of_p(p) = 300 + 40 * log(p₀ / p)
        p_of_θ(θ) = p₀ * exp(-(θ - 300) / 40)
        temperature = isentropic_test_temperature(spectral_grid, coordinates, pₛ, (k, p) -> θ_of_p(p), κ, p₀)

        p_top = pressure(1, 1000.0e2, coordinates)
        p_bottom = pressure(nlayers, 1000.0e2, coordinates)

        # 1: above the model top, 2: below the lowest layer but above ground, 3: below ground
        θ = [θ_of_p(p_top / 2), θ_of_p((p_bottom + 1000.0e2) / 2), θ_of_p(1010.0e2)]
        @test p_of_θ(θ[2]) > p_bottom
        out = zeros(Float64, spectral_grid.grid, length(θ))

        # temperature as the field to interpolate, as TemperatureOutput would
        top = temperature[:, 1]
        bottom = temperature[:, nlayers]

        interpolate_isentropic_layers!(out, temperature, temperature, pₛ, θ, coordinates, κ, p₀, ConstantExtrapolation())
        @test out[:, 1] == top                      # outer-most layers held constant
        @test out[:, 2] == bottom
        @test out[:, 3] == bottom

        # dry adiabatic descent from the lowest layer to p*, extrapolated with θ linear in ln(p),
        # below ground the Lorenz condition p* = pₛ
        adiabatic = DryAdiabaticExtrapolation(κ)
        interpolate_isentropic_layers!(out, temperature, temperature, pₛ, θ, coordinates, κ, p₀, adiabatic)
        @test out[:, 1] == top                      # unchanged above the model top
        @test out[:, 2] ≈ bottom .* (p_of_θ(θ[2]) / p_bottom)^κ rtol = 1.0e-10
        @test out[:, 3] ≈ bottom .* (1000.0e2 / p_bottom)^κ rtol = 1.0e-10

        # a weakly stratified lowest layer pair extrapolates p* far below ground: still the
        # surface value, finite, and masked by a SubsurfaceMask
        θ_weak(k, p) = k == nlayers ? θ_of_p(pressure(nlayers - 1, 1000.0e2, coordinates)) - 1.0e-3 : θ_of_p(p)
        temperature_weak = isentropic_test_temperature(spectral_grid, coordinates, pₛ, θ_weak, κ, p₀)
        θ_deep = [θ_of_p(1010.0e2) - 30]
        out_deep = zeros(Float64, spectral_grid.grid, 1)
        interpolate_isentropic_layers!(out_deep, temperature_weak, temperature_weak, pₛ, θ_deep, coordinates, κ, p₀, adiabatic)
        @test out_deep[:, 1] ≈ temperature_weak[:, nlayers] .* (1000.0e2 / p_bottom)^κ rtol = 1.0e-10
        interpolate_isentropic_layers!(
            out_deep, temperature_weak, temperature_weak, pₛ, θ_deep, coordinates, κ, p₀,
            SubsurfaceMask(above_surface = adiabatic),
        )
        @test all(isnan, out_deep[:, 1])

        # masked exactly where p* is below ground
        interpolate_isentropic_layers!(
            out, temperature, temperature, pₛ, θ, coordinates, κ, p₀,
            SubsurfaceMask(above_surface = adiabatic),
        )
        @test out[:, 1] == top                      # above the model top, not masked
        @test out[:, 2] ≈ bottom .* (p_of_θ(θ[2]) / p_bottom)^κ rtol = 1.0e-10
        @test all(isnan, out[:, 3])                 # below ground, masked
    end

    @testset "Unstable columns" begin
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8, NF = Float64)
        coordinates = SigmaCoordinates(spectral_grid)
        nlayers = spectral_grid.nlayers
        pₛ = fill!(zeros(Float64, spectral_grid.grid), 1000.0e2)
        p = [pressure(k, 1000.0e2, coordinates) for k in 1:nlayers]
        ln_p = isentropic_test_field(spectral_grid, coordinates, pₛ, log)

        # superadiabatic below layer 5: θ crosses 302 K between layers 3, 4 and again 6, 7
        θ_column = [400.0, 350.0, 320.0, 300.0, 295.0, 300.0, 305.0, 310.0]
        temperature = isentropic_test_temperature(spectral_grid, coordinates, pₛ, (k, p) -> θ_column[k], κ, p₀)
        out = zeros(Float64, spectral_grid.grid, 1)
        interpolate_isentropic_layers!(out, ln_p, temperature, pₛ, [302.0], coordinates, κ, p₀)
        w = (302 - 320) / (300 - 320)               # first crossing from the top
        @test all(out[:, 1] .≈ (1 - w) * log(p[3]) + w * log(p[4]))

        # unstable outermost pairs: Lorenz condition, p* = pₛ below and p* = p_top above
        θ_column = [340.0, 350.0, 330.0, 320.0, 310.0, 300.0, 290.0, 295.0]
        temperature = isentropic_test_temperature(spectral_grid, coordinates, pₛ, (k, p) -> θ_column[k], κ, p₀)
        θ = [360.0, 280.0]                          # above every θ, below every θ in the column
        out = zeros(Float64, spectral_grid.grid, length(θ))
        adiabatic = DryAdiabaticExtrapolation(κ)
        interpolate_isentropic_layers!(out, temperature, temperature, pₛ, θ, coordinates, κ, p₀, adiabatic)
        @test out[:, 1] == temperature[:, 1]
        @test out[:, 2] ≈ temperature[:, nlayers] .* (1000.0e2 / p[nlayers])^κ rtol = 1.0e-10

        # p* = pₛ is at, not below, the surface: not masked
        interpolate_isentropic_layers!(
            out, temperature, temperature, pₛ, θ, coordinates, κ, p₀,
            SubsurfaceMask(above_surface = adiabatic),
        )
        @test all(isfinite, out)
        @test out[:, 2] ≈ temperature[:, nlayers] .* (1000.0e2 / p[nlayers])^κ rtol = 1.0e-10
    end

    @testset "Dimension mismatches" begin
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8)
        coordinates = SigmaCoordinates(spectral_grid)
        NF = spectral_grid.NF

        pₛ = fill!(zeros(NF, spectral_grid.grid), 1000.0e2)
        in_field = 250 .+ zeros(NF, spectral_grid.grid, spectral_grid.nlayers)
        θ = NF[300, 330]

        # number of isentropic layers doesn't match the output field
        out = zeros(NF, spectral_grid.grid, 3)
        @test_throws DimensionMismatch interpolate_isentropic_layers!(out, in_field, in_field, pₛ, θ, coordinates, κ, p₀)

        # output field on a different grid
        out = zeros(NF, FullGaussianGrid(spectral_grid.grid.nlat_half), length(θ))
        @test_throws DimensionMismatch interpolate_isentropic_layers!(out, in_field, in_field, pₛ, θ, coordinates, κ, p₀)

        # temperature with a different number of layers
        out = zeros(NF, spectral_grid.grid, length(θ))
        temperature = 250 .+ zeros(NF, spectral_grid.grid, 4)
        @test_throws DimensionMismatch interpolate_isentropic_layers!(out, in_field, temperature, pₛ, θ, coordinates, κ, p₀)

        # a single model layer has no layer pair
        spectral_grid1 = SpectralGrid(truncation = 15, nlayers = 1)
        in_field1 = 250 .+ zeros(NF, spectral_grid1.grid, 1)
        pₛ1 = fill!(zeros(NF, spectral_grid1.grid), 1000.0e2)
        out1 = zeros(NF, spectral_grid1.grid, length(θ))
        @test_throws ArgumentError interpolate_isentropic_layers!(
            out1, in_field1, in_field1, pₛ1, θ, SigmaCoordinates(spectral_grid1), κ, p₀,
        )
    end
end

@testset "Output on isentropic layers" begin

    @testset "IsentropicLayers construction" begin
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8)

        layers = IsentropicLayers(spectral_grid)
        @test layers.potential_temperature == SpeedyWeather.DEFAULT_ISENTROPIC_LAYERS
        @test SpeedyWeather.get_nlayers(layers, spectral_grid) == length(layers.potential_temperature)
        @test SpeedyWeather.vertical_dimension_name(layers) == "potential_temperature"

        # potential temperature and scratch field in the model's number format
        @test eltype(layers.potential_temperature) == spectral_grid.NF
        @test size(layers.scratch) == (SpeedyWeather.get_npoints(spectral_grid.grid), length(layers.potential_temperature))

        θ = [300, 315, 330]
        output = NetCDFOutput(spectral_grid, PrimitiveDry, layers = IsentropicLayers(spectral_grid, θ))
        @test size(output.field3D, 2) == length(θ)

        # non-monotonic or negative isentropic layers are rejected at construction
        @test_throws AssertionError IsentropicLayers(spectral_grid, [300, 350, 330])
        @test_throws AssertionError IsentropicLayers(spectral_grid, [-300, 330])

        # the output layers select the same step as the output writer
        model = PrimitiveDryModel(spectral_grid; output)
        simulation = initialize!(model)
        temperature = simulation.variables.grid.temperature
        @test SpeedyWeather.get_prognostic_step(temperature, model.time_stepping, output.layers) ==
            SpeedyWeather.get_prognostic_step(temperature, model.time_stepping, output)

        # models without temperature are rejected when the output is initialized
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 1)
        output = NetCDFOutput(spectral_grid, Barotropic, layers = IsentropicLayers(spectral_grid), path = mktempdir())
        model = BarotropicModel(spectral_grid; output)
        simulation = initialize!(model)
        @test_throws ArgumentError run!(simulation, steps = 1, output = true)
    end

    @testset "NetCDFOutput on isentropic layers" begin
        # full grid model so that the horizontal interpolation onto the (full) output grid
        # is an exact copy, which lets us compare the file against the interpolation directly
        spectral_grid = SpectralGrid(truncation = 15, nlayers = 8, Grid = FullGaussianGrid)
        NF = spectral_grid.NF
        θ = [290, 300, 315, 330]

        output = NetCDFOutput(
            spectral_grid, PrimitiveDry,
            path = mktempdir(), write_restart = false, interval = Hour(6),
            layers = IsentropicLayers(spectral_grid, θ),
        )
        # no bitrounding so that the comparison below is exact
        add!(output, SpeedyWeather.TemperatureOutput(keepbits = 23), SpeedyWeather.ZonalVelocityOutput(keepbits = 23))

        model = PrimitiveDryModel(spectral_grid; output)
        simulation = initialize!(model)
        run!(simulation, period = Hour(6), output = true)
        @test model.feedback.nans_detected == false

        NCDataset(SpeedyWeather.get_output_path(model.output)) do ds
            # the potential temperature dimension replaces the sigma layer dimension
            @test haskey(ds.dim, "potential_temperature")
            @test !haskey(ds.dim, "layer")
            @test ds["potential_temperature"][:] ≈ θ
            @test ds["potential_temperature"].attrib["units"] == "K"
            @test ds["potential_temperature"].attrib["standard_name"] == "air_potential_temperature"

            # 3D atmospheric variables are on it, 2D variables are not
            for name in ("temp", "u", "v", "vor")
                @test NCDatasets.dimnames(ds[name]) == ("lon", "lat", "potential_temperature", "time")
            end
            @test NCDatasets.dimnames(ds["mslp"]) == ("lon", "lat", "time")

            # values match interpolating the final state directly, temperature with its
            # dry-adiabatic extrapolation, u constant
            vars = simulation.variables
            step(var) = SpeedyWeather.get_prognostic_step(var, model.time_stepping, model.output)
            temperature = step(vars.grid.temperature)
            pₛ = vars.parameterizations.surface_pressure
            nlon, nlat = length(ds["lon"]), length(ds["lat"])

            for (name, field, offset) in (("temp", temperature, 273.15), ("u", step(vars.grid.u), 0))
                reference = zeros(NF, spectral_grid.grid, length(θ))
                interpolate_isentropic_layers!(
                    reference, field, temperature, pₛ, model.output.layers.potential_temperature,
                    model.geometry.vertical_coordinates, model.atmosphere.κ,
                    model.atmosphere.reference_pressure,
                    SpeedyWeather.output_extrapolation(model.output.variables[Symbol(name)]),
                )
                from_file = ds[name][:, :, :, end] .+ offset         # ˚C back to K
                @test from_file ≈ reshape(Array(reference.data), nlon, nlat, length(θ)) rtol = 1.0e-5
                @test all(isfinite, from_file)
            end
        end
    end
end
