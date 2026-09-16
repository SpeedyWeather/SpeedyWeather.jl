import SpeedyWeather: interpolate_pressure_levels!, pressure
import SpeedyWeather: LinearInPressure, LinearInLogPressure
import SpeedyWeather: ConstantExtrapolation, DryAdiabaticExtrapolation, SubsurfaceMask

# a field that is an exact function f of pressure, so that interpolating it onto pressure
# levels can be checked against f evaluated at those levels
function pressure_test_field(spectral_grid, coordinates, pₛ, f::Function)
    field = zeros(spectral_grid.NF, spectral_grid.grid, spectral_grid.nlayers)
    for k in 1:spectral_grid.nlayers, ij in eachgridpoint(field)
        field[ij, k] = f(pressure(k, pₛ[ij], coordinates))
    end
    return field
end

# maximum relative error of `out` against `f(p)`, for all grid points and pressure levels
# that are inside the model levels (outside is extrapolation, tested separately).
# Returns the error and the number of points checked.
function max_relative_error(out, p, pₛ, coordinates, nlayers, f::Function)
    error, n = zero(eltype(out)), 0
    for k in eachindex(p), ij in eachgridpoint(out)
        pressure(1, pₛ[ij], coordinates) <= p[k] <= pressure(nlayers, pₛ[ij], coordinates) || continue
        error = max(error, abs(out[ij, k] - f(p[k])) / abs(f(p[k])))
        n += 1
    end
    return error, n
end

@testset "Vertical interpolation onto pressure levels" begin
    @testset "Interpolation with $VerticalCoordinates" for VerticalCoordinates in
        (SigmaCoordinates, SigmaPressureCoordinates)

        spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
        coordinates = VerticalCoordinates(spectral_grid)
        NF = spectral_grid.NF
        nlayers = spectral_grid.nlayers

        # surface pressure varying in space, as it would be over orography
        pₛ = zeros(NF, spectral_grid.grid)
        for ij in eachgridpoint(pₛ)
            pₛ[ij] = 1000.0e2 - 200.0e2 * (ij / length(pₛ))
        end

        p = NF[200.0e2, 400.0e2, 600.0e2, 800.0e2]
        out = zeros(NF, spectral_grid.grid, length(p))

        # a field linear in p is reproduced exactly by linear interpolation in p,
        # a field linear in log(p) exactly by linear interpolation in log(p)
        for (interpolation, f) in (
                (LinearInPressure(), p -> 2p + 3),
                (LinearInLogPressure(), p -> 2log(p) + 3),
            )
            in_field = pressure_test_field(spectral_grid, coordinates, pₛ, f)
            interpolate_pressure_levels!(out, in_field, pₛ, p, coordinates, interpolation)
            error, n = max_relative_error(out, p, pₛ, coordinates, nlayers, f)
            @test n > 0                     # don't pass vacuously
            @test error < 10eps(NF)
        end

        # interpolating onto pressure levels that coincide with the model levels returns
        # the model level values, constant pₛ so that they are the same in every column
        in_field = zeros(NF, spectral_grid.grid, nlayers)
        for k in 1:nlayers, ij in eachgridpoint(in_field)
            in_field[ij, k] = k^2 + ij
        end

        pₛ_const = fill!(zeros(NF, spectral_grid.grid), 1000.0e2)
        p_levels = NF[pressure(k, NF(1000.0e2), coordinates) for k in 1:nlayers]
        out_levels = zeros(NF, spectral_grid.grid, nlayers)

        for interpolation in (LinearInPressure(), LinearInLogPressure())
            interpolate_pressure_levels!(
                out_levels, in_field, pₛ_const, p_levels, coordinates, interpolation,
            )
            @test all(≈(0, atol = 10eps(NF)), (out_levels .- in_field) ./ in_field)
        end
    end

    @testset "Extrapolation beyond the model levels" begin
        spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
        coordinates = SigmaCoordinates(spectral_grid)
        NF = spectral_grid.NF
        nlayers = spectral_grid.nlayers
        κ = NF(2 / 7)

        pₛ = fill!(zeros(NF, spectral_grid.grid), 1000.0e2)
        in_field = zeros(NF, spectral_grid.grid, nlayers)
        for k in 1:nlayers, ij in eachgridpoint(in_field)
            in_field[ij, k] = 200 + 10k
        end

        p_top = pressure(1, NF(1000.0e2), coordinates)
        p_bottom = pressure(nlayers, NF(1000.0e2), coordinates)
        @test p_bottom < 1000.0e2     # lowest full level is above the surface

        # 1: above the model top, 2: below the lowest level but above ground, 3: below ground
        p = NF[p_top / 2, (p_bottom + 1000.0e2) / 2, 1010.0e2]
        out = zeros(NF, spectral_grid.grid, length(p))

        top = in_field[:, 1]
        bottom = in_field[:, nlayers]

        @testset "constant" begin
            interpolate_pressure_levels!(
                out, in_field, pₛ, p, coordinates,
                LinearInLogPressure(), ConstantExtrapolation(),
            )
            @test out[:, 1] == top          # outer-most levels held constant
            @test out[:, 2] == bottom
            @test out[:, 3] == bottom
        end

        @testset "dry adiabatic" begin
            interpolate_pressure_levels!(
                out, in_field, pₛ, p, coordinates,
                LinearInLogPressure(), DryAdiabaticExtrapolation(κ),
            )
            @test out[:, 1] == top                              # unchanged above the model top
            @test out[:, 2] ≈ bottom .* (p[2] / p_bottom)^κ     # descends adiabatically
            @test out[:, 3] ≈ bottom .* (p[3] / p_bottom)^κ
            @test all(out[:, 2] .> bottom)                      # and is warmer below
        end

        @testset "subsurface mask" begin
            interpolate_pressure_levels!(
                out, in_field, pₛ, p, coordinates, LinearInLogPressure(),
                SubsurfaceMask(above_surface = DryAdiabaticExtrapolation(κ)),
            )
            @test out[:, 1] == top                              # above the model top, not masked
            @test out[:, 2] ≈ bottom .* (p[2] / p_bottom)^κ     # above ground, extrapolated
            @test all(isnan, out[:, 3])                         # below ground, masked
        end
    end

    @testset "Dimension mismatches" begin
        spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
        coordinates = SigmaCoordinates(spectral_grid)
        NF = spectral_grid.NF

        pₛ = fill!(zeros(NF, spectral_grid.grid), 1000.0e2)
        in_field = zeros(NF, spectral_grid.grid, spectral_grid.nlayers)
        p = NF[500.0e2, 850.0e2]

        # number of pressure levels doesn't match the output field
        out = zeros(NF, spectral_grid.grid, 3)
        @test_throws DimensionMismatch interpolate_pressure_levels!(
            out, in_field, pₛ, p, coordinates,
        )

        # output field on a different grid
        out = zeros(NF, FullGaussianGrid(spectral_grid.grid.nlat_half), length(p))
        @test_throws DimensionMismatch interpolate_pressure_levels!(
            out, in_field, pₛ, p, coordinates,
        )

        # more pressure levels than model layers is fine, the output field decides
        p_many = NF[100.0e2 * i for i in 1:10]
        out = zeros(NF, spectral_grid.grid, length(p_many))
        interpolate_pressure_levels!(out, in_field, pₛ, p_many, coordinates)
        @test all(isfinite, out)
    end
end
