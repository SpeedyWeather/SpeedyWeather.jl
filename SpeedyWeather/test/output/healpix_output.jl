using Zarr, Dates
import SpeedyWeather.RingGrids

"""Latitudes (˚N), longitudes (˚E) and 1-based ring indices of every HEALPix pixel of
resolution `nside` in RING order, from the canonical pixel-centre formulas of Górski et al.
(2005), eqs. 4/5 (polar caps) and 8/9 (equatorial belt), in their `z = cos(θ)` form.

Deliberately independent of RingGrids' own `get_latd`/`get_lond_per_ring` (which use the
`90 - acosd(...)` form and derive the ring lengths from `min(4j, 2nlat_half, 8nlat_half-4j)`)
so that the comparison actually tests the HEALPix layout rather than restating it."""
function healpix_ring_reference(nside::Integer)
    npix = 12nside^2
    latds, londs = zeros(npix), zeros(npix)
    rings = zeros(Int, npix)

    p = 1                                       # running flat (RING) index, 1-based
    for i in 1:(4nside - 1)                     # ring index, north to south
        if i < nside                            # north polar cap, eq 4 and 5
            nlon = 4i
            z = 1 - i^2 / (3nside^2)
            φ = [(180 / (2i)) * (j - 0.5) for j in 1:nlon]
        elseif i <= 3nside                      # equatorial belt, eq 8 and 9
            nlon = 4nside
            z = 4 / 3 - 2i / (3nside)
            s = mod(i - nside + 1, 2) == 1 ? 1 : 2      # half-offset on alternating rings
            φ = [(180 / (2nside)) * (j - s / 2) for j in 1:nlon]
        else                                    # south polar cap, mirrored from the north
            i′ = 4nside - i
            nlon = 4i′
            z = -(1 - i′^2 / (3nside^2))
            φ = [(180 / (2i′)) * (j - 0.5) for j in 1:nlon]
        end

        latds[p:(p + nlon - 1)] .= asind(z)
        londs[p:(p + nlon - 1)] .= φ
        rings[p:(p + nlon - 1)] .= i
        p += nlon
    end
    @assert p == npix + 1
    return londs, latds, rings
end

@testset "HEALPixOutput type and defaults" begin
    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)
    output = HEALPixOutput(spectral_grid)
    @test output isa SpeedyWeather.HEALPixOutput
    @test output.active == false
    @test output.filename == "output_healpix.zarr"

    # default resolution follows the model grid, output grid is always HEALPix on CPU
    grid = output.field2D.grid
    @test grid isa HEALPixGrid
    @test grid.nlat_half == spectral_grid.grid.nlat_half

    # nside and nlat_half are two spellings of the same resolution
    @test HEALPixOutput(spectral_grid, nside = 4).field2D.grid.nlat_half == 8
    @test HEALPixOutput(spectral_grid, nlat_half = 8).field2D.grid.nlat_half == 8
    @test RingGrids.get_npoints(HEALPixOutput(spectral_grid, nside = 4).field2D) == 12 * 4^2

    # HEALPixGrid is only defined for even nlat_half, odd requests round up
    @test HEALPixOutput(spectral_grid, nlat_half = 7).field2D.grid.nlat_half == 8

    # ambiguous resolution requests are rejected rather than silently resolved
    @test_throws ArgumentError HEALPixOutput(spectral_grid, nside = 4, nlat_half = 8)
    @test_throws ArgumentError HEALPixOutput(
        spectral_grid, nside = 4, output_grid = HEALPixGrid(8)
    )
end

@testset "HEALPixOutput output grid selection" begin
    # the output grid follows the *model's* grid type when the model already runs on one of
    # the two supported grids, so those simulations skip the interpolation by default
    for Grid in (HEALPixGrid, OctaHEALPixGrid)
        spectral_grid = SpectralGrid(truncation = 6, nlayers = 1; Grid)
        output = HEALPixOutput(spectral_grid, ShallowWater)
        @test output.field2D.grid isa Grid
        @test output.field2D.grid.nlat_half == spectral_grid.grid.nlat_half
        @test isnothing(output.interpolator)                # interpolation skipped
    end

    # any other model grid defaults to HEALPixGrid and does need an interpolator
    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)
    @test !(spectral_grid.grid isa Union{HEALPixGrid, OctaHEALPixGrid})
    @test HEALPixOutput(spectral_grid).field2D.grid isa HEALPixGrid
    @test !isnothing(HEALPixOutput(spectral_grid).interpolator)

    # `nside` is a HEALPix-only parameter and always selects a HEALPixGrid, even when the
    # model runs on OctaHEALPix (which then needs an interpolator again)
    octa_grid = SpectralGrid(truncation = 6, nlayers = 1, Grid = OctaHEALPixGrid)
    output = HEALPixOutput(octa_grid, nside = 4)
    @test output.field2D.grid isa HEALPixGrid
    @test output.field2D.grid.nlat_half == 8
    @test !isnothing(output.interpolator)

    # OctaHEALPixGrid takes odd nlat_half, unlike HEALPixGrid
    @test HEALPixOutput(spectral_grid, output_grid = OctaHEALPixGrid(7)).field2D.grid.nlat_half == 7
    @test HEALPixOutput(octa_grid, nlat_half = 7).field2D.grid.nlat_half == 7
end

@testset "HEALPixOutput flat write indices" begin
    # the x/y dims of a variable collapse into one flat `cell` dimension, so the flat
    # indices are the rectangular ones with one leading colon removed
    vor = SpeedyWeather.VorticityOutput()                       # 3D + time
    mslp = SpeedyWeather.MeanSeaLevelPressureOutput()           # 2D + time
    @test SpeedyWeather.get_flat_indices(3, vor) == (:, :, 3)
    @test SpeedyWeather.get_flat_indices(3, mslp) == (:, 3)
    @test SpeedyWeather.get_flat_indices(2:4, vor) == (:, :, 2:4)
    @test length(SpeedyWeather.get_flat_indices(3, vor)) ==
        length(SpeedyWeather.get_indices(3, vor)) - 1
end

@testset "HEALPixOutput coordinates are canonical HEALPix RING" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_healpixtests_coords_")
    nside = 4

    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)
    output = HEALPixOutput(
        spectral_grid, ShallowWater;
        path = tmp_output_path, write_restart = false, nside,
    )
    model = ShallowWaterModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true, period = Hour(6))

    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))
    npix = 12nside^2

    # the three flat coordinate vectors, one entry per grid point
    @test size(g["lat"]) == (npix,)
    @test size(g["lon"]) == (npix,)
    @test size(g["ring"]) == (npix,)
    @test g["cell"][:] == collect(1:npix)

    # they are the canonical HEALPix RING pixel centres
    londs_ref, latds_ref, rings_ref = healpix_ring_reference(nside)
    @test g["lat"][:] ≈ latds_ref
    @test g["lon"][:] ≈ londs_ref
    @test g["ring"][:] == rings_ref

    # RING invariants: 4nside-1 rings north to south, 4i / 4nside / 4(4nside-i) points each
    rings = g["ring"][:]
    @test extrema(rings) == (1, 4nside - 1)
    @test issorted(rings)                               # cells are grouped ring by ring
    ring_lengths = [count(==(i), rings) for i in 1:(4nside - 1)]
    @test ring_lengths == [
        i < nside ? 4i : i <= 3nside ? 4nside : 4(4nside - i) for i in 1:(4nside - 1)
    ]

    # latitude decreases monotonically from north to south, ring by ring
    ring_lats = [g["lat"][findfirst(==(i), rings)] for i in 1:(4nside - 1)]
    @test issorted(ring_lats, rev = true)
    @test all(-90 .< ring_lats .< 90)                   # HEALPix has no pixel at the poles

    # the store is self-describing for healpy/cuHPX consumers
    @test g.attrs["healpix_nside"] == nside
    @test g.attrs["healpix_npix"] == npix
    @test g.attrs["healpix_order"] == "RING"
end

@testset "HEALPixOutput for ShallowWaterModel, interpolating" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_healpixtests_sw_")
    period = Day(1)
    nside = 4

    # default octahedral Gaussian model grid, so the data has to be interpolated
    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)
    output = HEALPixOutput(
        spectral_grid, ShallowWater;
        path = tmp_output_path, write_restart = false, nside,
    )
    @test !isnothing(output.interpolator)   # interpolation is needed here

    model = ShallowWaterModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true; period)
    @test simulation.model.feedback.nans_detected == false

    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))
    npix = 12nside^2

    # all declared output variables and all coordinates made it into the store
    for var in values(output.variables)
        @test haskey(g.arrays, var.name)
    end
    for c in ("cell", "lat", "lon", "ring", "layer", "time")
        @test haskey(g.arrays, c)
    end

    expected_times = Int(period / output.interval) + 1
    @test length(g["time"][:]) == expected_times
    @test g["time"][1] == 0.0

    # data is flat: (cell, layer, time) for 3D, (cell, time) for 2D — one horizontal dim
    @test size(g["vor"]) == (npix, spectral_grid.nlayers, expected_times)
    @test size(g["eta"]) == (npix, expected_times)

    # `_ARRAY_DIMENSIONS` is row-major (xarray order), i.e. reversed w.r.t. the Julia shape
    @test g["vor"].attrs["_ARRAY_DIMENSIONS"] == ["time", "layer", "cell"]
    @test g["eta"].attrs["_ARRAY_DIMENSIONS"] == ["time", "cell"]
    @test g["lat"].attrs["_ARRAY_DIMENSIONS"] == ["cell"]
    @test g["ring"].attrs["_ARRAY_DIMENSIONS"] == ["cell"]

    # values are finite (not stuck on the fill value)
    @test all(isfinite, g["vor"][:, :, :])
    @test all(isfinite, g["eta"][:, :])
end

@testset "HEALPixOutput skips interpolation on a HEALPix simulation" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_healpixtests_skip_")
    period = Hour(12)

    # model already runs on a HEALPixGrid, so no interpolation should be needed
    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1, Grid = HEALPixGrid)
    @test spectral_grid.grid isa HEALPixGrid

    output = HEALPixOutput(
        spectral_grid, ShallowWater;
        path = tmp_output_path, write_restart = false,
    )
    # the feature: same grid in and out ⇒ no interpolator is built at all
    @test isnothing(output.interpolator)
    @test RingGrids.grids_match(output.field2D.grid, spectral_grid.grid)

    # keep all mantissa bits so the write path is lossless and we can compare exactly;
    # any interpolation would blend neighbouring cells and break the equality below
    output.variables[:u].keepbits = 23

    model = ShallowWaterModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true; period)
    @test simulation.model.feedback.nans_detected == false

    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))

    # coordinates are the model grid's own, not a regridded set
    londs, latds = RingGrids.get_londlatds(spectral_grid.grid)
    @test g["lat"][:] ≈ latds
    @test g["lon"][:] ≈ londs
    @test g["ring"][:] == RingGrids.whichring(spectral_grid.grid)

    # and the last written slice is bit-identical to the model's own grid field (grid.u
    # carries a leapfrog step dimension, pick the same step the writer does)
    u = simulation.variables.grid.u
    u_model = Array(parent(SpeedyWeather.get_prognostic_step(u, model.time_stepping, output)))
    @test g["u"][:, :, end] == Float32.(u_model)
end

@testset "HEALPixOutput on an OctaHEALPixGrid" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_healpixtests_octa_")
    period = Hour(12)

    # OctaHEALPixGrid uses the *same* flat store convention as HEALPixGrid, but is a
    # SpeedyWeather-specific tessellation that healpy/cuHPX cannot read, so it must not
    # advertise any healpix_* metadata
    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1, Grid = OctaHEALPixGrid)
    output = HEALPixOutput(
        spectral_grid, ShallowWater;
        path = tmp_output_path, write_restart = false,
    )
    @test output.field2D.grid isa OctaHEALPixGrid
    @test isnothing(output.interpolator)        # model already on this grid
    output.variables[:u].keepbits = 23          # lossless, so we can compare exactly

    model = ShallowWaterModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true; period)
    @test simulation.model.feedback.nans_detected == false

    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))
    grid = spectral_grid.grid
    npix = RingGrids.get_npoints(grid)
    @test npix == 4 * grid.nlat_half^2          # OctaHEALPix has 4nlat_half^2 points

    # same store convention as a HEALPixGrid store: flat cell dim, flat lat/lon/ring
    @test size(g["vor"]) == (npix, spectral_grid.nlayers, Int(period / output.interval) + 1)
    @test g["vor"].attrs["_ARRAY_DIMENSIONS"] == ["time", "layer", "cell"]
    @test g["cell"][:] == collect(1:npix)
    for c in ("lat", "lon", "ring")
        @test size(g[c]) == (npix,)
        @test g[c].attrs["_ARRAY_DIMENSIONS"] == ["cell"]
    end

    # coordinates are the model grid's own, in the same north-to-south ring order
    londs, latds = RingGrids.get_londlatds(grid)
    @test g["lat"][:] ≈ latds
    @test g["lon"][:] ≈ londs
    @test g["ring"][:] == RingGrids.whichring(grid)
    @test extrema(g["ring"][:]) == (1, RingGrids.get_nlat(grid))
    @test issorted(g["ring"][:])

    # same metadata key set as a HEALPixGrid store, with `nside` generalized to the side
    # length of one square face so that npix == nfaces * nside^2 holds here too
    @test g.attrs["grid"] == "OctaHEALPixGrid"
    @test g.attrs["npix"] == npix
    @test g.attrs["nlat_half"] == grid.nlat_half
    @test g.attrs["equal_area"] == true
    @test g.attrs["healpix_order"] == "RING"
    @test g.attrs["healpix_npix"] == npix
    @test g.attrs["healpix_nside"] == grid.nlat_half        # OctaHEALPix faces are nlat_half wide
    @test g.attrs["healpix_nfaces"] == 4                    # vs 12 for standard HEALPix
    @test g.attrs["healpix_npix"] == g.attrs["healpix_nfaces"] * g.attrs["healpix_nside"]^2
    @test g["cell"].attrs["ordering"] == "RING"
    @test g["cell"].attrs["nside"] == grid.nlat_half
    @test g["cell"].attrs["nfaces"] == 4

    # `healpix_nfaces` (and `grid`) is what tells a reader this is not standard HEALPix:
    # assuming the usual 12 faces would give the wrong npix and silently misread the store
    @test g.attrs["healpix_nfaces"] != 12
    @test 12 * g.attrs["healpix_nside"]^2 != npix
    @test haskey(g.attrs, "note")                           # spells out the incompatibility

    # and the data still round-trips losslessly through the no-interpolation path
    u = simulation.variables.grid.u
    u_model = Array(parent(SpeedyWeather.get_prognostic_step(u, model.time_stepping, output)))
    @test g["u"][:, :, end] == Float32.(u_model)
end

@testset "HEALPixOutput HEALPixGrid store advertises HEALPix metadata" begin
    # counterpart to the OctaHEALPix test above: a true HEALPixGrid store *must* carry the
    # healpix_* keys healpy/cuHPX look for
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_healpixtests_meta_")
    nside = 4

    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)
    output = HEALPixOutput(
        spectral_grid, ShallowWater;
        path = tmp_output_path, write_restart = false, nside,
    )
    model = ShallowWaterModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true, period = Hour(6))

    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))
    @test g.attrs["grid"] == "HEALPixGrid"
    @test g.attrs["healpix_nside"] == nside
    @test g.attrs["healpix_npix"] == 12nside^2
    @test g.attrs["healpix_nfaces"] == 12                   # the standard tessellation
    @test g.attrs["healpix_order"] == "RING"
    @test g["cell"].attrs["ordering"] == "RING"
    @test g["cell"].attrs["nside"] == nside
    @test g["cell"].attrs["nfaces"] == 12
    @test !haskey(g.attrs, "note")                          # nothing to warn about here

    # npix == nfaces * nside^2 is the invariant that holds across both output grids
    @test g.attrs["healpix_npix"] == g.attrs["healpix_nfaces"] * g.attrs["healpix_nside"]^2
end

@testset "HEALPixOutput time chunking" begin
    period = Day(1)
    nside = 2

    # a time_chunk that does not divide the number of outputs leaves a trailing partial
    # chunk, which has to be flushed on finalize! — compare against an unbuffered run
    function run_healpix(time_chunk)
        tmp = mktempdir(pwd(), prefix = "tmp_healpixtests_chunk_")
        spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)
        output = HEALPixOutput(
            spectral_grid, ShallowWater;
            path = tmp, write_restart = false, nside, time_chunk,
        )
        model = ShallowWaterModel(spectral_grid; output)
        simulation = initialize!(model)
        run!(simulation, output = true; period)
        return Zarr.zopen(joinpath(model.output.run_path, model.output.filename))
    end

    reference = run_healpix(1)
    chunked = run_healpix(3)                    # 5 outputs = one full chunk + 2 buffered

    @test size(chunked["vor"]) == size(reference["vor"])
    @test chunked["vor"][:, :, :] == reference["vor"][:, :, :]
    @test chunked["time"][:] == reference["time"][:]
    @test all(isfinite, chunked["vor"][:, :, end])      # trailing partial chunk flushed
end

@testset "HEALPixOutput for PrimitiveWetModel with soil layers" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_healpixtests_pw_")
    period = Day(1)
    nside = 4

    spectral_grid = SpectralGrid(truncation = 6, nlayers = 4)
    output = HEALPixOutput(
        spectral_grid, PrimitiveWet;
        path = tmp_output_path, write_restart = false, nside,
    )
    model = PrimitiveWetModel(spectral_grid; output)
    simulation = initialize!(model)
    run!(simulation, output = true; period)
    @test simulation.model.feedback.nans_detected == false

    g = Zarr.zopen(joinpath(model.output.run_path, model.output.filename))
    npix = 12nside^2
    expected_times = Int(period / output.interval) + 1

    @test haskey(g.arrays, "temp")
    @test haskey(g.arrays, "humid")
    @test haskey(g.arrays, "soil_layer")
    @test size(g["temp"]) == (npix, spectral_grid.nlayers, expected_times)

    # mslp, u10/v10 and tsurf derive a 2D field from 3D model state in their own `output!`
    # methods; make sure they all reach the flat backend. Skip the t=0 snapshot for mslp —
    # physics hasn't run yet so the surface-air temperature buffer is still uninitialized.
    @test size(g["mslp"]) == (npix, expected_times)
    @test all(isfinite, g["mslp"][:, 2:end])
end
