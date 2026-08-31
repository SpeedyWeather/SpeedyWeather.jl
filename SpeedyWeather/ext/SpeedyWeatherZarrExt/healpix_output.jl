# `HEALPixOutput`: writes output onto a `HEALPixGrid` keeping the horizontal dimension flat,
# with arrays shaped (cell, vertical, time). Cells are in HEALPix RING order so the store can
# be read by healpy/cuHPX without reordering. Everything that is not specific to that layout —
# time-chunk buffering, the time axis, coordinate helpers — lives in shared.jl.

"""$(TYPEDSIGNATURES)
Resolution parameter `nlat_half` of the output HEALPix grid from the user's `nside` or
`nlat_half` keyword, defaulting to the model grid's own `nlat_half`. `HEALPixGrid` is only
defined for even `nlat_half` (`nside = nlat_half ÷ 2`), so an odd value is rounded up."""
function healpix_nlat_half(model_nlat_half::Integer; nside = nothing, nlat_half = nothing)
    isnothing(nside) || isnothing(nlat_half) ||
        throw(ArgumentError("HEALPixOutput: pass either `nside` or `nlat_half`, not both."))

    isnothing(nside) || return 2 * Int(nside)
    n = Int(something(nlat_half, model_nlat_half))
    n > 0 || throw(ArgumentError("HEALPixOutput: nlat_half=$n must be positive."))
    iseven(n) && return n

    @info "HEALPixOutput: HEALPixGrid is only defined for even nlat_half, rounding $n up to $(n + 1) (nside=$((n + 1) ÷ 2))."
    return n + 1
end

"""$(TYPEDSIGNATURES)
Constructor for [`HEALPixOutput`](@ref) (extension version, available once `Zarr.jl` is
loaded). Builds the output `HEALPixGrid` and the scratch fields to write from, and pre-fills
`output.variables` with the defaults for `Model`.

The output resolution is given either as `nside`, as `nlat_half = 2nside`, or as a ready-made
`output_grid::HEALPixGrid`, and defaults to the model grid's `nlat_half` (rounded up to the
nearest even number, which `HEALPixGrid` requires). An interpolator onto that grid is only
built when it is actually needed: if the model already runs on the very same HEALPix grid,
`output.interpolator` stays `nothing` and the data is copied straight out, skipping
interpolation entirely."""
function HEALPixOutput(
        SG::SpectralGrid,
        Model::Type{<:AbstractModel} = Barotropic;
        nlayers_soil = DEFAULT_NLAYERS_SOIL,
        nside = nothing,
        nlat_half = nothing,
        output_grid::Union{HEALPixGrid, Nothing} = nothing,
        output_NF::DataType = DEFAULT_OUTPUT_NF,
        interval::Period = Second(DEFAULT_OUTPUT_INTERVAL),
        compressor = nothing,
        kwargs...
    )

    # OUTPUT GRID: from `output_grid`, or built from `nside`/`nlat_half`/the model's own
    # resolution. Always on CPU, we write from the host. The `::HEALPixGrid` annotation
    # rejects other (Octa)HEALPix-like grids, which healpy/cuHPX cannot read.
    output_grid = if isnothing(output_grid)
        HEALPixGrid(healpix_nlat_half(SG.grid.nlat_half; nside, nlat_half), CPU())
    else
        isnothing(nside) && isnothing(nlat_half) || throw(
            ArgumentError(
                "HEALPixOutput: pass either `output_grid` or `nside`/`nlat_half`, not both."
            )
        )
        on_architecture(CPU(), output_grid)
    end

    # INPUT GRID (but on CPU)
    input_grid = on_architecture(CPU(), SG.grid)

    # SKIP INTERPOLATION if the model already runs on this very HEALPix grid: `grids_match`
    # compares the (nonparametric) grid type and nlat_half, which is exactly the condition
    # under which `RingGrids.interpolate!` would degenerate to a copy anyway. Not building
    # the interpolator also saves precomputing its stencil indices and weights.
    interpolator = if grids_match(output_grid, input_grid)
        nothing
    else
        RingGrids.interpolator(output_grid, input_grid, NF = DEFAULT_OUTPUT_NF)
    end

    # CREATE HEALPIX FIELDS TO WRITE OUT FROM
    (; nlayers) = SG
    land_fraction = Field(output_NF, output_grid)
    field2D = Field(output_NF, output_grid)
    field3D = Field(output_NF, output_grid, nlayers)
    field3Dland = Field(output_NF, output_grid, nlayers_soil)

    # Concrete type parameters, see the ZarrOutput constructor for the compressor/group ones.
    # `Itp` is `Nothing` when interpolation is skipped, which makes the branch in
    # `interpolate_output!` a compile-time constant.
    C = typeof(resolve_compressor(compressor))
    Z = Zarr.ZGroup{Zarr.DirectoryStore}

    interval_sec = Second(interval)
    DT = DateTime
    S = typeof(interval_sec)
    F2 = typeof(field2D)
    F3 = typeof(field3D)
    Itp = typeof(interpolator)

    output = HEALPixOutput{F2, F3, Itp, DT, S, C, Z}(;
        interval = interval_sec,
        interpolator,
        land_fraction,
        field2D,
        field3D,
        field3Dland,
        compressor,
        kwargs...
    )

    add_default!(output.variables, Model)
    return output
end

# to dispatch over the dataset type
SpeedyWeather.dataset_type(::HEALPixOutput) = Zarr.ZGroup

# flat layout: the x/y dims of a variable collapse into the single `cell` dimension
zarr_write_indices(::HEALPixOutput, i, variable::AbstractOutputVariable) =
    get_flat_indices(i, variable)

"""$(TYPEDSIGNATURES)
Initialize `HEALPixOutput` by creating a Zarr group on disk and storing the initial
conditions of `vars`. The store keeps the horizontal dimension flat: variables are shaped
`(cell, layer/soil_layer, time)` with the leading `cell` dimension running over all
`12nside²` HEALPix cells in RING order, i.e. exactly the layout of a `Field`'s data.
Alongside `cell` (the flat index `1:npix`) the store holds the per-cell coordinate vectors
`lat`, `lon` and `ring`, plus the `layer`, `soil_layer` and `time` coordinates."""
function initialize!(
        output::HEALPixOutput,
        vars::Variables,
        model::AbstractModel,
    )
    output.active || return nothing

    assert_soil_layers(output, model)

    # SHARED INITIALIZATION (run folder, output frequency, counters, callbacks)
    initialize!(output.core, output, model)

    # Total number of output snapshots: IC + one per `output_every_n_steps`.
    n_outputs = vars.prognostic.clock.n_time_steps ÷ output.output_every_n_steps + 1

    # The Zarr store is a *directory*, not a single file.
    (; run_path, filename) = output
    store_path = joinpath(run_path, filename)
    output.overwrite && isdir(store_path) && rm(store_path; recursive = true)

    # global attributes so the store is self-describing for healpy/cuHPX consumers
    grid = output.field2D.grid
    g = Zarr.zgroup(
        store_path; attrs = Dict(
            "healpix_nside" => RingGrids.nside_healpix(grid.nlat_half),
            "healpix_npix" => get_npoints(grid),
            "healpix_order" => "RING",
        )
    )
    output.zarr_group = g

    # COORDINATES: the flat cell index and its lat/lon/ring, plus the vertical ones
    write_healpix_coordinates!(g, output, model)

    # TIME: full-length, chunked by `output.time_chunk`.
    create_time_axis!(g, output, n_outputs)
    output!(output, vars.prognostic.clock.time)   # write initial time

    # VARIABLES, remove output variables not existent in simulation.variables
    simulation = Simulation(vars, model)
    prune_nonexisting_variables!(output, simulation)

    # then define every output variable in the Zarr store and write initial conditions
    for (key, var) in output.variables
        define_variable!(g, output, var, n_outputs, eltype(output.field2D))
        output!(output, var, simulation)
    end

    # calculate land fraction on output grid (interpolated, or copied when skipping)
    if hasproperty(model, :land_sea_mask)
        land_fraction_cpu = on_architecture(CPU(), model.land_sea_mask.land_fraction)
        SpeedyWeather.interpolate_output!(output, output.land_fraction, land_fraction_cpu)
    end

    # consolidate the store metadata (.zmetadata) for faster opening with xarray etc.;
    # the schema is complete at this point and all later writes only touch chunk files
    Zarr.consolidate_metadata(g)

    return nothing
end

"""$(TYPEDSIGNATURES)
Write the coordinate arrays of the flat HEALPix store into the Zarr group `g`:

- `cell`: the flat index `1:npix` of every grid point, the store's horizontal dimension.
  HEALPix RING order, so cell `ij` is RING pixel `ij-1` in healpy/cuHPX's 0-based indexing.
- `lat`, `lon`: latitude (˚N) and longitude (˚E) of every cell, one flat vector each.
- `ring`: the 1-based latitude ring index of every cell, north to south.
- `layer`, `soil_layer`: the vertical coordinates, as for [`ZarrOutput`](@ref).

`lat`, `lon` and `ring` are all defined along `cell` and are therefore *non-dimension*
coordinates in the xarray sense: a consumer indexes into `cell` and reads off where on the
sphere that cell sits, with no knowledge of SpeedyWeather's grid machinery needed."""
function write_healpix_coordinates!(g::Zarr.ZGroup, output::HEALPixOutput, model::AbstractModel)
    grid = output.field2D.grid
    npix = get_npoints(grid)
    nside = RingGrids.nside_healpix(grid.nlat_half)

    # flat coordinates of every grid point, in ring order (0-360˚E, then north to south)
    londs, latds = get_londlatds(grid)
    rings = collect(whichring(grid))            # ring index j of every grid point ij

    σ = convert.(eltype(latds), on_architecture(CPU(), model.geometry.σ_levels_full))
    soil_indices = collect(1:get_soil_layers(model))

    write_coordinate!(
        g, "cell", collect(1:npix);
        attrs = Dict(
            "units" => "1", "_ARRAY_DIMENSIONS" => ["cell"],
            "long_name" => "HEALPix cell index (RING ordering, 1-based)",
            "nside" => nside, "ordering" => "RING",
        )
    )
    write_coordinate!(
        g, "lat", collect(latds);
        attrs = Dict(
            "units" => "degrees_north", "long_name" => "latitude",
            "standard_name" => "latitude", "_ARRAY_DIMENSIONS" => ["cell"],
        )
    )
    write_coordinate!(
        g, "lon", collect(londs);
        attrs = Dict(
            "units" => "degrees_east", "long_name" => "longitude",
            "standard_name" => "longitude", "_ARRAY_DIMENSIONS" => ["cell"],
        )
    )
    write_coordinate!(
        g, "ring", rings;
        attrs = Dict(
            "units" => "1", "_ARRAY_DIMENSIONS" => ["cell"],
            "long_name" => "latitude ring index of the cell, 1 (north) to $(get_nlat(grid)) (south)",
        )
    )
    write_coordinate!(
        g, "layer", collect(σ);
        attrs = Dict("units" => "1", "long_name" => "sigma layer", "_ARRAY_DIMENSIONS" => ["layer"])
    )
    write_coordinate!(
        g, "soil_layer", collect(soil_indices);
        attrs = Dict("units" => "1", "long_name" => "soil layer index", "_ARRAY_DIMENSIONS" => ["soil_layer"])
    )
    return nothing
end

"""$(TYPEDSIGNATURES)
Define a Zarr array for output `var` in the Zarr group `g` in the flat HEALPix layout:
`(cell, vertical, time)` with the dimensions `var` doesn't have collapsed away. The time
axis is pre-allocated to its final length `n_outputs`; unwritten chunks read back as the
array's `fill_value`."""
function define_variable!(
        g::Zarr.ZGroup,
        output::HEALPixOutput,
        var::AbstractOutputVariable,
        n_outputs::Int,
        output_NF::Type{<:AbstractFloat} = DEFAULT_OUTPUT_NF,
    )
    # hook for custom output variables to define their own (vertical) dimension
    define_dimension!(g, var)

    # x and y collapse into the single flat `cell` dimension, so a variable without a
    # horizontal dimension has no place in this layout
    @assert var.dims_xyzt[1] && var.dims_xyzt[2] "HEALPixOutput only supports output " *
        "variables with a horizontal dimension, but $(var.name) has dims_xyzt=$(var.dims_xyzt)."

    missing_value = hasfield(typeof(var), :missing_value) ? var.missing_value : DEFAULT_MISSING_VALUE

    ncells = get_npoints(output.field2D.grid)
    nz = get_nlayers(output, var)
    full_shape = (ncells, nz, n_outputs)

    # Chunking: 0 (default) or any non-positive value ⇒ full extent. Otherwise clamp to the
    # dimension size so users can't request chunks larger than the array (Zarr needs chunk ≤ shape).
    cc = output.cell_chunk > 0 ? min(output.cell_chunk, ncells) : ncells
    cz = output.vertical_chunk > 0 ? min(output.vertical_chunk, nz) : nz
    full_chunks = (cc, cz, max(output.time_chunk, 1))

    # the vertical dimension depends on the variable, e.g. "layer" or "soil_layer"
    all_dims = ("cell", vertical_dimension(var), "time")

    # pick out the active dims: cell is always on (asserted above), the vertical and time
    # ones follow the variable's z/t flags
    active = (true, is3D(var), hastime(var))
    shape = Tuple(d for (d, on) in zip(full_shape, active) if on)
    chunks = Tuple(c for (c, on) in zip(full_chunks, active) if on)
    dims = String[string(d) for (d, on) in zip(all_dims, active) if on]

    # `_ARRAY_DIMENSIONS` is read row-major by xarray, our shape/chunks are column-major
    reverse!(dims)

    return zcreate_output_variable!(g, output, var, shape, chunks, dims, output_NF, missing_value)
end
