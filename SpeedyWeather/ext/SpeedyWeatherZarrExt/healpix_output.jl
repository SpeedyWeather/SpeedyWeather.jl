# `HEALPixOutput`: writes output onto a `HEALPixGrid`or `OctaHEALPixGridkeeping the horizontal dimension flat,
# with arrays shaped (cell, vertical, time). Cells are in HEALPix RING order so the store can
# be read by healpy/cuHPX without reordering. Everything that is not specific to that layout —
# time-chunk buffering, the time axis, coordinate helpers — lives in shared.jl.

"""Grids `HEALPixOutput` can write, both equal-area and ring-ordered north to south."""
const OUTPUT_HEALPIX_GRIDS = Union{HEALPixGrid, OctaHEALPixGrid}

"""$(TYPEDSIGNATURES)
Resolution parameter `nlat_half` for an output grid of type `Grid`, from the user's
`nlat_half` keyword or the model grid's own `nlat_half`. `HEALPixGrid` is only defined for
even `nlat_half` (`nside = nlat_half ÷ 2`), so an odd value is rounded up; `OctaHEALPixGrid`
takes any positive `nlat_half`."""
function healpix_nlat_half(
        Grid::Type{<:OUTPUT_HEALPIX_GRIDS}, model_nlat_half::Integer; nlat_half = nothing
    )
    n = Int(something(nlat_half, model_nlat_half))
    n > 0 || throw(ArgumentError("HEALPixOutput: nlat_half=$n must be positive."))
    (Grid <: HEALPixGrid && isodd(n)) || return n

    @info "HEALPixOutput: HEALPixGrid is only defined for even nlat_half, rounding $n up to $(n + 1) (nside=$((n + 1) ÷ 2))."
    return n + 1
end

"""$(TYPEDSIGNATURES)
Resolve the output grid of a [`HEALPixOutput`](@ref) from the keyword arguments and the model
grid `model_grid`, always on CPU as we write from the host. In order of precedence:

1. an explicit `output_grid`, which must be a `HEALPixGrid` or an `OctaHEALPixGrid` and
   cannot be combined with `nside`/`nlat_half`;
2. `nside`, which is a HEALPix-specific parameter and therefore always selects a
   `HEALPixGrid` of `nlat_half = 2nside`;
3. otherwise the *model's own grid type* when the model already runs on one of the supported
   grids — so a HEALPix or OctaHEALPix simulation is written on its own grid and the
   interpolation is skipped — and a `HEALPixGrid` for every other model grid. The resolution
   is `nlat_half` if given, else the model's."""
function healpix_output_grid(model_grid::AbstractGrid; nside = nothing, nlat_half = nothing, output_grid = nothing)
    if !isnothing(output_grid)
        isnothing(nside) && isnothing(nlat_half) || throw(
            ArgumentError(
                "HEALPixOutput: pass either `output_grid` or `nside`/`nlat_half`, not both."
            )
        )
        return on_architecture(CPU(), output_grid)
    end

    if !isnothing(nside)
        isnothing(nlat_half) || throw(
            ArgumentError("HEALPixOutput: pass either `nside` or `nlat_half`, not both.")
        )
        return HEALPixGrid(2 * Int(nside), CPU())
    end

    # follow the model's grid type when it is already one we can write, so that a HEALPix or
    # OctaHEALPix simulation skips the interpolation by default
    Grid = model_grid isa OUTPUT_HEALPIX_GRIDS ? RingGrids.nonparametric_type(model_grid) : HEALPixGrid
    return Grid(healpix_nlat_half(Grid, model_grid.nlat_half; nlat_half), CPU())
end

"""$(TYPEDSIGNATURES)
Number of square faces the grid's sphere is tessellated into: 12 for the standard
`HEALPixGrid` (Górski et al. 2005, Nθ=3, Nφ=4), 4 for the `OctaHEALPixGrid` (Nθ=1, Nφ=4)."""
healpix_nfaces(::HEALPixGrid) = 12
healpix_nfaces(::OctaHEALPixGrid) = 4

"""$(TYPEDSIGNATURES)
Resolution parameter `nside` of an output grid: the side length, in grid points, of one of its
square faces, so that `npix = nfaces * nside^2` holds for both supported grids. For a standard
`HEALPixGrid` this is `nlat_half ÷ 2` and recovers the usual `npix = 12nside²`; for an
`OctaHEALPixGrid`, whose 4 faces are `nlat_half × nlat_half`, it is `nlat_half` itself and
gives `npix = 4nside²`."""
healpix_nside(grid::HEALPixGrid) = RingGrids.nside_healpix(grid.nlat_half)
healpix_nside(grid::OctaHEALPixGrid) = grid.nlat_half

"""$(TYPEDSIGNATURES)
Global Zarr attributes describing the output `grid`, so the store is self-describing.

Both supported grids carry the same key set, so a reader needs no special-casing:
`healpix_nside` / `healpix_npix` / `healpix_order` alongside the grid-independent `grid`,
`npix`, `nlat_half`, `nrings`, `cell_ordering` and `equal_area`. `nside` generalizes as the
side length of one square face (see [`healpix_nside`](@ref)) so that
`healpix_npix == healpix_nfaces * healpix_nside^2` for both, and `healpix_order` is `"RING"`
for both because cells run ring by ring north to south in either case.

**`healpix_nfaces` is the key that tells them apart**, and matters: a reader that assumes the
standard 12-face tessellation will compute `12nside²` and silently misread an `OctaHEALPixGrid`
store, whose `healpix_npix` is `4nside²`. `OctaHEALPixGrid` is a SpeedyWeather-specific member
of the HEALPix family, not the tessellation healpy and cuHPX implement, so its stores also
carry an explicit `note` saying so; check `grid` or `healpix_nfaces` before handing a store to
either tool."""
function healpix_store_attributes(grid::OUTPUT_HEALPIX_GRIDS)
    nside = healpix_nside(grid)
    nfaces = healpix_nfaces(grid)
    npix = get_npoints(grid)
    @assert npix == nfaces * nside^2 "$(RingGrids.nonparametric_type(grid)) with " *
        "nlat_half=$(grid.nlat_half): npix=$npix but nfaces*nside^2=$(nfaces * nside^2)."

    attrs = Dict{String, Any}(
        "grid" => string(nameof(RingGrids.nonparametric_type(grid))),
        "npix" => npix,
        "nlat_half" => grid.nlat_half,
        "nrings" => get_nlat(grid),
        "cell_ordering" => "ring-major, north to south, west to east within each ring",
        "equal_area" => true,
        "healpix_nside" => nside,
        "healpix_npix" => npix,
        "healpix_nfaces" => nfaces,
        "healpix_order" => "RING",
    )

    grid isa HEALPixGrid && return attrs

    attrs["note"] = "OctaHEALPixGrid is a SpeedyWeather-specific member of the HEALPix " *
        "family (4 faces, no equatorial belt, npix = 4*nside^2) and is NOT the standard " *
        "12-face HEALPix tessellation: healpy and cuHPX cannot read it. The store layout, " *
        "coordinates and attribute names follow the same convention as a HEALPixGrid store."
    return attrs
end

"""$(TYPEDSIGNATURES)
Constructor for [`HEALPixOutput`](@ref). Builds the output `HEALPixGrid` and the scratch fields 
to write from, and pre-fills `output.variables` with the defaults for `Model`.

The output grid is a `HEALPixGrid` or an `OctaHEALPixGrid`; both are written with the same
flat, ring-ordered store convention, but only a `HEALPixGrid` store is standard HEALPix and
readable by healpy/cuHPX (see [`healpix_store_attributes`](@ref)).

Its resolution is given either as `nside` (the HEALPix parameter, as used e.g. by cuHPX; it
always selects a `HEALPixGrid` of `nlat_half = 2nside`), as `nlat_half`, or as a ready-made
`output_grid`. By default the output follows the *model's* grid: a HEALPix or OctaHEALPix
simulation is written on its own grid at its own resolution — so the interpolation is skipped
— and any other model grid is written on a `HEALPixGrid` of the model's `nlat_half`, rounded
up to the nearest even number as `HEALPixGrid` requires.

An interpolator onto that grid is only built when it is actually needed: if the model already
runs on the very same grid, `output.interpolator` stays `nothing` and the data is copied
straight out, skipping interpolation entirely."""
function HEALPixOutput(
        SG::SpectralGrid,
        Model::Type{<:AbstractModel} = Barotropic;
        nlayers_soil = DEFAULT_NLAYERS_SOIL,
        nside = nothing,
        nlat_half = nothing,
        output_grid::Union{OUTPUT_HEALPIX_GRIDS, Nothing} = nothing,
        output_NF::DataType = DEFAULT_OUTPUT_NF,
        interval::Period = Second(DEFAULT_OUTPUT_INTERVAL),
        compressor = nothing,
        kwargs...
    )

    # OUTPUT GRID: from `output_grid`, from `nside`/`nlat_half`, or following the model's own
    # grid type and resolution (which is what makes the interpolation skippable by default)
    output_grid = healpix_output_grid(SG.grid; nside, nlat_half, output_grid)

    # INPUT GRID (but on CPU)
    input_grid = on_architecture(CPU(), SG.grid)

    # SKIP INTERPOLATION if the model already runs on this very grid: `grids_match`
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

    # global attributes so the store is self-describing; a HEALPixGrid store additionally
    # advertises the healpix_* keys healpy/cuHPX look for, an OctaHEALPixGrid one must not
    g = Zarr.zgroup(store_path; attrs = healpix_store_attributes(output.field2D.grid))
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
Write the coordinate arrays of the flat store into the Zarr group `g`:

- `cell`: the flat index `1:npix` of every grid point, the store's horizontal dimension,
  running ring by ring north to south. On a `HEALPixGrid` that is the standard HEALPix RING
  order, so cell `ij` is RING pixel `ij-1` in healpy/cuHPX's 0-based indexing.
- `lat`, `lon`: latitude (˚N, from the equator) and longitude (˚E, in `[0, 360)`) of every
  cell, one flat vector each. This matches healpy's `pix2ang(..., lonlat=true)` convention;
  healpy's *default* is colatitude in radians, i.e. `θ = deg2rad(90 - lat)`.
- `ring`: the 1-based latitude ring index of every cell, north to south.
- `layer`, `soil_layer`: the vertical coordinates, as for [`ZarrOutput`](@ref).

`lat`, `lon` and `ring` are all defined along `cell` and are therefore *non-dimension*
coordinates in the xarray sense: a consumer indexes into `cell` and reads off where on the
sphere that cell sits, with no knowledge of SpeedyWeather's grid machinery needed. This is
the same for both supported grids; only the tessellation behind the coordinates differs."""
function write_healpix_coordinates!(g::Zarr.ZGroup, output::HEALPixOutput, model::AbstractModel)
    grid = output.field2D.grid
    npix = get_npoints(grid)
    grid_name = string(nameof(RingGrids.nonparametric_type(grid)))   # bare name, see above

    # flat coordinates of every grid point, in ring order (0-360˚E, then north to south)
    londs, latds = get_londlatds(grid)
    rings = collect(whichring(grid))            # ring index j of every grid point ij

    σ = convert.(eltype(latds), on_architecture(CPU(), model.geometry.σ_levels_full))
    soil_indices = collect(1:get_soil_layers(model))

    write_coordinate!(
        g, "cell", collect(1:npix);
        attrs = Dict{String, Any}(
            "units" => "1", "_ARRAY_DIMENSIONS" => ["cell"],
            "long_name" => "$grid_name cell index, ring-major north to south, 1-based",
            "ordering" => "RING",
            "nside" => healpix_nside(grid),
            "nfaces" => healpix_nfaces(grid),
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
