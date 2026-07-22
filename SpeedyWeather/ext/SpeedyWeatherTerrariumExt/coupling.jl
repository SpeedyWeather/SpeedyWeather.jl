
"""
$(TYPEDSIGNATURES)

Construct a [`TerrariumLand`](@ref) from a pre-built Terrarium model for usage 
as a SpeedyWeather land model. Based on the SpeedyWeather model type either a 
dry or wet land model is used."""
SpeedyWeather.LandModel(
    spectral_grid::SpectralGrid,
    model::Terrarium.AbstractModel{NF}; kwargs...
) where {NF} =
    TerrariumLand(spectral_grid, model; kwargs...)

"""$(TYPEDEF)

Variable dimension for the full Terrarium land state. `allocate` for an
`AbstractVariable{TerrariumVars}` calls `Terrarium.initialize` on `model.land`
so that the entire Terrarium state lives inside SpeedyWeather's `Variables`
container at `vars.prognostic.land.terrarium`. The `model.land` must therefore
be an [`AbstractTerrariumLandModel`](@ref)."""
struct TerrariumVars <: SpeedyWeather.AbstractVariableDim end

function SpeedyWeather.allocate(::SpeedyWeather.AbstractVariable{TerrariumVars}, model::AbstractModel)
    land = model.land
    # Construct a Terrarium clock with a DateTime placeholder so its time type is
    # DateTime; the actual initial datetime is synced from the SpeedyWeather clock
    # in `initialize!(vars, land, model)` once the user-provided `time` has been
    # written there.
    return Terrarium.initialize(
        land.model;
        clock = Terrarium.Clock(time = SpeedyWeather.DEFAULT_DATE),
        boundary_conditions = land.boundary_conditions,
        input_variables = land.input_variables,
        fields = land.fields,
    )
end

"""$(TYPEDEF)

Abstract base type for SpeedyWeather land components that delegate land
surface dynamics to a Terrarium land model. Subtypes must provide a
Terrarium model + integrator configuration and a SpeedyWeather
[`LandGeometry`](@ref). The default Terrarium-coupling interface is defined
on this abstract type:

* [`variables`](@ref) declares `vars.prognostic.land.terrarium` (the full
  Terrarium state) plus thin grid-side mirrors `:soil_temperature` and
  `:soil_moisture` that the rest of SpeedyWeather (radiation, surface
  fluxes, output) reads from.
* [`initialize!`](@ref) initializes the Terrarium timestepper against the
  freshly-allocated state and seeds the soil mirrors.
* [`timestep!`](@ref) pushes SpeedyWeather atmospheric forcings into the
  Terrarium inputs, runs the Terrarium integrator for the duration of the
  SpeedyWeather step using `land.Δt` as the sub-step, and
  copies soil temperatures, moisture and surface fluxes back into the
  SpeedyWeather variables.
"""
abstract type AbstractTerrariumLandModel <: SpeedyWeather.AbstractLand end

@inline SpeedyWeather.get_nlayers(land::AbstractTerrariumLandModel) = land.geometry.nlayers

"""$(TYPEDSIGNATURES)

Boolean land mask of the underlying Terrarium [`ColumnRingGrid`](@ref). Terrarium
only allocates columns where this mask is `true`, so every Terrarium state field
has length `sum(mask)` (the number of land columns)."""
@inline land_mask(land::AbstractTerrariumLandModel) = land.model.grid.mask.data

"""$(TYPEDSIGNATURES)

Boolean land mask (`Field{Bool}`) derived from a SpeedyWeather [`AbstractLandSeaMask`](@ref):
`true` wherever the fractional `land_fraction` exceeds `threshold`. Pass the result to
`Terrarium.ColumnRingGrid` so the Terrarium land columns are derived from the *same*
land-sea mask object the atmospheric model uses, instead of a separately loaded copy. The
default `threshold = 0` makes the Terrarium mask a superset of the SpeedyWeather land, which
is required for consistent coupling (see [`check_mask_consistency`](@ref))."""
land_mask(land_sea_mask::SpeedyWeather.AbstractLandSeaMask; threshold = 0) =
    land_sea_mask.land_fraction .> threshold

"""$(TYPEDSIGNATURES)

Convenience constructor building a Terrarium [`ColumnRingGrid`](@ref) directly from a
SpeedyWeather [`AbstractLandSeaMask`](@ref), so the Terrarium land columns and the
atmospheric land-sea mask share a single source. Equivalent to passing
`land_mask(land_sea_mask; threshold)` as the mask; `threshold` selects the fractional
land cutoff (default `0`, i.e. any land)."""
Terrarium.ColumnRingGrid(
    arch::Terrarium.AbstractArchitecture,
    NF::Type{<:AbstractFloat},
    vert::Terrarium.AbstractVerticalSpacing,
    rings::SpeedyWeather.RingGrids.AbstractGrid,
    land_sea_mask::SpeedyWeather.AbstractLandSeaMask;
    threshold = 0,
) = Terrarium.ColumnRingGrid(arch, NF, vert, rings, land_mask(land_sea_mask; threshold))

"""$(TYPEDSIGNATURES)

Fill the SpeedyWeather soil mirror variables (`soil_temperature`, and `soil_moisture`
if present) with fallback values at every grid point *outside* the Terrarium land mask.
Terrarium only owns the mask points; every other point keeps its allocation value (`0`)
otherwise, which shows up as unphysical 0 K soil temperatures at ocean/coastal cells and
smears into the output."""
function fill_fallback!(vars::Variables, land::AbstractTerrariumLandModel)
    mask = land_mask(land)
    NF = eltype(vars.prognostic.land.soil_temperature)
    st = vars.prognostic.land.soil_temperature.data
    @views st[.!mask, :] .= NF(land.ocean_temperature)
    if haskey(vars.prognostic.land, :soil_moisture)
        sm = vars.prognostic.land.soil_moisture.data
        @views sm[.!mask, :] .= NF(land.ocean_moisture)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)

Warn if the SpeedyWeather land-sea mask and the Terrarium land mask disagree. Terrarium
allocates a soil column only where its Boolean mask is `true`, whereas SpeedyWeather weights
land fluxes by the fractional `land_fraction`. A grid point with `land_fraction > 0` but no
Terrarium column receives no soil state from Terrarium and hence contributes no
land surface fluxes even though it is (partially) land. A column at
a point with `land_fraction == 0` only wastes compute. The recommended construction is
`mask = land_fraction .> 0` so the Terrarium mask is a superset of the SpeedyWeather land."""
function check_mask_consistency(land::AbstractTerrariumLandModel, land_sea_mask)
    mask = SpeedyWeather.on_architecture(SpeedyWeather.CPU(), land_mask(land))
    land_fraction = SpeedyWeather.on_architecture(SpeedyWeather.CPU(), land_sea_mask.land_fraction.data)

    is_speedy_land = land_fraction .> 0
    orphaned = is_speedy_land .& .!mask         # land in SpeedyWeather, no Terrarium column
    if any(orphaned)
        @warn "Terrarium land mask misses $(count(orphaned)) grid point(s) with " *
            "land_fraction > 0 (max land_fraction = $(maximum(land_fraction[orphaned]))). " *
            "These cells get fallback soil values and contribute zero land surface fluxes. " *
            "Build the Terrarium mask as `land_sea_mask.land_fraction .> 0` to cover all land."
    end

    superfluous = .!is_speedy_land .& mask      # Terrarium column over pure ocean
    if any(superfluous)
        @info "Terrarium allocates $(count(superfluous)) soil column(s) over grid points " *
            "with land_fraction == 0. These are integrated but never coupled back (wasted compute)."
    end
    return nothing
end

# extend the allocation free copies here to work with the Oceananigans Fields as well
@inline SpeedyWeather.RingGrids.copy_unmasked!(
    dest::Terrarium.Oceananigans.AbstractField,
    src::SpeedyWeather.AbstractField,
    indices) = SpeedyWeather.RingGrids.copy_unmasked!(interior(dest), src, indices)

@inline SpeedyWeather.RingGrids.copy_unmasked!(
    dest::SpeedyWeather.AbstractField, 
    src::Terrarium.Oceananigans.AbstractField,
    indices) = SpeedyWeather.RingGrids.copy_unmasked!(dest, interior(src), indices)

"""$(TYPEDEF)

The canonical [`AbstractTerrariumLandModel`](@ref) implementation for wet
atmospheric models: wraps a full Terrarium land model (soil + surface energy
balance + vegetation) and all the metadata required to construct its initial
state. Pass as `land=` to `PrimitiveWetModel`.

$(TYPEDFIELDS)"""
struct TerrariumLand{
        NF,
        LG <: LandGeometry,
        TM <: Terrarium.AbstractModel{NF},
        TS <: Terrarium.AbstractTimeStepper,
        BC <: NamedTuple,
        IV <: Tuple,
        IN <: NamedTuple,
        FL <: NamedTuple,
        MI,
    } <: AbstractTerrariumLandModel
    "SpeedyWeather spectral grid"
    spectral_grid::SpectralGrid
    "SpeedyWeather land geometry (a single effective surface layer)"
    geometry::LG
    "Underlying Terrarium land model"
    model::TM
    "Terrarium time stepper used inside each SpeedyWeather step"
    timestepper::TS
    "Boundary conditions forwarded to `Terrarium.initialize`"
    boundary_conditions::BC
    "Additional input variables forwarded to `Terrarium.initialize`"
    input_variables::IV
    "Field initializers forwarded to the on-the-fly `ModelIntegrator`"
    initializers::IN
    "Preconstructed Terrarium fields forwarded to `Terrarium.initialize`"
    fields::FL
    "Terrarium-internal sub-step (seconds) used to integrate within each SpeedyWeather step"
    Δt::NF
    "Indices of the common land sea mask, used for allocation-free copying between SpeedyWeather and Terrarium"
    mask_indices::MI
    "Fallback soil temperature [K] for grid points outside the Terrarium land mask (ocean-only cells)"
    ocean_temperature::NF
    "Fallback soil moisture (saturation fraction) [1] for grid points outside the Terrarium land mask (ocean-only cells)"
    ocean_moisture::NF
end

"""$(TYPEDSIGNATURES)

Construct a [`TerrariumLand`](@ref) from a Terrarium land model and a
SpeedyWeather spectral grid."""
function TerrariumLand(
        spectral_grid::SpectralGrid,
        model::Terrarium.AbstractModel{NF};
        timestepper::Terrarium.AbstractTimeStepper = ForwardEuler(NF),
        boundary_conditions::NamedTuple = (;),
        input_variables::Tuple = (),
        initializers::NamedTuple = (;),
        fields::NamedTuple = (;),
        Δt::Real = 300,
        ocean_temperature::Real = 285,
        ocean_moisture::Real = 0,
    ) where {NF}
    field_grid = Terrarium.get_field_grid(model.grid)
    Δz_arr = Terrarium.on_architecture(Terrarium.CPU(), field_grid.z.Δᵃᵃᶜ)
    # The Oceananigans vertical spacing is an OffsetArray; take the last entry
    # as the layer thickness for the SpeedyWeather LandGeometry. It's not actually used, but
    # we set it here for consistency.
    geometry = LandGeometry(1, NF[Δz_arr[end]])
    # `unmasked_indices` treats `true` as masked-out; `model.grid.mask` is `true` at land
    # points, so invert it to get the indices of the (unmasked) land columns.
    mask_indices = RingGrids.unmasked_indices(.!model.grid.mask)
    return TerrariumLand(
        spectral_grid, geometry, model, timestepper,
        boundary_conditions, input_variables, initializers, fields, NF(Δt), mask_indices,
        NF(ocean_temperature), NF(ocean_moisture),
    )
end

"""$(TYPEDSIGNATURES)

Construct a [`TerrariumLand`](@ref) from a pre-built Terrarium
`ModelIntegrator`. The SpeedyWeather spectral grid is built from the Terrarium
ring grid; extra keyword arguments are forwarded to the `SpectralGrid`
constructor."""
function TerrariumLand(
        integrator::ModelIntegrator{NF, Arch, Grid};
        spectral_grid_kwargs...,
    ) where {NF, Arch, Grid <: ColumnRingGrid}
    spectral_grid = SpectralGrid(integrator.model.grid.rings; NF, spectral_grid_kwargs...)
    return TerrariumLand(
        spectral_grid, integrator.model;
        timestepper = integrator.timestepper,
        initializers = integrator.initializers,
    )
end

function SpeedyWeather.variables(::AbstractTerrariumLandModel)
    return (
        # The full Terrarium state, owned by SpeedyWeather's Variables tree.
        SpeedyWeather.PrognosticVariable(
            name = :terrarium, dims = TerrariumVars(),
            namespace = :land, desc = "Terrarium land state",
        ),
        # Thin grid-side mirrors of surface soil temperature / moisture so the
        # rest of SpeedyWeather (longwave/shortwave radiation, surface fluxes,
        # output writers) can read them. They are kept in sync from the
        # Terrarium state inside `initialize!` and `timestep!`.
        SpeedyWeather.PrognosticVariable(
            name = :soil_temperature, dims = SpeedyWeather.Grid2D(),
            units = "K", desc = "Soil temperature mirrored from Terrarium",
            namespace = :land,
        ),
        SpeedyWeather.PrognosticVariable(
            name = :soil_moisture, dims = SpeedyWeather.Grid2D(),
            units = "1", desc = "Soil moisture (saturation fraction) mirrored from Terrarium",
            namespace = :land,
        ),
    )
end

# wet land model
function SpeedyWeather.initialize!(
        vars::Variables,
        land::AbstractTerrariumLandModel,
        model::PrimitiveWetModel,
    )
    state = vars.prognostic.land.terrarium
    NF = eltype(vars.prognostic.land.soil_temperature)
    mask = land_mask(land)

    # Sync the Terrarium clock's initial time with the SpeedyWeather clock,
    # which was set from the `time` kwarg of `initialize!(model; time=...)`.
    state.clock.time = vars.prognostic.clock.time

    # initialize the "ModelIntegrator"
    integrator = ModelIntegrator(
        state.clock, land.model, InputSources(NF),
        state, land.initializers, land.timestepper,
    )
    Terrarium.initialize!(integrator)

    Tsoil = interior(state.temperature)[:, 1, end] .+ NF(273.15)
    sat = interior(state.saturation_water_ice)[:, 1, end]

    @assert length(mask) == length(vars.prognostic.land.soil_temperature) "Terrarium land mask (length = $(length(mask))) does not span the full SpeedyWeather ring grid (length = $(length(vars.prognostic.land.soil_temperature)))."
    @assert count(mask) == length(Tsoil) "Number of Terrarium land columns (length = $(length(Tsoil))) does not match the number of land points in the mask (count = $(count(mask)))."

    # warn if the Terrarium mask and the SpeedyWeather land-sea mask disagree
    check_mask_consistency(land, model.land_sea_mask)

    # fill ocean/non-Terrarium points with fallback values first, then seed the land columns
    fill_fallback!(vars, land)
    vars.prognostic.land.soil_temperature[mask] .= Tsoil
    vars.prognostic.land.soil_moisture[mask] .= sat
    return nothing
end

function SpeedyWeather.timestep!(
        vars::Variables,
        land::AbstractTerrariumLandModel,
        model::PrimitiveWetModel,
    )
    state = vars.prognostic.land.terrarium
    tmodel = land.model
    consts = tmodel.constants
    NF = eltype(state)
    mask = land_mask(land)
    indices = land.mask_indices
    nlayers = model.spectral_grid.nlayers

    # Atmospheric forcings on the lowest model level / surface
    # choose step dimension depending on atmospheric time stepper
    # and read like parameterization via DummyParameterization
    l = SpeedyWeather.which_prognostic_step(vars.grid.temperature, model.time_stepping, SpeedyWeather.DummyParameterization())
    Tair = RingGrids.field_view(vars.grid.temperature, :, nlayers, l)
    humid = RingGrids.field_view(vars.grid.humidity, :, nlayers, l)
    pres = vars.parameterizations.surface_pressure
    wind = vars.parameterizations.surface_wind_speed
    rain = vars.parameterizations.rain_rate
    snow = vars.parameterizations.snow_rate
    Rsd = vars.parameterizations.surface_shortwave_down
    Rld = vars.parameterizations.surface_longwave_down

    # Push forcings into Terrarium inputs (copy_unmasked! avoids allocations)
    inputs = state.inputs
    RingGrids.copy_unmasked!(inputs.air_temperature, Tair, indices)
    Terrarium.set!(inputs.air_temperature, inputs.air_temperature - NF(273.15))   # K -> °C
    RingGrids.copy_unmasked!(inputs.air_pressure, pres, indices)
    RingGrids.copy_unmasked!(inputs.specific_humidity, humid, indices)
    RingGrids.copy_unmasked!(inputs.rainfall, rain, indices)
    RingGrids.copy_unmasked!(inputs.snowfall, snow, indices)
    RingGrids.copy_unmasked!(inputs.windspeed, wind, indices)
    RingGrids.copy_unmasked!(inputs.surface_shortwave_down, Rsd, indices)
    RingGrids.copy_unmasked!(inputs.surface_longwave_down, Rld, indices)

    # Constructing ModelIntegrator is allocation-free: it is an immutable struct
    # of references so no data is copied.  `InputSources` is empty intentionally:
    # we own the input-update cycle above (via set!) and do not want Terrarium's
    # update_inputs! to overwrite those values during the substeps.
    integrator = ModelIntegrator(
        state.clock, tmodel, InputSources(NF),
        state, land.initializers, land.timestepper,
    )
    Terrarium.run!(integrator; period = vars.prognostic.clock.Δt, Δt = land.Δt)

    # Push Terrarium output back into SpeedyWeather coupling variables.
    # The Terrarium state itself is mutated in place above and remains in
    # `vars.prognostic.land.terrarium`; the soil mirrors are refreshed so
    # SpeedyWeather radiation / surface flux components see current values.
    vars.prognostic.land.soil_temperature[mask] .= interior(state.skin_temperature) .+ NF(273.15)
    vars.prognostic.land.soil_moisture[mask] .= @view interior(state.saturation_water_ice)[:, 1, end]
    if haskey(vars.prognostic.land, :sensible_heat_flux)
        RingGrids.copy_unmasked!(vars.prognostic.land.sensible_heat_flux, state.sensible_heat_flux, indices)
    end
    if haskey(vars.prognostic.land, :surface_humidity_flux)
        RingGrids.copy_unmasked!(vars.prognostic.land.surface_humidity_flux, state.latent_heat_flux, indices)
        vars.prognostic.land.surface_humidity_flux.data ./= consts.thermodynamics.latent_heat_vaporization
    end
    if haskey(vars.parameterizations, :surface_longwave_up)
        RingGrids.copy_unmasked!(vars.parameterizations.surface_longwave_up, state.surface_longwave_up, indices)
    end
    if haskey(vars.parameterizations, :surface_shortwave_up)
        RingGrids.copy_unmasked!(vars.parameterizations.surface_shortwave_up, state.surface_shortwave_up, indices)
    end
    return nothing
end

# Dry land model
function SpeedyWeather.initialize!(
        vars::Variables,
        land::TerrariumLand,
        model::PrimitiveDryModel,
    )
    state = vars.prognostic.land.terrarium
    NF = eltype(vars.prognostic.land.soil_temperature)
    mask = land_mask(land)

    # Sync the Terrarium clock's initial time with the SpeedyWeather clock,
    # which was set from the `time` kwarg of `initialize!(model; time=...)`.
    state.clock.time = vars.prognostic.clock.time

    # warn if the Terrarium mask and the SpeedyWeather land-sea mask disagree
    check_mask_consistency(land, model.land_sea_mask)

    # fill ocean/non-Terrarium points with fallback values first, then seed the land columns
    fill_fallback!(vars, land)
    vars.prognostic.land.soil_temperature[mask] .= @view(interior(state.temperature)[:, 1, end]) .+ NF(273.15)
    return nothing
end

function SpeedyWeather.timestep!(
        vars::Variables,
        land::TerrariumLand,
        model::PrimitiveDryModel,
    )
    state = vars.prognostic.land.terrarium
    land_model = land.model
    NF = eltype(state)
    mask = land_mask(land)

    # Only air temperature is needed; convert K -> °C
    # choose step dimension depending on atmospheric time stepper
    # and read like parameterization via DummyParameterization
    l = SpeedyWeather.which_prognostic_step(vars.grid.temperature, model.time_stepping, SpeedyWeather.DummyParameterization())
    Tair = vars.grid.temperature[mask, end, l]
    inputs = state.inputs
    Terrarium.set!(inputs.air_temperature, Tair)
    Terrarium.set!(inputs.air_temperature, inputs.air_temperature - NF(273.15))

    # Same reasoning as in TerrariumLand.timestep!: free to construct, empty
    # InputSources(NF) so SpeedyWeather owns the input-update cycle.
    integrator = ModelIntegrator(
        state.clock, land_model, InputSources(NF),
        state, land.initializers, land.timestepper,
    )
    Terrarium.run!(integrator; period = vars.prognostic.clock.Δt, Δt = land.Δt)

    # Surface soil temperature from the bottom of the column (last z-index)
    vars.prognostic.land.soil_temperature[mask] .= @view(interior(state.temperature)[:, 1, end]) .+ NF(273.15)
    return nothing
end