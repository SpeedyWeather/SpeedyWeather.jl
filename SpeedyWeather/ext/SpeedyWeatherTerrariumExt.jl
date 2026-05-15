module SpeedyWeatherTerrariumExt

using SpeedyWeather
using Terrarium
using DocStringExtensions

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

"""$(TYPEDEF)

The canonical [`AbstractTerrariumLandModel`](@ref) implementation for wet
atmospheric models: wraps a full Terrarium land model (soil + surface energy
balance + vegetation) and all the metadata required to construct its initial
state. Pass as `land=` to `PrimitiveWetModel`.

$(TYPEDFIELDS)"""
struct TerrariumWetLand{
        NF,
        LG <: LandGeometry,
        TM <: Terrarium.AbstractModel{NF},
        TS <: Terrarium.AbstractTimeStepper,
        BC <: NamedTuple,
        IV <: Tuple,
        IN <: NamedTuple,
        FL <: NamedTuple,
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
end

"""$(TYPEDSIGNATURES)

Construct a [`TerrariumWetLand`](@ref) from a Terrarium land model and a
SpeedyWeather spectral grid. The single effective land layer thickness is
derived from the bottom-most layer of the Terrarium grid."""
function TerrariumWetLand(
        spectral_grid::SpectralGrid,
        model::Terrarium.AbstractModel{NF};
        timestepper::Terrarium.AbstractTimeStepper = ForwardEuler(NF),
        boundary_conditions::NamedTuple = (;),
        input_variables::Tuple = (),
        initializers::NamedTuple = (;),
        fields::NamedTuple = (;),
        Δt::Real = 300,
    ) where {NF}
    field_grid = Terrarium.get_field_grid(model.grid)
    Δz_arr = Terrarium.on_architecture(Terrarium.CPU(), field_grid.z.Δᵃᵃᶜ)
    # The Oceananigans vertical spacing is an OffsetArray; take the last entry
    # as the bottom-most layer thickness for the SpeedyWeather LandGeometry.
    geometry = LandGeometry(1, NF[Δz_arr[end]])
    return TerrariumWetLand(
        spectral_grid, geometry, model, timestepper,
        boundary_conditions, input_variables, initializers, fields, NF(Δt),
    )
end

"""$(TYPEDSIGNATURES)

Construct a [`TerrariumWetLand`](@ref) from a pre-built Terrarium
`ModelIntegrator`. The SpeedyWeather spectral grid is built from the Terrarium
ring grid; extra keyword arguments are forwarded to the `SpectralGrid`
constructor."""
function TerrariumWetLand(
        integrator::ModelIntegrator{NF, Arch, Grid};
        spectral_grid_kwargs...,
    ) where {NF, Arch, Grid <: ColumnRingGrid}
    spectral_grid = SpectralGrid(integrator.model.grid.rings; NF, spectral_grid_kwargs...)
    return TerrariumWetLand(
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
            name = :soil_temperature, dims = SpeedyWeather.Land2D(),
            units = "K", desc = "Soil temperature mirrored from Terrarium",
            namespace = :land,
        ),
        SpeedyWeather.PrognosticVariable(
            name = :soil_moisture, dims = SpeedyWeather.Land2D(),
            units = "1", desc = "Soil moisture (saturation fraction) mirrored from Terrarium",
            namespace = :land,
        ),
    )
end

function SpeedyWeather.initialize!(
        vars::Variables,
        land::AbstractTerrariumLandModel,
        ::PrimitiveEquation,
    )
    state = vars.prognostic.land.terrarium
    NF = eltype(vars.prognostic.land.soil_temperature)

    # Sync the Terrarium clock's initial time with the SpeedyWeather clock,
    # which was set from the `time` kwarg of `initialize!(model; time=...)`.
    state.clock.time = vars.prognostic.clock.time

    # initialize the "ModelIntegrator"
    integrator = ModelIntegrator(
        state.clock, land.model, InputSources(),
        state, land.initializers, land.timestepper,
    )
    Terrarium.initialize!(integrator)

    Tsoil = interior(state.temperature)[:, 1, end] .+ NF(273.15)
    sat = interior(state.saturation_water_ice)[:, 1, end]
    vars.prognostic.land.soil_temperature .= Tsoil
    vars.prognostic.land.soil_moisture .= sat
    return nothing
end

function SpeedyWeather.timestep!(
        vars::Variables,
        land::AbstractTerrariumLandModel,
        ::PrimitiveEquation,
    )
    state = vars.prognostic.land.terrarium
    tmodel = land.model
    consts = tmodel.constants
    NF = eltype(state)

    # Atmospheric forcings on the lowest model level / surface
    Tair = @view vars.grid.temperature[:, end]
    humid = @view vars.grid.humidity[:, end]
    pres = vars.grid.pressure                                # log surface pressure
    wind = vars.parameterizations.surface_wind_speed
    rain = vars.parameterizations.rain_rate
    snow = vars.parameterizations.snow_rate
    Rsd = vars.parameterizations.surface_shortwave_down
    Rld = vars.parameterizations.surface_longwave_down

    # Push forcings into Terrarium inputs (set! avoids allocating new fields)
    inputs = state.inputs
    Terrarium.set!(inputs.air_temperature, Tair)
    Terrarium.set!(inputs.air_temperature, inputs.air_temperature - NF(273.15))   # K -> °C
    Terrarium.set!(inputs.air_pressure, pres)
    Terrarium.set!(inputs.air_pressure, exp(inputs.air_pressure))                  # log(Pa) -> Pa
    Terrarium.set!(inputs.specific_humidity, humid)
    Terrarium.set!(inputs.rainfall, rain)
    Terrarium.set!(inputs.snowfall, snow)
    Terrarium.set!(inputs.windspeed, wind)
    Terrarium.set!(inputs.surface_shortwave_down, Rsd)
    Terrarium.set!(inputs.surface_longwave_down, Rld)

    # Constructing ModelIntegrator is allocation-free: it is an immutable struct
    # of references so no data is copied.  InputSources() is empty intentionally:
    # we own the input-update cycle above (via set!) and do not want Terrarium's
    # update_inputs! to overwrite those values during the substeps.
    integrator = ModelIntegrator(
        state.clock, tmodel, InputSources(),
        state, land.initializers, land.timestepper,
    )
    Terrarium.run!(integrator; period = vars.prognostic.clock.Δt, Δt = land.Δt)

    # Push Terrarium output back into SpeedyWeather coupling variables.
    # The Terrarium state itself is mutated in place above and remains in
    # `vars.prognostic.land.terrarium`; the soil mirrors are refreshed so
    # SpeedyWeather radiation / surface flux components see current values.
    vars.prognostic.land.soil_temperature .= state.skin_temperature .+ NF(273.15)
    vars.prognostic.land.soil_moisture .= interior(state.saturation_water_ice)[:, 1, end]
    if haskey(vars.prognostic.land, :sensible_heat_flux)
        vars.prognostic.land.sensible_heat_flux .= state.sensible_heat_flux
    end
    if haskey(vars.prognostic.land, :surface_humidity_flux)
        vars.prognostic.land.surface_humidity_flux .= state.latent_heat_flux ./ consts.thermodynamics.latent_heat_vaporization
    end
    if haskey(vars.parameterizations, :surface_longwave_up)
        vars.parameterizations.surface_longwave_up .= state.surface_longwave_up
    end
    if haskey(vars.parameterizations, :surface_shortwave_up)
        vars.parameterizations.surface_shortwave_up .= state.surface_shortwave_up
    end
    return nothing
end

"""$(TYPEDEF)

[`AbstractTerrariumLandModel`](@ref) for dry atmospheric models: wraps a
Terrarium soil model (heat conduction only, no surface energy balance or
moisture). The only atmospheric forcing passed in is near-surface air
temperature; soil moisture is not tracked. Pass as `land=` to
`PrimitiveDryModel`.

$(TYPEDFIELDS)"""
struct TerrariumDryLand{
        NF,
        LG <: LandGeometry,
        TM <: Terrarium.AbstractModel{NF},
        TS <: Terrarium.AbstractTimeStepper,
        BC <: NamedTuple,
        IV <: Tuple,
        IN <: NamedTuple,
        FL <: NamedTuple,
    } <: AbstractTerrariumLandModel
    "SpeedyWeather spectral grid"
    spectral_grid::SpectralGrid
    "SpeedyWeather land geometry (a single effective surface layer)"
    geometry::LG
    "Underlying Terrarium soil model"
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
end

"""$(TYPEDSIGNATURES)

Construct a [`TerrariumDryLand`](@ref) from a Terrarium soil model and a
SpeedyWeather spectral grid."""
function TerrariumDryLand(
        spectral_grid::SpectralGrid,
        model::Terrarium.AbstractModel{NF};
        timestepper::Terrarium.AbstractTimeStepper = ForwardEuler(NF),
        boundary_conditions::NamedTuple = (;),
        input_variables::Tuple = (),
        initializers::NamedTuple = (;),
        fields::NamedTuple = (;),
        Δt::Real = 300,
    ) where {NF}
    field_grid = Terrarium.get_field_grid(model.grid)
    Δz_arr = Terrarium.on_architecture(Terrarium.CPU(), field_grid.z.Δᵃᵃᶜ)
    geometry = LandGeometry(1, NF[Δz_arr[end]])
    return TerrariumDryLand(
        spectral_grid, geometry, model, timestepper,
        boundary_conditions, input_variables, initializers, fields, NF(Δt),
    )
end

"""$(TYPEDSIGNATURES)

Construct a [`TerrariumDryLand`](@ref) from a pre-built Terrarium
`ModelIntegrator`. The SpeedyWeather spectral grid is built from the Terrarium
ring grid; extra keyword arguments are forwarded to the `SpectralGrid`
constructor."""
function TerrariumDryLand(
        integrator::ModelIntegrator{NF, Arch, Grid};
        spectral_grid_kwargs...,
    ) where {NF, Arch, Grid <: ColumnRingGrid}
    spectral_grid = SpectralGrid(integrator.model.grid.rings; NF, spectral_grid_kwargs...)
    return TerrariumDryLand(
        spectral_grid, integrator.model;
        timestepper = integrator.timestepper,
        initializers = integrator.initializers,
    )
end

# Dry land: only soil_temperature mirror, no moisture.
function SpeedyWeather.variables(::TerrariumDryLand)
    return (
        SpeedyWeather.PrognosticVariable(
            name = :terrarium, dims = TerrariumVars(),
            namespace = :land, desc = "Terrarium land state",
        ),
        SpeedyWeather.PrognosticVariable(
            name = :soil_temperature, dims = SpeedyWeather.Land3D(),
            units = "K", desc = "Soil temperature mirrored from Terrarium",
            namespace = :land,
        ),
    )
end

function SpeedyWeather.initialize!(
        vars::Variables,
        land::TerrariumDryLand,
        ::PrimitiveEquation,
    )
    state = vars.prognostic.land.terrarium
    NF = eltype(vars.prognostic.land.soil_temperature)

    # Sync the Terrarium clock's initial time with the SpeedyWeather clock,
    # which was set from the `time` kwarg of `initialize!(model; time=...)`.
    state.clock.time = vars.prognostic.clock.time

    vars.prognostic.land.soil_temperature .= interior(state.temperature)[:, 1, end] .+ NF(273.15)
    return nothing
end

function SpeedyWeather.timestep!(
        vars::Variables,
        land::TerrariumDryLand,
        ::PrimitiveEquation,
    )
    state = vars.prognostic.land.terrarium
    tmodel = land.model
    NF = eltype(state)

    # Only air temperature is needed; convert K -> °C
    Tair = @view vars.grid.temperature[:, end]
    inputs = state.inputs
    Terrarium.set!(inputs.air_temperature, Tair)
    Terrarium.set!(inputs.air_temperature, inputs.air_temperature - NF(273.15))

    # Same reasoning as in TerrariumWetLand.timestep!: free to construct, empty
    # InputSources() so SpeedyWeather owns the input-update cycle.
    integrator = ModelIntegrator(
        state.clock, tmodel, InputSources(),
        state, land.initializers, land.timestepper,
    )
    Terrarium.run!(integrator; period = vars.prognostic.clock.Δt, Δt = land.Δt)

    # Surface soil temperature from the bottom of the column (last z-index)
    vars.prognostic.land.soil_temperature .= interior(state.temperature)[:, 1, end] .+ NF(273.15)
    return nothing
end

end # module
