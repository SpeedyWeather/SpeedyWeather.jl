module SpeedyWeatherTerrariumExt

using SpeedyWeather
using Terrarium
using DocStringExtensions

import SpeedyWeather: AbstractVariableDim, AbstractVariable, AbstractModel,
    AbstractLand, AbstractWetLand, SpectralGrid, LandGeometry,
    PrognosticVariable, Land3D, Variables,
    PrimitiveEquation, get_nlayers, model_class, model_type,
    variables, initialize!, timestep!
import SpeedyWeather.RingGrids

import Terrarium: Clock, ColumnRingGrid, ForwardEuler, InputSources,
    ModelIntegrator, interior, on_architecture, set!

"""$(TYPEDEF)

Variable dimension for the full Terrarium land state. `Base.zero` for an
`AbstractVariable{TerrariumVars}` calls `Terrarium.initialize` on `model.land`
so that the entire Terrarium state lives inside SpeedyWeather's `Variables`
container at `vars.prognostic.land.terrarium`. The `model.land` must therefore
be an [`AbstractTerrariumLandModel`](@ref)."""
struct TerrariumVars <: AbstractVariableDim end

function Base.zero(::AbstractVariable{TerrariumVars}, model::AbstractModel)
    land = model.land
    return Terrarium.initialize(
        terrarium_model(land);
        clock = terrarium_clock(land),
        boundary_conditions = terrarium_boundary_conditions(land),
        input_variables = terrarium_input_variables(land),
        fields = terrarium_fields(land),
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
  SpeedyWeather step using `terrarium_substep(land)` as the sub-step, and
  copies soil temperatures, moisture and surface fluxes back into the
  SpeedyWeather variables.

Subtypes are expected to expose the following accessors (default
implementations forward to fields of the same name on `land`):

* [`terrarium_model`](@ref)`(land)` → `Terrarium.AbstractModel`
* [`terrarium_timestepper`](@ref)`(land)`
* [`terrarium_clock`](@ref)`(land)` → `Terrarium.Clock`
* [`terrarium_boundary_conditions`](@ref)`(land)` → `NamedTuple`
* [`terrarium_input_variables`](@ref)`(land)` → `Tuple`
* [`terrarium_initializers`](@ref)`(land)` → `NamedTuple`
* [`terrarium_fields`](@ref)`(land)` → `NamedTuple`
* [`terrarium_substep`](@ref)`(land)` → seconds (`Real`)
* [`set_terrarium_timestepper!`](@ref)`(land, ts)` (only needed if the
  timestepper requires state-dependent initialization)
"""
abstract type AbstractTerrariumLandModel <: AbstractLand end

@inline get_nlayers(land::AbstractTerrariumLandModel) = land.geometry.nlayers

"""$(TYPEDSIGNATURES) Return the underlying Terrarium `AbstractModel`."""
terrarium_model(land::AbstractTerrariumLandModel) = land.model

"""$(TYPEDSIGNATURES) Return the Terrarium time stepper."""
terrarium_timestepper(land::AbstractTerrariumLandModel) = land.timestepper

"""$(TYPEDSIGNATURES) Return the Terrarium [`Clock`](@ref Terrarium.Clock)."""
terrarium_clock(land::AbstractTerrariumLandModel) = land.clock

"""$(TYPEDSIGNATURES) Return the Terrarium boundary conditions `NamedTuple`."""
terrarium_boundary_conditions(land::AbstractTerrariumLandModel) = land.boundary_conditions

"""$(TYPEDSIGNATURES) Return the tuple of additional Terrarium input variables."""
terrarium_input_variables(land::AbstractTerrariumLandModel) = land.input_variables

"""$(TYPEDSIGNATURES) Return the Terrarium field initializers `NamedTuple`."""
terrarium_initializers(land::AbstractTerrariumLandModel) = land.initializers

"""$(TYPEDSIGNATURES) Return the preconstructed Terrarium fields `NamedTuple`."""
terrarium_fields(land::AbstractTerrariumLandModel) = land.fields

"""$(TYPEDSIGNATURES) Return the Terrarium-internal sub-step in seconds."""
terrarium_substep(land::AbstractTerrariumLandModel) = land.Δt

"""$(TYPEDSIGNATURES) Store a (possibly state-initialized) Terrarium time stepper on `land`."""
set_terrarium_timestepper!(land::AbstractTerrariumLandModel, ts) = setfield!(land, :timestepper, ts)

"""$(TYPEDEF)

The canonical [`AbstractTerrariumLandModel`](@ref) implementation: a
SpeedyWeather wet-land component that wraps a Terrarium land model and all
of the metadata required to construct its initial state. A
`TerrariumWetLand` is passed as the `land=` keyword argument of a
`PrimitiveWetModel`.

$(TYPEDFIELDS)"""
mutable struct TerrariumWetLand{
        NF,
        LG <: LandGeometry,
        TM <: Terrarium.AbstractModel{NF},
        TS <: Terrarium.AbstractTimeStepper,
        CK <: Clock,
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
    "Initial Terrarium clock"
    clock::CK
    "Boundary conditions forwarded to `Terrarium.initialize`"
    boundary_conditions::BC
    "Additional input variables forwarded to `Terrarium.initialize`"
    input_variables::IV
    "Field initializers forwarded to the on-the-fly `ModelIntegrator`"
    initializers::IN
    "Preconstructed Terrarium fields forwarded to `Terrarium.initialize`"
    fields::FL
    "Terrarium-internal sub-step (seconds) used to integrate within each SpeedyWeather step"
    Δt::Float64
end

"""$(TYPEDSIGNATURES)

Construct a [`TerrariumWetLand`](@ref) from a Terrarium land model and a
SpeedyWeather spectral grid. The single effective land layer thickness is
derived from the bottom-most layer of the Terrarium grid."""
function TerrariumWetLand(
        spectral_grid::SpectralGrid,
        model::Terrarium.AbstractModel{NF};
        timestepper::Terrarium.AbstractTimeStepper = ForwardEuler(NF),
        clock::Clock = Clock(time = zero(NF)),
        boundary_conditions::NamedTuple = (;),
        input_variables::Tuple = (),
        initializers::NamedTuple = (;),
        fields::NamedTuple = (;),
        Δt::Real = 300,
    ) where {NF}
    field_grid = Terrarium.get_field_grid(model.grid)
    Δz_arr = on_architecture(Terrarium.CPU(), field_grid.z.Δᵃᵃᶜ)
    # The Oceananigans vertical spacing is an OffsetArray; take the last entry
    # as the bottom-most layer thickness for the SpeedyWeather LandGeometry.
    geometry = LandGeometry(1, NF[Δz_arr[end]])
    return TerrariumWetLand(
        spectral_grid, geometry, model, timestepper, clock,
        boundary_conditions, input_variables, initializers, fields, Float64(Δt),
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
        clock = integrator.clock,
        initializers = integrator.initializers,
    )
end

function variables(::AbstractTerrariumLandModel)
    return (
        # The full Terrarium state, owned by SpeedyWeather's Variables tree.
        PrognosticVariable(
            name = :terrarium, dims = TerrariumVars(),
            namespace = :land, desc = "Terrarium land state",
        ),
        # Thin grid-side mirrors of surface soil temperature / moisture so the
        # rest of SpeedyWeather (longwave/shortwave radiation, surface fluxes,
        # output writers) can read them. They are kept in sync from the
        # Terrarium state inside `initialize!` and `timestep!`.
        PrognosticVariable(
            name = :soil_temperature, dims = Land3D(),
            units = "K", desc = "Soil temperature mirrored from Terrarium",
            namespace = :land,
        ),
        PrognosticVariable(
            name = :soil_moisture, dims = Land3D(),
            units = "1", desc = "Soil moisture (saturation fraction) mirrored from Terrarium",
            namespace = :land,
        ),
    )
end

function initialize!(
        vars::Variables,
        land::AbstractTerrariumLandModel,
        ::PrimitiveEquation,
    )
    state = vars.prognostic.land.terrarium
    NF = eltype(vars.prognostic.land.soil_temperature)
    # Terrarium time steppers may need access to the state (e.g. Heun's
    # intermediate stage). ForwardEuler.initialize is the identity so this is
    # a no-op for the common case.
    ts = Terrarium.initialize(terrarium_timestepper(land), terrarium_model(land), state)
    set_terrarium_timestepper!(land, ts)

    Tsoil = interior(state.temperature)[:, 1, end] .+ NF(273.15)
    sat = interior(state.saturation_water_ice)[:, 1, end]
    vars.prognostic.land.soil_temperature .= Tsoil
    vars.prognostic.land.soil_moisture .= sat
    return nothing
end

function timestep!(
        vars::Variables,
        land::AbstractTerrariumLandModel,
        ::PrimitiveEquation,
    )
    state = vars.prognostic.land.terrarium
    tmodel = terrarium_model(land)
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
    set!(inputs.air_temperature, Tair)
    set!(inputs.air_temperature, inputs.air_temperature - NF(273.15))   # K -> °C
    set!(inputs.air_pressure, pres)
    set!(inputs.air_pressure, exp(inputs.air_pressure))                  # log(Pa) -> Pa
    set!(inputs.specific_humidity, humid)
    set!(inputs.rainfall, rain)
    set!(inputs.snowfall, snow)
    set!(inputs.windspeed, wind)
    set!(inputs.surface_shortwave_down, Rsd)
    set!(inputs.surface_longwave_down, Rld)

    # Wrap the SpeedyWeather-owned state in an integrator and step it forward
    # for the duration of the SpeedyWeather step using `terrarium_substep`.
    integrator = ModelIntegrator(
        state.clock, tmodel, InputSources(),
        state, terrarium_initializers(land), terrarium_timestepper(land),
    )
    Terrarium.run!(integrator; period = vars.prognostic.clock.Δt, Δt = terrarium_substep(land))

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
        vars.prognostic.land.surface_humidity_flux .= state.latent_heat_flux ./ consts.Llg
    end
    if haskey(vars.parameterizations, :surface_longwave_up)
        vars.parameterizations.surface_longwave_up .= state.surface_longwave_up
    end
    if haskey(vars.parameterizations, :surface_shortwave_up)
        vars.parameterizations.surface_shortwave_up .= state.surface_shortwave_up
    end
    return nothing
end

end # module
