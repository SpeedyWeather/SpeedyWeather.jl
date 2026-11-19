abstract type AbstractLongwave <: AbstractRadiation end
abstract type AbstractLongwaveRadiativeTransfer <: AbstractLongwave end

export UniformCooling

"""
Uniform cooling following Paulius and Garner, 2006. JAS. https://doi.org/10.1175/JAS3705.1
imposing a default temperature tendency of -1.5K/day (=1K/16hours for a `time_scale` of 16 hours)
on every level except for the stratosphere (diagnosed as `temp < temp_min`) where a
relaxation term with `time_scale_stratosphere` towards `temp_stratosphere` is applied.

    dT/dt = -1.5K/day for T > 207.5K else (200K-T) / 5 days

Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct UniformCooling{NF} <: AbstractLongwave
    "[OPTION] time scale of cooling, default = -1.5K/day = -1K/16hrs"
    time_scale::Second = Hour(16)

    "[OPTION] temperature [K] below which stratospheric relaxation is applied"
    @param temp_min::NF = 207.5 (bounds = Positive,)

    "[OPTION] target temperature [K] of stratospheric relaxation"
    @param temp_stratosphere::NF = 200 (bounds = Positive,)

    "[OPTION] time scale of stratospheric relaxation"
    time_scale_stratosphere::Second = Day(5)
end

Adapt.@adapt_structure UniformCooling
UniformCooling(SG::SpectralGrid; kwargs...) = UniformCooling{SG.NF}(; kwargs...)
initialize!(radiation::UniformCooling, model::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds function parameterization!(ij, vars, longwave::UniformCooling, model)

    T = get_prognostic_step(vars.grid.temperature, model.time_stepping, longwave)
    dTdt = get_tendency_step(vars.tendencies.grid.temperature, model.time_stepping, longwave)
    (; temp_min, temp_stratosphere) = longwave
    nlayers = size(T, 2)

    NF = eltype(T)
    cooling = -inv(convert(NF, Second(longwave.time_scale).value))
    τ⁻¹ = inv(convert(NF, Second(longwave.time_scale_stratosphere).value))

    for k in 1:nlayers
        # Paulius and Garner, 2006, eq (1) and (2)
        dTdt[ij, k] += ifelse(T[ij, k] > temp_min, cooling, (temp_stratosphere - T[ij, k]) * τ⁻¹)
    end
    return nothing
end

# JEEVANJEE TEMPERATURE FLUX
export JeevanjeeRadiation
"""
Temperature flux longwave radiation from Jeevanjee and Romps, 2018,
following Seeley and Wordsworth, 2023, eq (1)

    dF/dT = α*(T_t - T)

with `F` the upward temperature flux between two layers with temperature
difference `dT`, `α = 0.025 W/m²/K²` and `T_t = 200K` a prescribed tropopause
temperature. Flux into the lowermost layer is 0. Flux out of uppermost
layer also 0, but `dT/dt = (T_t - T)/τ` is added to relax the uppermost
layer towards the tropopause temperature `T_t` with time scale `τ = 24h`
(Seeley and Wordsworth, 2023 use 6h, which is unstable a low resolutions here).
Fields are
$(TYPEDFIELDS)"""
@parameterized @kwdef struct JeevanjeeRadiation{NF} <: AbstractLongwave
    "[OPTION] Radiative forcing constant (W/m²/K²)"
    @param α::NF = 0.025 (bounds = Nonnegative,)

    "[OPTION] Emissivity of the atmosphere [1]"
    @param emissivity_atmosphere::NF = 0 (bounds = 0 .. 1,)

    "[OPTION] Emissivity for surface flux over ocean [1]"
    @param emissivity_ocean::NF = 0.65 (bounds = 0 .. 1,)

    "[OPTION] Emissivity for surface flux over land [1]"
    @param emissivity_land::NF = 0.6 (bounds = 0 .. 1,)

    "[OPTION] Tropopause temperature [K]"
    @param temp_tropopause::NF = 200 (bounds = Positive,)

    "[OPTION] Tropopause relaxation time scale to temp_tropopause"
    time_scale::Second = Hour(24)
end

Adapt.@adapt_structure JeevanjeeRadiation
JeevanjeeRadiation(SG::SpectralGrid; kwargs...) = JeevanjeeRadiation{SG.NF}(; kwargs...)
initialize!(::JeevanjeeRadiation, ::PrimitiveEquation) = nothing

# the sea ice model, or nothing for models without one
get_sea_ice(model) = hasproperty(model, :sea_ice) ? model.sea_ice : nothing

"""$(TYPEDSIGNATURES)
Temperature [K] at which the ocean part of grid point `ij` emits longwave radiation, dispatched on the sea ice
model `model.sea_ice`. Over open water that is the sea surface temperature `sst`. `scheme` is the radiation
component that reads the prognostic step. Used by all radiation schemes that emit from the ocean surface, so
that a sea ice model with an ice surface temperature of its own only needs a new method here."""
@propagate_inbounds ocean_surface_emission_temperature(ij, vars, sst, scheme, model) =
    ocean_surface_emission_temperature(ij, vars, sst, get_sea_ice(model), scheme, model)

"""$(TYPEDSIGNATURES)
Sea ice (concentration ℵ) without an ice surface temperature of its own (`PrescribedSeaIce`, `ThermodynamicSeaIce`,
or `sea_ice = nothing` with a manually set concentration) is assumed to have the surface air temperature but not to
be warmer than the water underneath, i.e. `(1 - ℵ) * sst + ℵ * min(sst, surface_air_temperature)`. Without a sea ice
concentration this is the sea surface temperature. The surface air temperature is computed in the `:boundary_layer`
(`SurfaceCondition`), which runs before `:radiation` (see `model.parameterizations`), so it is from the same time step.
A custom `parameterizations` tuple with sea ice has to keep `:boundary_layer` before `:radiation`."""
@propagate_inbounds function ocean_surface_emission_temperature(ij, vars, sst, ::Union{Nothing, AbstractSeaIce}, scheme, model)
    progn_ocean = vars.prognostic.ocean
    ts = time_stepper(model.time_stepping, :ocean)
    sea_surface_temperature = sst[ij]
    sea_ice_concentration = haskey(progn_ocean, :sea_ice_concentration) ?
        get_prognostic_step(progn_ocean.sea_ice_concentration, ts, scheme)[ij] : zero(sea_surface_temperature)
    surface_air_temperature = haskey(vars.parameterizations, :surface_air_temperature) ?
        vars.parameterizations.surface_air_temperature[ij] : sea_surface_temperature
    # sea ice emits at the surface air temperature, but never warmer than the water underneath
    return (1 - sea_ice_concentration) * sea_surface_temperature +
        sea_ice_concentration * min(sea_surface_temperature, surface_air_temperature)
end

# function barrier
@propagate_inbounds function parameterization!(ij, vars, longwave::JeevanjeeRadiation, model)

    T = get_prognostic_step(vars.grid.temperature, model.time_stepping, longwave)
    dTdt = get_tendency_step(vars.tendencies.grid.temperature, model.time_stepping, longwave)
    pₛ = vars.parameterizations.surface_pressure[ij]            # surface pressure [Pa]
    nlayers = size(T, 2)

    (; α) = longwave
    τ⁻¹ = inv(convert(eltype(T), Second(longwave.time_scale).value))
    ϵ_ocean = longwave.emissivity_ocean
    ϵ_land = longwave.emissivity_land
    ϵ = longwave.emissivity_atmosphere
    σ = model.atmosphere.stefan_boltzmann
    cₚ = model.atmosphere.heat_capacity
    Tₜ = longwave.temp_tropopause

    land_fraction = model.land_sea_mask.land_fraction[ij]
    sst = get_prognostic_step(vars.prognostic.ocean.sea_surface_temperature, time_stepper(model.time_stepping, :ocean), longwave)
    lst = get_prognostic_step(vars.prognostic.land.soil_temperature, time_stepper(model.time_stepping, :land), longwave) # TODO use skin temperature?

    # extension to Jeevanjee: Include temperature flux (Stefan-Boltzmann)
    # between surface and lowermost air temperature
    # but zero flux if land/sea not available
    Fₖ_ocean = ϵ_ocean * σ * ocean_surface_emission_temperature(ij, vars, sst, longwave, model)^4  # [W/m²], sea ice aware
    vars.parameterizations.ocean.surface_longwave_up[ij] = Fₖ_ocean     # for ocean model forcing

    Fₖ_land = ϵ_land * σ * lst[ij, 1]^4                                 # [W/m²]
    vars.parameterizations.land.surface_longwave_up[ij] = Fₖ_land       # for land model forcing

    Fₖ_down = ϵ * σ * T[ij, nlayers]^4
    vars.parameterizations.surface_longwave_down[ij] = Fₖ_down          # for surface energy budget
    dTdt[ij, nlayers] -= surface_flux_to_tendency(Fₖ_down / cₚ, pₛ, model)   # out of layer k

    # land-sea mask weighted combined flux from land and ocean
    Fₖ::eltype(T) = (1 - land_fraction) * Fₖ_ocean + land_fraction * Fₖ_land
    dTdt[ij, nlayers] += surface_flux_to_tendency(Fₖ / cₚ, pₛ, model)   # convert [W/m²] / cₚ to K·Pa/s
    vars.parameterizations.surface_longwave_up[ij] = Fₖ                 # [W/m²] store for output/diagnostics

    # integrate from surface up
    for k in nlayers:-1:2
        # Seeley and Wordsworth, 2023 eq (1)
        Fₖ += (T[ij, k - 1] - T[ij, k]) * α * (Tₜ - T[ij, k])           # upward flux from layer k into k-1
        dTdt[ij, k] -= flux_to_tendency(Fₖ / cₚ, pₛ, k, model)          # out of layer k
        dTdt[ij, k - 1] += flux_to_tendency(Fₖ / cₚ, pₛ, k - 1, model)  # into layer k-1
    end

    # Relax the uppermost level towards prescribed "tropopause temperature"
    dTdt[ij, 1] += (Tₜ - T[ij, 1]) * τ⁻¹

    # for diagnostic, use Fₖ as the outgoing longwave radiation although it's technically into the
    # uppermost layer from below (not out of it)
    return vars.parameterizations.outgoing_longwave[ij] = Fₖ
end

function variables(::AbstractLongwave)
    return (
        ParameterizationVariable(:surface_longwave_down, Grid2D(), desc = "Surface longwave radiation down", units = "W/m^2"),
        ParameterizationVariable(:surface_longwave_up, Grid2D(), desc = "Surface longwave radiation up over ocean", units = "W/m^2", namespace = :ocean),
        ParameterizationVariable(:surface_longwave_up, Grid2D(), desc = "Surface longwave radiation up over land", units = "W/m^2", namespace = :land),
        ParameterizationVariable(:surface_longwave_up, Grid2D(), desc = "Surface longwave radiation up", units = "W/m^2"),
        ParameterizationVariable(:outgoing_longwave, Grid2D(), desc = "Top-of-atmosphere longwave radiation up", units = "W/m^2"),
    )
end

export OneBandLongwave
"""One-band longwave radiation scheme with transmissivity and radiative transfer components.
Following Frierson, 2006, JAS, https://doi.org/10.1175/JAS3753.1. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct OneBandLongwave{T, R} <: AbstractLongwave
    @component transmissivity::T
    @component radiative_transfer::R
end

Adapt.@adapt_structure OneBandLongwave

# primitive wet model version
function OneBandLongwave(
        SG::SpectralGrid;
        transmissivity = ByrneOGormanLongwaveTransmissivity(SG),
        radiative_transfer = OneBandLongwaveRadiativeTransfer(SG, stratospheric_emissivity = 0.05),
    )
    return OneBandLongwave(transmissivity, radiative_transfer)
end

# primitive dry model version (not a type only a function)
export OneBandGreyLongwave
function OneBandGreyLongwave(
        SG::SpectralGrid;
        transmissivity = ConstantLongwaveTransmissivity(SG),
        radiative_transfer = OneBandLongwaveRadiativeTransfer(SG),
    )
    return OneBandLongwave(transmissivity, radiative_transfer)
end

Base.show(io::IO, M::OneBandLongwave) = Base.show(io, M, values = false)

# initialize one after another
function initialize!(radiation::OneBandLongwave, model::PrimitiveEquation)
    initialize!(radiation.transmissivity, model)
    initialize!(radiation.radiative_transfer, model)
    return nothing
end

@propagate_inbounds function parameterization!(ij, vars, radiation::OneBandLongwave, model)
    # pass on array that was used to compute transmissivity (scratch array)
    t = transmissivity!(ij, vars, radiation.transmissivity, model)
    longwave_radiative_transfer!(ij, vars, t, radiation.radiative_transfer, model)
    return nothing
end

export OneBandLongwaveRadiativeTransfer

"""Radiative transfer solver for OneBandLongwave radiation scheme.
Computes longwave radiative transfer with upward and downward beams including surface fluxes.
Optionally with an additional stratospheric emission directly to space following the
stratospheric correction in Fortran SPEEDY (`epslw`): the stratosphere [0, `stratosphere_pressure`]
emits `stratospheric_emissivity * σT⁴`, distributed on the layers by their pressure overlap,
representing emission in the centre of the CO₂ band that a grey scheme with a transparent
stratosphere lacks, and needed to balance stratospheric ozone heating.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct OneBandLongwaveRadiativeTransfer{NF} <: AbstractLongwaveRadiativeTransfer
    "[OPTION] Emissivity for surface flux over ocean [1]"
    @param emissivity_ocean::NF = 0.98 (bounds = 0 .. 1,)

    "[OPTION] Emissivity for surface flux over land [1]"
    @param emissivity_land::NF = 0.98 (bounds = 0 .. 1,)

    "[OPTION] Additional emissivity of the stratosphere, emitting directly to space [1]"
    @param stratospheric_emissivity::NF = 0 (bounds = 0 .. 1,)

    "[OPTION] Lower boundary of the stratosphere for the stratospheric emission [Pa]"
    stratosphere_pressure::NF = 14000
end

Adapt.@adapt_structure OneBandLongwaveRadiativeTransfer

# generator function
OneBandLongwaveRadiativeTransfer(SG::SpectralGrid; kwargs...) = OneBandLongwaveRadiativeTransfer{SG.NF}(; kwargs...)
initialize!(::OneBandLongwaveRadiativeTransfer, ::PrimitiveEquation) = nothing

@propagate_inbounds function longwave_radiative_transfer!(
        ij,
        vars,
        transmissivity,
        longwave::OneBandLongwaveRadiativeTransfer,
        model,
    )
    T = get_prognostic_step(vars.grid.temperature, model.time_stepping, longwave)
    NF = eltype(T)
    dTdt = get_tendency_step(vars.tendencies.grid.temperature, model.time_stepping, longwave)
    pₛ = vars.parameterizations.surface_pressure[ij]                    # surface pressure [Pa]
    nlayers = size(T, 2)

    ϵ_ocean = longwave.emissivity_ocean
    ϵ_land = longwave.emissivity_land
    σ = model.atmosphere.stefan_boltzmann
    cₚ = model.atmosphere.heat_capacity

    land_fraction = model.land_sea_mask.land_fraction[ij]
    sst = get_prognostic_step(vars.prognostic.ocean.sea_surface_temperature, time_stepper(model.time_stepping, :ocean), longwave)
    # TODO use skin temperature?
    lst = get_prognostic_step(vars.prognostic.land.soil_temperature, time_stepper(model.time_stepping, :land), longwave)            

    U_ocean = ϵ_ocean * σ * ocean_surface_emission_temperature(ij, vars, sst, longwave, model)^4  # [W/m²], sea ice aware
    vars.parameterizations.ocean.surface_longwave_up[ij] = U_ocean      # for ocean model forcing

    U_land = ϵ_land * σ * lst[ij, 1]^4                                  # [W/m²], 1 for top layer
    vars.parameterizations.land.surface_longwave_up[ij] = U_land        # for land model forcing

    # land-sea mask weighted combined flux from land and ocean (surface boundary condition)
    U::NF = (1 - land_fraction) * U_ocean + land_fraction * U_land
    dTdt[ij, nlayers] += surface_flux_to_tendency(U / cₚ, pₛ, model)     # convert [W/m²] / cₚ to K·Pa/s
    vars.parameterizations.surface_longwave_up[ij] = U                  # [W/m²] store for output/diagnostics

    # UPWARD BEAM
    for k in nlayers:-1:2
        t = transmissivity[ij, k]
        U = U * t + (1 - t) * σ * T[ij, k]^4
        dTdt[ij, k] -= flux_to_tendency(U / cₚ, pₛ, k, model)           # out of layer k
        dTdt[ij, k - 1] += flux_to_tendency(U / cₚ, pₛ, k - 1, model)   # into layer k-1
    end

    # Outgoing longwave radiation at TOA
    t = transmissivity[ij, 1]
    U = U * t + (1 - t) * σ * T[ij, 1]^4
    dTdt[ij, 1] -= flux_to_tendency(U / cₚ, pₛ, 1, model)               # out of layer 1

    # Additional stratospheric emission directly to space, distributed by pressure overlap
    ϵ_strat = longwave.stratospheric_emissivity
    p_strat = longwave.stratosphere_pressure
    coord = model.geometry.vertical_coordinates
    for k in 1:nlayers                                                  # overlap = 0 below the stratosphere
        p_top = pressure_above(k, pₛ, coord)
        overlap = max(0, min(pressure_below(k, pₛ, coord), p_strat) - p_top) / p_strat
        E = ϵ_strat * overlap * σ * T[ij, k]^4
        dTdt[ij, k] -= flux_to_tendency(E / cₚ, pₛ, k, model)
        U += E
    end
    vars.parameterizations.outgoing_longwave[ij] = U

    # DOWNWARD BEAM
    D::NF = 0               # top boundary condition (no longwave coming from space)
    for k in 1:(nlayers - 1)
        t = transmissivity[ij, k]
        D = D * t + (1 - t) * σ * T[ij, k]^4
        dTdt[ij, k] -= flux_to_tendency(D / cₚ, pₛ, k, model)           # out of layer k
        dTdt[ij, k + 1] += flux_to_tendency(D / cₚ, pₛ, k + 1, model)   # into layer k+1
    end

    # Surface downward longwave radiation
    t = transmissivity[ij, nlayers]
    D = D * t + (1 - t) * σ * T[ij, nlayers]^4
    vars.parameterizations.surface_longwave_down[ij] = D
    dTdt[ij, nlayers] -= surface_flux_to_tendency(D / cₚ, pₛ, model)    # out of layer k

    return nothing
end
