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
@propagate_inbounds function parameterization!(ij, diagn, progn, longwave::UniformCooling, model)

    T = diagn.grid.temp_grid_prev
    dTdt = diagn.tendencies.temp_tend_grid
    (; temp_min, temp_stratosphere) = longwave
    nlayers = size(T, 2)

    NF = eltype(T)
    cooling = -inv(convert(NF, Second(longwave.time_scale).value))
    τ⁻¹ = inv(convert(NF, Second(longwave.time_scale_stratosphere).value))

    for k in 1:nlayers
        # Paulius and Garner, 2006, eq (1) and (2)
        dTdt[ij, k] += T[ij, k] > temp_min ? cooling : (temp_stratosphere - T[ij, k]) * τ⁻¹
    end
    return
end

## JEEVANJEE TEMPERATURE FLUX

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
    emissivity_atmosphere::NF = 0

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

# function barrier
@propagate_inbounds function parameterization!(ij, diagn, progn, longwave::JeevanjeeRadiation, model)

    T = diagn.grid.temp_grid                            # to match Seeley, 2023 notation
    dTdt = diagn.tendencies.temp_tend_grid
    pₛ = diagn.grid.pres_grid_prev[ij]                  # surface pressure [Pa]
    nlayers = size(T, 2)

    (; α) = longwave
    τ⁻¹ = inv(convert(eltype(T), Second(longwave.time_scale).value))
    ϵ_ocean = longwave.emissivity_ocean
    ϵ_land = longwave.emissivity_land
    ϵ = longwave.emissivity_atmosphere
    σ = model.atmosphere.stefan_boltzmann
    cₚ = model.atmosphere.heat_capacity
    Tₜ = longwave.temp_tropopause

    land_fraction = model.land_sea_mask.mask[ij]
    sst = progn.ocean.sea_surface_temperature[ij]
    lst = progn.land.soil_temperature[ij, 1]            # TODO use skin temperature?

    # extension to Jeevanjee: Include temperature flux (Stefan-Boltzmann)
    # between surface and lowermost air temperature
    # but zero flux if land/sea not available
    Fₖ_ocean = isfinite(sst) ? ϵ_ocean * σ * sst^4 : zero(sst)      # [W/m²]
    diagn.physics.ocean.surface_longwave_up[ij] = Fₖ_ocean      # for ocean model forcing

    Fₖ_land = isfinite(lst) ? ϵ_land * σ * lst^4 : zero(lst)        # [W/m²]
    diagn.physics.land.surface_longwave_up[ij] = Fₖ_land        # for land model forcing

    Fₖ_down = ϵ * σ * T[ij, nlayers]^4
    diagn.physics.surface_longwave_down[ij] = Fₖ_down          # for surface energy budget
    dTdt[ij, nlayers] -= surface_flux_to_tendency(Fₖ_down / cₚ, pₛ, model)   # out of layer k

    # land-sea mask weighted combined flux from land and ocean
    local Fₖ::eltype(T)
    Fₖ = (1 - land_fraction) * Fₖ_ocean + land_fraction * Fₖ_land
    dTdt[ij, nlayers] += surface_flux_to_tendency(Fₖ / cₚ, pₛ, model) # convert [W/m²] / cₚ to K·Pa/s
    diagn.physics.surface_longwave_up[ij] = Fₖ                      # [W/m²] store for output/diagnostics

    # integrate from surface up
    for k in nlayers:-1:2
        # Seeley and Wordsworth, 2023 eq (1)
        Fₖ += (T[ij, k - 1] - T[ij, k]) * α * (Tₜ - T[ij, k])         # upward flux from layer k into k-1
        dTdt[ij, k] -= flux_to_tendency(Fₖ / cₚ, pₛ, k, model)    # out of layer k
        dTdt[ij, k - 1] += flux_to_tendency(Fₖ / cₚ, pₛ, k - 1, model)    # into layer k-1
    end

    # Relax the uppermost level towards prescribed "tropopause temperature"
    dTdt[ij, 1] += (Tₜ - T[ij, 1]) * τ⁻¹

    # for diagnostic, use Fₖ as the outgoing longwave radiation although it's technically into the
    # uppermost layer from below (not out of it)
    return diagn.physics.outgoing_longwave[ij] = Fₖ
end

function variables(::AbstractLongwave)
    return (
        DiagnosticVariable(name = :surface_longwave_down, dims = Grid2D(), desc = "Surface longwave radiation down", units = "W/m^2"),
        DiagnosticVariable(name = :surface_longwave_up, dims = Grid2D(), desc = "Surface longwave radiation up over ocean", units = "W/m^2", namespace = :ocean),
        DiagnosticVariable(name = :surface_longwave_up, dims = Grid2D(), desc = "Surface longwave radiation up over land", units = "W/m^2", namespace = :land),
        DiagnosticVariable(name = :surface_longwave_up, dims = Grid2D(), desc = "Surface longwave radiation up", units = "W/m^2"),
        DiagnosticVariable(name = :outgoing_longwave, dims = Grid2D(), desc = "TOA Longwave radiation up", units = "W/m^2"),
    )
end

export OneBandLongwave
"""One-band longwave radiation scheme with transmissivity and radiative transfer components.
Following Frierson, 2006, JAS, https://doi.org/10.1175/JAS3706.1. Fields are $(TYPEDFIELDS)"""
struct OneBandLongwave{T, R} <: AbstractLongwave
    transmissivity::T
    radiative_transfer::R
end

Adapt.@adapt_structure OneBandLongwave

# primitive wet model version
function OneBandLongwave(
        SG::SpectralGrid;
        transmissivity = FriersonLongwaveTransmissivity(SG),
        radiative_transfer = OneBandLongwaveRadiativeTransfer(SG),
    )
    return OneBandLongwave(transmissivity, radiative_transfer)
end

# primitive dry model version
export OneBandGreyLongwave
function OneBandGreyLongwave(
        SG::SpectralGrid;
        transmissivity = TransparentLongwaveTransmissivity(SG),
        radiative_transfer = OneBandLongwaveRadiativeTransfer(SG),
    )
    return OneBandLongwave(transmissivity, radiative_transfer)
end

function Base.show(io::IO, M::OneBandLongwave)
    println(io, "OneBandLongwave <: AbstractLongwave")
    properties = propertynames(M)
    n = length(properties)
    for (i, key) in enumerate(properties)
        val = getfield(M, key)
        s = i == n ? "└" : "├"  # choose ending └ for last property
        p = i == n ? print : println
        p(io, "$s $key: $(typeof(val))")
    end
    return
end

# initialize one after another
function initialize!(radiation::OneBandLongwave, model::PrimitiveEquation)
    initialize!(radiation.transmissivity, model)
    return initialize!(radiation.radiative_transfer, model)
end

@propagate_inbounds function parameterization!(ij, diagn, progn, radiation::OneBandLongwave, model)
    # pass on array that was used to compute transmissivity (scratch array)
    t = transmissivity!(ij, diagn, progn, radiation.transmissivity, model)
    return longwave_radiative_transfer!(ij, diagn, t, progn, radiation.radiative_transfer, model)
end

export OneBandLongwaveRadiativeTransfer

"""Radiative transfer solver for OneBandLongwave radiation scheme.
Computes longwave radiative transfer with upward and downward beams including surface fluxes.
Fields are $(TYPEDFIELDS)"""
@kwdef struct OneBandLongwaveRadiativeTransfer{NF} <: AbstractLongwaveRadiativeTransfer
    "[OPTION] Emissivity for surface flux over ocean [1]"
    emissivity_ocean::NF = 1

    "[OPTION] Emissivity for surface flux over land [1]"
    emissivity_land::NF = 1
end

Adapt.@adapt_structure OneBandLongwaveRadiativeTransfer

# generator function
OneBandLongwaveRadiativeTransfer(SG::SpectralGrid; kwargs...) = OneBandLongwaveRadiativeTransfer{SG.NF}(; kwargs...)
initialize!(::OneBandLongwaveRadiativeTransfer, ::PrimitiveEquation) = nothing

@propagate_inbounds function longwave_radiative_transfer!(
        ij,
        diagn,
        transmissivity,
        progn,
        longwave::OneBandLongwaveRadiativeTransfer,
        model,
    )
    T = diagn.grid.temp_grid_prev
    NF = eltype(T)
    dTdt = diagn.tendencies.temp_tend_grid
    pₛ = diagn.grid.pres_grid_prev[ij]                  # surface pressure [Pa]
    nlayers = size(T, 2)

    ϵ_ocean = longwave.emissivity_ocean
    ϵ_land = longwave.emissivity_land
    σ = model.atmosphere.stefan_boltzmann
    cₚ = model.atmosphere.heat_capacity

    land_fraction = model.land_sea_mask.mask[ij]
    sst = progn.ocean.sea_surface_temperature[ij]
    lst = progn.land.soil_temperature[ij, 1]            # TODO use skin temperature?

    U_ocean = isfinite(sst) ? ϵ_ocean * σ * sst^4 : zero(sst)       # [W/m²]
    diagn.physics.ocean.surface_longwave_up[ij] = U_ocean           # for ocean model forcing

    U_land = isfinite(lst) ? ϵ_land * σ * lst^4 : zero(lst)         # [W/m²]
    diagn.physics.land.surface_longwave_up[ij] = U_land             # for land model forcing

    # land-sea mask weighted combined flux from land and ocean (surface boundary condition)
    U::NF = (1 - land_fraction) * U_ocean + land_fraction * U_land
    dTdt[ij, nlayers] += surface_flux_to_tendency(U / cₚ, pₛ, model)     # convert [W/m²] / cₚ to K·Pa/s
    diagn.physics.surface_longwave_up[ij] = U                           # [W/m²] store for output/diagnostics

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
    diagn.physics.outgoing_longwave[ij] = U

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
    diagn.physics.surface_longwave_down[ij] = D
    dTdt[ij, nlayers] -= surface_flux_to_tendency(D / cₚ, pₛ, model)    # out of layer k

    return nothing
end
