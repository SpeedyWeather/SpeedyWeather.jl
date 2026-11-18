# function barrier for all sea ice
@propagate_inbounds sea_ice_timestep!(vars::Variables, model::PrimitiveEquation) =
    timestep!(vars, model.sea_ice, model)

function variables(::AbstractSeaIce)
    return (
        # prescribed sea ice has no step dimension as its not time stepped, dynamic sea ice defines it with one
        PrognosticVariable(:sea_ice_concentration, GridXY(), namespace = :ocean, desc = "Sea ice concentration", units = "1"),
    )
end

function variables(::AbstractDynamicSeaIce, model::AbstractModel)
    nsteps = get_nsteps(model.time_stepping, :ocean)
    pg = nsteps.prognostic
    tg = nsteps.tendency
    return (
        PrognosticVariable(:sea_ice_concentration, OceanXYT(pg), namespace = :ocean, desc = "Sea ice concentration", units = "1"),
        TendencyVariable(:sea_ice_concentration, OceanXYT(tg), namespace = :ocean, desc = "Tendency of sea ice concentration", units = "1/s"),
    )
end

export PrescribedSeaIce

"""Prescribed sea ice that declares the necessary allocations for the sea ice concentration,
but all dynamics are expected to be externally set. Used for coupling to external sea ice models."""
struct PrescribedSeaIce <: AbstractPrescribedSeaIce end
PrescribedSeaIce(::SpectralGrid) = PrescribedSeaIce() # added constructor, just to be consistent with call signatures
initialize!(vars::Variables, ::PrescribedSeaIce, model) = nothing
timestep!(vars::Variables, ::PrescribedSeaIce, model) = nothing

export ThermodynamicSeaIce

"""Thermodynamic sea ice model with the sea ice volume per unit ocean area `V` [m] as prognostic variable.
The sea ice concentration ℵ and the ice thickness `h` follow from `V`: the thickness is tied to the concentration,
`h(ℵ) = h_thin + (h_pack - h_thin)ℵ`, so that low concentrations represent thin, new ice and full cover thick pack
ice, with `V = ℵ h(ℵ)`, see `sea_ice_concentration`. Beyond full cover (`V > h_pack`) the ice thickens. Heat is
conducted through the ice (zero-layer Semtner, 1976), `F_c = k (T_f - T_i) / h`, with an ice surface temperature
`T_i` diagnosed from the ice surface energy balance in the longwave radiation, see
`ocean_surface_emission_temperature`. Freezing, bottom melt and top melt all change `V` by the energy they use,
so that the ocean heat content minus the latent heat stored in sea ice is conserved. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct ThermodynamicSeaIce{NF} <: AbstractDynamicSeaIce
    "[OPTION] Freezing temperature of sea water [K]"
    @param freezing_temperature::NF = 273.15 - 1.8 (bounds = Positive,)

    "[OPTION] Melting temperature of the sea ice surface, upper bound of the ice surface temperature [K]"
    @param melting_temperature::NF = 273.15 (bounds = Positive,)

    "[OPTION] Thermal conductivity of sea ice [W/m/K]"
    @param conductivity::NF = 2 (bounds = Positive,)

    "[OPTION] Ice thickness at vanishing concentration, i.e. of new, thin ice [m]"
    @param thin_ice_thickness::NF = 0.2 (bounds = Positive,)

    "[OPTION] Ice thickness at full cover, i.e. of pack ice [m], at least `thin_ice_thickness`"
    @param pack_ice_thickness::NF = 2.5 (bounds = Positive,)

    "[OPTION] Density of sea ice [kg/m³]"
    @param density::NF = 917 (bounds = Positive,)

    "[OPTION] Latent heat of fusion of sea ice [J/kg]"
    @param latent_heat_fusion::NF = 3.34e5 (bounds = Positive,)

    "[OPTION] Basal melt coefficient, ice-ocean heat flux per Kelvin above freezing [W/m²/K]"
    @param basal_melt_coefficient::NF = 6.1 (bounds = Nonnegative,)

    "[OPTION] Heat capacity per area [J/m²/K] with which ice freezes over a prescribed (not interactive) ocean, 50 m of water"
    @param prescribed_ocean_heat_capacity::NF = 4184 * 1000 * 50 (bounds = Positive,)
end

# immutable and adapted: a core component of the model, so the radiation can dispatch on it inside (GPU) kernels
Adapt.@adapt_structure ThermodynamicSeaIce
ThermodynamicSeaIce(SG::SpectralGrid; kwargs...) = ThermodynamicSeaIce{SG.NF}(; kwargs...)
initialize!(::ThermodynamicSeaIce, ::AbstractModel) = nothing

function variables(::ThermodynamicSeaIce, model::AbstractModel)
    nsteps = get_nsteps(model.time_stepping, :ocean)
    pg = nsteps.prognostic
    tg = nsteps.tendency
    return (
        # the ice volume is the prognostic variable, time stepped with the top melt as tendency
        PrognosticVariable(:sea_ice_volume, OceanXYT(pg), namespace = :ocean, desc = "Sea ice volume per unit ocean area", units = "m"),
        TendencyVariable(:sea_ice_volume, OceanXYT(tg), namespace = :ocean, desc = "Tendency of sea ice volume", units = "m/s"),
        # the concentration follows from the volume (no tendency, so not time stepped), read by albedo, fluxes, radiation
        PrognosticVariable(:sea_ice_concentration, OceanXYT(pg), namespace = :ocean, desc = "Sea ice concentration", units = "1"),
        # diagnosed in the longwave radiation from the ice surface energy balance
        ParameterizationVariable(:sea_ice_surface_temperature, Grid2D(), namespace = :ocean, desc = "Sea ice surface temperature", units = "K"),
        ParameterizationVariable(:sea_ice_conductive_flux, Grid2D(), namespace = :ocean, desc = "Heat conducted from the ocean through the sea ice, per ocean area", units = "W/m^2"),
        ParameterizationVariable(:sea_ice_top_melt, Grid2D(), namespace = :ocean, desc = "Heat flux melting sea ice at its surface, per ocean area", units = "W/m^2"),
        # written by the albedo, for the ice surface energy balance
        ParameterizationVariable(:sea_ice_albedo, Grid2D(), namespace = :ocean, desc = "Albedo of the sea ice part of the ocean", units = "1"),
    )
end

"""$(TYPEDSIGNATURES)
Sea ice thickness [m] at concentration `ℵ`, linear between `thin_ice_thickness` at ℵ = 0 and
`pack_ice_thickness` at ℵ = 1."""
@inline sea_ice_thickness(ℵ, sea_ice::ThermodynamicSeaIce) =
    sea_ice_thickness(ℵ, sea_ice.thin_ice_thickness, sea_ice.pack_ice_thickness)
@inline sea_ice_thickness(ℵ, h_thin, h_pack) = h_thin + (max(h_pack, h_thin) - h_thin) * ℵ   # pack at least thin

"""$(TYPEDSIGNATURES)
Sea ice thickness [m] of the ice volume per unit ocean area `V` [m]: `h(ℵ(V))` up to full cover, beyond which
(`V > pack_ice_thickness`) the ice thickens, i.e. `h = V` at full cover."""
@inline sea_ice_thickness_from_volume(V, sea_ice::ThermodynamicSeaIce) =
    sea_ice_thickness_from_volume(V, sea_ice.thin_ice_thickness, sea_ice.pack_ice_thickness)
@inline sea_ice_thickness_from_volume(V, h_thin, h_pack) =
    sea_ice_thickness(sea_ice_concentration(V, h_thin, h_pack), h_thin, h_pack) + max(V - max(h_pack, h_thin), zero(V))

"""$(TYPEDSIGNATURES)
Sea ice volume per unit ocean area [m] at concentration `ℵ ∈ [0, 1]`, `V = ℵ h(ℵ)`."""
@inline sea_ice_volume(ℵ, sea_ice::ThermodynamicSeaIce) = sea_ice_volume(ℵ, sea_ice.thin_ice_thickness, sea_ice.pack_ice_thickness)
@inline sea_ice_volume(ℵ, h_thin, h_pack) = ℵ * sea_ice_thickness(ℵ, h_thin, h_pack)

"""$(TYPEDSIGNATURES)
Sea ice concentration [1] of the ice volume per unit ocean area `V` [m], the exact inverse of `V = ℵ h(ℵ)`,
written such that it never divides by `pack_ice_thickness - thin_ice_thickness` (which may be zero) and its
denominator is at least `2thin_ice_thickness`. One for `V ≥ pack_ice_thickness`, zero for `V ≤ 0`. Branch-free."""
@inline sea_ice_concentration(V, sea_ice::ThermodynamicSeaIce) =
    sea_ice_concentration(V, sea_ice.thin_ice_thickness, sea_ice.pack_ice_thickness)
@inline function sea_ice_concentration(V, h_thin, h_pack)
    h_pack = max(h_pack, h_thin)    # pack ice at least as thick as thin ice, else the sqrt could be of a negative
    V = clamp(V, zero(V), h_pack)
    ℵ = 2V / (h_thin + sqrt(h_thin^2 + 4 * (h_pack - h_thin) * V))
    return min(ℵ, one(ℵ))
end

"""$(TYPEDSIGNATURES)
Ice surface temperature `T_i` [K] from the zero-layer (Semtner, 1976) surface energy balance of sea ice

    ϵσT_i⁴ = absorbed_radiation + conductance (T_f - T_i)

with `absorbed_radiation` [W/m²] the absorbed shortwave and downward longwave radiation, `conductance = k/h` [W/m²/K]
that of the ice
and `T_f` the freezing temperature at the ice base, solved with a fixed number (branch-free) of Newton steps
from `T_start`, then capped at the melting temperature `T_m`. Returns `T_i`, the conductive flux `F_c = conductance (T_f - T_i)`
[W/m²] (positive upward, cooling the ocean) and the surplus `max(absorbed_radiation + F_c - ϵσT_i⁴, 0)` [W/m²] that
melts ice at the top when `T_i` is capped. Four Newton steps converge also for the thinnest ice (largest conductance)."""
@inline function sea_ice_surface_balance(absorbed_radiation, conductance, ϵσ, T_freeze, T_melt, T_start)
    # f(T) = ϵσT⁴ + G(T - T_f) - A is convex and increasing for T > 0, Newton converges from any T > 0
    T = T_start
    for _ in 1:4
        f = ϵσ * T^4 + conductance * (T - T_freeze) - absorbed_radiation
        dfdT = 4ϵσ * T^3 + conductance      # > 0, as the conductance > 0
        T -= f / dfdT
    end
    T = min(T, T_melt)                      # the ice surface cannot be warmer than melting
    F_c = conductance * (T_freeze - T)                # conducted from the ice base (ocean) to the ice surface
    # energy left at the surface to melt ice at T = T_melt, not negative (the Newton residual stays with the ocean)
    surplus = max(absorbed_radiation + F_c - ϵσ * T^4, zero(T))
    return T, F_c, surplus
end

# initialize the ice volume from the concentration (zero unless set before) after the ocean initial conditions
function initialize!(vars::Variables, sea_ice::ThermodynamicSeaIce, model::PrimitiveEquation)
    V = vars.prognostic.ocean.sea_ice_volume.data
    ℵ = vars.prognostic.ocean.sea_ice_concentration.data
    NF = eltype(V)
    (; thin_ice_thickness, pack_ice_thickness) = sea_ice
    ocean = model.land_sea_mask.land_fraction.data .< 1     # at least partially ocean
    V .= ifelse.(ocean, sea_ice_volume.(clamp.(ℵ, zero(NF), one(NF)), thin_ice_thickness, pack_ice_thickness), zero(NF))
    ℵ .= sea_ice_concentration.(V, thin_ice_thickness, pack_ice_thickness)
    return nothing
end

function timestep!(vars::Variables, sea_ice::ThermodynamicSeaIce, model::PrimitiveEquation)
    # the top melt [W/m²] was computed by the longwave radiation (zero without radiation)
    dV = get_tendency_step(vars.tendencies.ocean.sea_ice_volume, time_stepper(model.time_stepping, :ocean), sea_ice)
    top_melt = vars.parameterizations.ocean.sea_ice_top_melt
    (; land_fraction) = model.land_sea_mask

    @boundscheck size(dV) == size(top_melt) || throw(BoundsError())
    ρL⁻¹ = inv(sea_ice.density * sea_ice.latent_heat_fusion)       # [m³/J] ice volume per energy

    launch!(architecture(dV), LinearWorkOrder, size(dV), sea_ice_kernel!, dV, top_melt, land_fraction, ρL⁻¹)
    return nothing
end

@kernel inbounds = true function sea_ice_kernel!(dV, top_melt, land_fraction, ρL⁻¹)
    ij = @index(Global, Linear)     # every grid point ij
    if land_fraction[ij] < 1        # at least partially ocean
        dV[ij] = -top_melt[ij] * ρL⁻¹               # melt at the top [m/s]
    end
end

# applied after the time step: the ocean freezes ice or melts it from below as exact energy adjustments
function Base.filter!(vars::Variables, sea_ice::ThermodynamicSeaIce, model::PrimitiveEquation)
    V = vars.prognostic.ocean.sea_ice_volume.data
    ℵ = vars.prognostic.ocean.sea_ice_concentration.data
    NF = eltype(V)
    ocean = model.land_sea_mask.land_fraction.data .< 1     # at least partially ocean

    if haskey(vars.prognostic.ocean, :sea_surface_temperature)
        sst = vars.prognostic.ocean.sea_surface_temperature.data
        Δt = time_step(model.time_stepping, :ocean, vars.prognostic.clock)
        ρL = sea_ice.density * sea_ice.latent_heat_fusion   # [J/m³] latent heat stored per ice volume
        T_f = sea_ice.freezing_temperature
        γ = sea_ice.basal_melt_coefficient

        # only an interactive ocean (with an SST tendency) pays and receives the latent heat
        interactive = haskey(vars.tendencies.ocean, :sea_surface_temperature)
        heat_capacity = convert(NF, interactive ? ocean_heat_capacity(model.ocean, sea_ice) : sea_ice.prescribed_ocean_heat_capacity)
        sst_change = NF(interactive)                        # 1 to change the SST, 0 to keep it

        # 1. top melt beyond the ice present (V < 0): that heat warms the ocean instead
        sst .+= ifelse.(ocean, sst_change .* ρL .* max.(.-V, zero(NF)) ./ heat_capacity, zero(NF))
        V .= max.(V, zero(NF))

        # 2. freezing: the ocean's heat deficit below freezing becomes ice
        heat_deficit = ifelse.(ocean, heat_capacity .* max.(T_f .- sst, zero(NF)), zero(NF))     # [J/m²]
        if interactive
            # the ocean pays for the ice and returns to freezing, no upper bound (the ice thickens)
            V .+= heat_deficit ./ ρL
            sst .+= heat_deficit ./ heat_capacity
        else
            # a prescribed SST stays below freezing, so the same deficit would freeze again every step:
            # grow ice towards full cover (as main's concentration clamp) but never thicken it beyond
            V .= max.(V, min.(V .+ heat_deficit ./ ρL, sea_ice.pack_ice_thickness))
        end

        # 3. bottom melt: a basal flux above freezing melts ice, limited to the ice present, the ocean pays for it
        ΔV = ifelse.(ocean, min.(γ .* max.(sst .- T_f, zero(NF)) .* Δt ./ ρL, V), zero(NF))
        V .-= ΔV
        sst .-= sst_change .* ρL .* ΔV ./ heat_capacity
    end

    (; thin_ice_thickness, pack_ice_thickness) = sea_ice
    V .= max.(V, zero(NF))
    ℵ .= sea_ice_concentration.(V, thin_ice_thickness, pack_ice_thickness)
    return nothing
end

"""$(TYPEDSIGNATURES)
Heat capacity per area [J/m²/K] of the ocean with which `sea_ice` freezes and melts. That of the mixed layer for a
`SlabOcean`, otherwise `prescribed_ocean_heat_capacity` of the sea ice model."""
ocean_heat_capacity(ocean, sea_ice::ThermodynamicSeaIce) =
    hasproperty(ocean, :heat_capacity_mixed_layer) ? ocean.heat_capacity_mixed_layer : sea_ice.prescribed_ocean_heat_capacity

# SEA ICE IN THE LONGWAVE RADIATION, see `ocean_surface_emission_temperature` in longwave_radiation.jl
needs_surface_longwave_down(::ThermodynamicSeaIce) = true

"""$(TYPEDSIGNATURES)
`ThermodynamicSeaIce`: the ocean part of grid point `ij` is split into sea ice (concentration ℵ) and open water.
The ice surface temperature `T_i` is solved from the ice surface energy balance (see `sea_ice_surface_balance`)
with the absorbed shortwave (ice albedo), the downward longwave and the heat conducted through the ice of thickness
`h(V)`. Stores `T_i`, the conductive flux and the top melt (both per ocean area) and returns the flux-equivalent
temperature `T* = ((1 - ℵ) sst⁴ + ℵ T_i⁴)^(1/4)`, so that `emissivity * σ * T*⁴` is exactly the emission of the
open water and the ice together. Without ice (ℵ = 0) this is exactly `sst`. Note: unlike the other methods, this one
also writes these sea ice diagnostics, which the slab ocean and the sea ice time step use."""
@propagate_inbounds function ocean_surface_emission_temperature(
        ij, vars, sst, surface_longwave_down, emissivity, sea_ice::ThermodynamicSeaIce, scheme, model
    )
    progn_ocean = vars.prognostic.ocean
    params_ocean = vars.parameterizations.ocean
    ts = time_stepper(model.time_stepping, :ocean)
    sea_surface_temperature = sst[ij]
    NF = typeof(sea_surface_temperature)

    ℵ = get_prognostic_step(progn_ocean.sea_ice_concentration, ts, scheme)[ij]
    V = get_prognostic_step(progn_ocean.sea_ice_volume, ts, scheme)[ij]
    surface_shortwave_down = haskey(vars.parameterizations, :surface_shortwave_down) ?
        vars.parameterizations.surface_shortwave_down[ij] : zero(NF)
    surface_air_temperature = haskey(vars.parameterizations, :surface_air_temperature) ?
        vars.parameterizations.surface_air_temperature[ij] : sea_surface_temperature

    # ice surface energy balance: absorbed radiation, conductance k/h of the ice, emission ϵσT⁴
    absorbed_radiation = (1 - params_ocean.sea_ice_albedo[ij]) * surface_shortwave_down + surface_longwave_down
    conductance = sea_ice.conductivity / sea_ice_thickness_from_volume(V, sea_ice)
    ϵσ = emissivity * model.atmosphere.stefan_boltzmann
    T_start = min(sea_surface_temperature, surface_air_temperature)
    T_i, F_c, surplus = sea_ice_surface_balance(
        absorbed_radiation, conductance, ϵσ, sea_ice.freezing_temperature, sea_ice.melting_temperature, T_start
    )

    params_ocean.sea_ice_surface_temperature[ij] = T_i
    # [W/m²] per ocean area, zero (not 0·NaN) without ice, e.g. over land with a NaN sea surface temperature
    params_ocean.sea_ice_conductive_flux[ij] = ifelse(ℵ > 0, ℵ * F_c, zero(NF))     # cools the ocean
    params_ocean.sea_ice_top_melt[ij] = ifelse(ℵ > 0, ℵ * surplus, zero(NF))        # melts ice at the top

    # flux-equivalent emission temperature of open water and ice, exactly the SST without ice
    T⁴ = (1 - ℵ) * sea_surface_temperature^4 + ℵ * T_i^4
    return ifelse(ℵ > 0, sqrt(sqrt(T⁴)), sea_surface_temperature)
end
