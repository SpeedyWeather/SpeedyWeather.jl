# function barrier for all sea ice
@propagate_inbounds sea_ice_timestep!(vars::Variables, model::PrimitiveEquation) =
    timestep!(vars, model.sea_ice, model)

function variables(::AbstractSeaIce)
    return (
        # all prognostic variables need a step dimension, by default use only one
        PrognosticVariable(:sea_ice_concentration, Grid3D(1), namespace = :ocean, desc = "Sea ice concentration", units = "1"),
    )
end

export PrescribedSeaIce

"""Prescribed sea ice that declares the necessary allocations for the sea ice concentration,
but all dynamics are expected to be externally set. Used for coupling to external sea ice models."""
struct PrescribedSeaIce <: AbstractSeaIce end
PrescribedSeaIce(::SpectralGrid) = PrescribedSeaIce() # added constructor, just to be consistent with call signatures
initialize!(vars::Variables, ::PrescribedSeaIce, model) = nothing
timestep!(vars::Variables, ::PrescribedSeaIce, model) = nothing

export ThermodynamicSeaIce

"""Thermodynamic sea ice model using sea ice concentration as
prognostic variable, also modifies sea surface temperature as
cooling below freezing grows sea ice. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef mutable struct ThermodynamicSeaIce{NF} <: AbstractSeaIce
    "[OPTION] Freezing temperature of sea water [K]"
    @param freezing_temperature::NF = 273.15 - 1.8 (bounds = Positive,)

    "[OPTION] Melting rate of sea ice [m²/m²/s/K]"
    @param melt_rate::NF = 1.0e-7 (bounds = Nonnegative,)

    "[OPTION] Freezing rate of sea ice [m²/m²/K]"
    @param freeze_rate::NF = 0.5 (bounds = Nonnegative,)
end

ThermodynamicSeaIce(SG::SpectralGrid; kwargs...) = ThermodynamicSeaIce{SG.NF}(; kwargs...)
initialize!(::ThermodynamicSeaIce, ::AbstractModel) = nothing
initialize!(vars::Variables, sea_ice_model::ThermodynamicSeaIce, model::PrimitiveEquation) = nothing

function variables(::ThermodynamicSeaIce, model::AbstractModel)
    nsteps = get_nsteps(model.time_stepping, model)
    pg = nsteps.prognostic_grid
    tg = nsteps.tendency_grid
    return (
        PrognosticVariable(:sea_ice_concentration, GridXYT(pg), namespace = :ocean, desc = "Sea ice concentration", units = "1"),
        TendencyVariable(:sea_ice_concentration, GridXYT(tg), namespace = :ocean, desc = "Sea ice concentration", units = "1"),
    )
end

# leapfrog when possible
@inline which_prognostic_step(var, ::AbstractLeapfrog, ::ThermodynamicSeaIce) = 2

function timestep!(vars::Variables, sea_ice_model::ThermodynamicSeaIce, model::PrimitiveEquation)
    # escape immediately if there is no sea surface temperature
    haskey(vars.prognostic.ocean, :sea_surface_temperature) || return nothing

    # if ocean does not have an SST tendency use scratch array to write into the void
    dsst = haskey(vars.tendencies.ocean, :sea_surface_temperature) ?
        get_tendency_step(vars.tendencies.ocean.sea_surface_temperature, model.time_stepping, sea_ice_model) :
        vars.scratch.grid.a_2D

    # sea ice concentration written as \aleph yay!
    dℵ = get_tendency_step(vars.tendencies.ocean.sea_ice_concentration, model.time_stepping, sea_ice_model)  
    sst = get_prognostic_step(vars.prognostic.ocean.sea_surface_temperature, model.time_stepping, model.ocean)
    
    (; Δt) = model.time_stepping
    (; land_fraction) = model.land_sea_mask

    m = sea_ice_model.melt_rate             # melt rate [m²/m²/s/K]
    f = sea_ice_model.freeze_rate           # freeze rate [m²/m²/K]
    T_freeze = sea_ice_model.freezing_temperature
    parameters = (; m, f, T_freeze, Δt)

    launch!(
        architecture(dℵ), LinearWorkOrder, size(dℵ), sea_ice_kernel!,
        dℵ, dsst, sst, land_fraction, parameters
    )
    return nothing
end

@kernel inbounds = true function sea_ice_kernel!(dℵ, dsst, sst, land_fraction, parameters)
    ij = @index(Global, Linear)     # every grid point ij

    if land_fraction[ij] < 1        # at least partially ocean
        (; m, f, T_freeze, Δt) = parameters

        # ice-sst flux as a relaxation term wrt to freezing, with different melt/freeze rates
        # include 1/Δt here as SST below freezing is proportional to Δt
        dT = sst[ij] - T_freeze                     # uncorrected difference to freezing temperature
        dsst[ij] -= min(dT, 0)/Δt                   # increase SST to freezing temperature
        F = -m * max(dT, 0) - f / Δt * min(dT, 0)   # melt if above freezing, freeze if below
        dℵ[ij] = F                                  # sea ice tendency
    end
end

# applied after the time stepping for any kind of "hacky" correction
function Base.filter!(vars::Variables, sea_ice_model::ThermodynamicSeaIce, model::PrimitiveEquation)
    # clamp sea ice concentration in [0, 1]
    ℵ = vars.prognostic.ocean.sea_ice_concentration
    ℵ .= max.(min.(ℵ, 1), 0)
    return nothing
end