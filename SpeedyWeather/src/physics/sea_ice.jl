abstract type AbstractSeaIce <: AbstractModelComponent end

# function barrier for all sea ice
@propagate_inbounds function sea_ice_timestep!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        model::PrimitiveEquation
    )
    return timestep!(progn, diagn, model.sea_ice, model)
end

export ThermodynamicSeaIce

"""Thermodynamic sea ice model using sea ice concentration as
prognostic variable, also modifies sea surface temperature as
cooling below freezing grows sea ice. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef mutable struct ThermodynamicSeaIce{NF} <: AbstractSeaIce
    "[OPTION] Freezing temperature of sea water [K]"
    @param temp_freeze::NF = 273.15 - 1.8 (bounds=Positive,)

    "[OPTION] Melting rate of sea ice [m²/m²/s/K]"
    @param melt_rate::NF = 1.0e-6 (bounds=Nonnegative,)

    "[OPTION] Freezing rate of sea ice [m²/m²/K]"
    @param freeze_rate::NF = 0.1 (bounds=Nonnegative,)
end

ThermodynamicSeaIce(SG::SpectralGrid; kwargs...) = ThermodynamicSeaIce{SG.NF}(; kwargs...)
initialize!(::ThermodynamicSeaIce, ::AbstractModel) = nothing

function variables(::ThermodynamicSeaIce)
    return (
        PrognosticVariable(name = :sea_ice_concentration, dims = Grid2D(), namespace = :ocean, desc = "Sea ice concentration", units = "1"),
    )
end

# don't affect concentration (may be set with set!)
function initialize!(
        ocean::PrognosticVariablesOcean,
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        sea_ice_model::ThermodynamicSeaIce,
        model::PrimitiveEquation,
    ) where {PrognosticVariablesOcean}
    return nothing
end

function timestep!(
        progn::PrognosticVariables,
        diagn::DiagnosticVariables,
        sea_ice_model::ThermodynamicSeaIce,
        model::PrimitiveEquation
    )

    sst = progn.ocean.sea_surface_temperature
    ℵ = progn.ocean.sea_ice_concentration   # sea ice concentration [0, 1] as \aleph yay!

    Δt = model.time_stepping.Δt_sec
    (; mask) = model.land_sea_mask

    m = sea_ice_model.melt_rate             # melt rate [m²/m²/s/K]
    f = sea_ice_model.freeze_rate           # freeze rate [m²/m²/K]
    temp_freeze = sea_ice_model.temp_freeze
    parameters = (; m, f, temp_freeze, Δt)

    launch!(
        architecture(ℵ), LinearWorkOrder, size(ℵ), sea_ice_kernel!,
        ℵ, sst, mask, parameters
    )
    return nothing
end

@kernel inbounds = true function sea_ice_kernel!(ℵ, sst, mask, parameters)
    ij = @index(Global, Linear)    # every grid point ij

    if mask[ij] < 1 && isfinite(sst[ij])        # at least partially ocean, SST not NaN (=masked)

        (; m, f, temp_freeze, Δt) = parameters

        # ice-sst flux as a relaxation term wrt to freezing, with different melt/freeze rates
        # include 1/Δt here as SST below freezing is proportional to Δt
        dT = sst[ij] - temp_freeze              # uncorrected difference to freezing temperature
        F = -m * max(dT, 0) - f / Δt * min(dT, 0)     # melt if above freezing, freeze if below
        sst[ij] = max(sst[ij], temp_freeze)     # cap sst at freezing

        # Euler forward step to update sea ice concentration, cap between [0, 1]
        ℵ[ij] += Δt * F
        ℵ[ij] = max(min(ℵ[ij], 1), 0)
    end
end
