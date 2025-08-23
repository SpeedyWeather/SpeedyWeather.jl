abstract type AbstractSeaIce <: AbstractModelComponent end

# function barrier for all oceans
function sea_ice_timestep!( progn::PrognosticVariables,
                            diagn::DiagnosticVariables,
                            model::PrimitiveEquation)
    timestep!(progn, diagn, model.sea_ice, model)
end

export ThermodynamicSeaIce
@kwdef mutable struct ThermodynamicSeaIce{NF} <: AbstractSeaIce
    "[OPTION] Freezing temperature of sea water [K]"
    temp_freeze::NF = 273.15-1.8

    "[OPTION] Melting rate of sea ice [m²/m²/s/K]"
    melt_rate::NF = 5e-5

    "[OPTION] Freezing rate of sea ice [m²/m²/K]"
    freeze_rate::NF = 0.12
end

ThermodynamicSeaIce(SG::SpectralGrid; kwargs...) = ThermodynamicSeaIce{SG.NF}(;kwargs...)
initialize!(::ThermodynamicSeaIce, ::AbstractModel) = nothing

# don't affect concentration (may be set with set!)
function initialize!(
    ocean::PrognosticVariablesOcean,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    sea_ice_model::ThermodynamicSeaIce,
    model::PrimitiveEquation,
)
    return nothing
end

function timestep!( progn::PrognosticVariables,
                    diagn::DiagnosticVariables,
                    sea_ice_model::ThermodynamicSeaIce,
                    model::PrimitiveEquation)
    
    sst = progn.ocean.sea_surface_temperature
    ℵ = progn.ocean.sea_ice_concentration   # sea ice concentration [0, 1] as \aleph yay!
    
    Δt = model.time_stepping.Δt_sec
    (; mask) = model.land_sea_mask

    m = sea_ice_model.melt_rate             # melt rate [m²/m²/s/K]
    f_Δt = sea_ice_model.freeze_rate / Δt   # include 1/Δt here as SST below freezing is proportional to Δt
    temp_freeze = sea_ice_model.temp_freeze

    launch!(architecture(ℵ), LinearWorkOrder, size(ℵ), sea_ice_kernel!, ℵ, sst, mask, temp_freeze, m, f_Δt, Δt)

    return nothing
end

@kernel inbounds=true function sea_ice_kernel!(ℵ, sst, mask, @Const(temp_freeze), @Const(m), @Const(f_Δt), @Const(Δt))
    ij = @index(Global, Linear)    # every grid point ij

    if mask[ij] < 1 && isfinite(sst[ij])        # at least partially ocean, SST not NaN (=masked)

        # ice-sst flux as a relaxation term wrt to freezing, with different melt/freeze rates
        dT = sst[ij] - temp_freeze              # uncorrected difference to freezing temperature
        F = -m*max(dT, 0) - f_Δt*min(dT, 0)     # melt if above freezing, freeze if below
        sst[ij] = max(sst[ij], temp_freeze)     # cap sst at freezing

        # Euler forward step to update sea ice concentration, cap between [0, 1]
        ℵ[ij] += Δt*F
        ℵ[ij] = max(min(ℵ[ij], 1), 0)
    end
end