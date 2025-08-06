abstract type AbstractSeaIce <: AbstractModelComponent end

export NoSeaIce
struct NoSeaIce <: AbstractSeaIce end
NoSeaIce(::SpectralGrid) = NoSeaIce()
initialize!(::NoSeaIce, ::AbstractModel) = nothing

# set all concentration to zero
function initialize!(
    ocean::PrognosticVariablesOcean,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    sea_ice_model::NoSeaIce,
    model::PrimitiveEquation,
)
    ocean.sea_ice_concentration .= 0
end

# function barrier for all oceans
function sea_ice_timestep!( progn::PrognosticVariables,
                            diagn::DiagnosticVariables,
                            model::PrimitiveEquation)
    sea_ice_timestep!(progn, diagn, model.sea_ice, model)
end

# NoSeaIce does not do anything 
function sea_ice_timestep!( progn::PrognosticVariables,
                            diagn::DiagnosticVariables,
                            sea_ice_model::NoSeaIce,
                            model::PrimitiveEquation)
    return nothing
end

export ThermodynamicSeaIce
@kwdef mutable struct ThermodynamicSeaIce{NF} <: AbstractSeaIce
    "[OPTION] Freezing temperature of sea water [K]"
    temp_freeze::NF = 273.15-1.8

    "[OPTION] Melting rate of sea ice [m²/m²/s/K]"
    melt_rate::NF = 5e-5

    "[OPTION] Freezing rate of sea ice [m²/m²/s/K]"
    freeze_rate::NF = 5e-5
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

function sea_ice_timestep!( progn::PrognosticVariables,
                            diagn::DiagnosticVariables,
                            sea_ice_model::ThermodynamicSeaIce,
                            model::PrimitiveEquation)
    
    sst = progn.ocean.sea_surface_temperature
    ice = progn.ocean.sea_ice_concentration

    Δt = model.time_stepping.Δt_sec
    (; mask) = model.land_sea_mask

    m = sea_ice_model.melt_rate
    f = sea_ice_model.freeze_rate
    temp_freeze = sea_ice_model.temp_freeze

    # Euler forward step
    @inbounds for ij in eachgridpoint(sst, ice, mask)
        if mask[ij] < 1     # at least partially ocean

            # ice-sst flux as a relaxation term wrt to freezing, with different melt/freeze rates
            dT = sst[ij] - temp_freeze              # uncorrected difference to freezing temperature
            F = -m*max(dT, 0) - f*min(dT, 0)        # melt if above freezing, freeze if below
            sst[ij] = max(sst[ij], temp_freeze)     # cap sst at freezing

            # update sea ice concentration, cap between [0, 1]
            ice[ij] += Δt*F
            ice[ij] = max(min(ice[ij], 1), 0)
        end
    end

    return nothing
end