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
    capacity::NF = 1e-4
    temp_freeze::NF = 273.15-1.8
    growth::NF = 3e-1
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

# NoSeaIce does not do anything 
function sea_ice_timestep!( progn::PrognosticVariables,
                            diagn::DiagnosticVariables,
                            sea_ice_model::ThermodynamicSeaIce,
                            model::PrimitiveEquation)
    
    sst = progn.ocean.sea_surface_temperature
    ice = progn.ocean.sea_ice_concentration

    Δt = model.time_stepping.Δt_sec
    (; mask) = model.land_sea_mask

    γ = sea_ice_model.capacity
    g = sea_ice_model.growth
    temp_freeze = sea_ice_model.temp_freeze

    # Euler forward step
    @inbounds for ij in eachgridpoint(sst, ice, mask)
        if mask[ij] < 1     # at least partially ocean

            F = -γ*(sst[ij] - temp_freeze)          # ice-sst flux as a relaxation term
            F += -min(sst[ij] - temp_freeze, 0)/Δt  # move flux below freezing to ice growth
            sst[ij] = max(sst[ij], temp_freeze)     # cap sst at freezing

            # update concentration, cap between [0, 1]
            ice[ij] += Δt*F*g
            ice[ij] = max(min(ice[ij], 1), 0)
        end
    end

    return nothing
end