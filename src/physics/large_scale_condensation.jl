abstract type AbstractCondensation <: AbstractParameterization end

# no condensation
large_scale_condensation!(::ColumnVariables, ::Nothing, ::PrimitiveEquation) = nothing

export ImplicitCondensation
"""
Large scale condensation as with implicit precipitation.
$(TYPEDFIELDS)"""
@kwdef mutable struct ImplicitCondensation{NF<:AbstractFloat} <: AbstractCondensation
    "[OPTION] Relative humidity threshold [1 = 100%] to trigger condensation"
    relative_humidity_threshold::NF = 0.95

    "[OPTION] Reevaporation efficiency [1/(kg/kg)], 0 for no reevaporation"
    reevaporation::NF = 30

    "[OPTION] Convert precipitation below freezing to snow?"
    snow::Bool = true

    "[OPTION] Freezing temperature for snow fall [K]"
    freezing_threshold::NF = 263

    "[OPTION] Melting temperature for snow fall [K]"
    melting_threshold::NF = 278

    "[OPTION] Time scale in multiples of time step Δt, the larger the less immediate"
    time_scale::NF = 3
end

ImplicitCondensation(SG::SpectralGrid; kwargs...) = ImplicitCondensation{SG.NF}(; kwargs...)

# nothing to initialize with this scheme
initialize!(scheme::ImplicitCondensation, model::PrimitiveEquation) = nothing

# do nothing fall back for primitive dry 
function large_scale_condensation!( 
    column::ColumnVariables,
    model::PrimitiveEquation,
)
    return nothing
end

# function barrier for all AbstractCondensation
function large_scale_condensation!( 
    column::ColumnVariables,
    model::PrimitiveWet,
)
    # dispatch by large scale condensation type
    large_scale_condensation!(column, model.large_scale_condensation, model)
end

# function barrier for ImplicitCondensation to unpack model
function large_scale_condensation!( 
    column::ColumnVariables,
    condensation::ImplicitCondensation,
    model::PrimitiveWet,
)
    saturation_humidity!(column, model.clausius_clapeyron)
    large_scale_condensation!(column, condensation,
        model.clausius_clapeyron, model.geometry, model.planet, model.atmosphere, model.time_stepping)
end

"""
$(TYPEDSIGNATURES)
Large-scale condensation for a `column` by relaxation back to 100%
relative humidity. Calculates the tendencies for specific humidity
and temperature from latent heat release and integrates the
large-scale precipitation vertically for output."""
function large_scale_condensation!(
    column::ColumnVariables,
    condensation::ImplicitCondensation,
    clausius_clapeyron::AbstractClausiusClapeyron,
    geometry::Geometry,
    planet::AbstractPlanet,
    atmosphere::AbstractAtmosphere,
    time_stepping::AbstractTimeStepper,
)

    (; pres, temp, humid) = column          # prognostic vars (from previous time step for numerical stability)
    (; temp_tend, humid_tend) = column      # tendencies to write into
    (; sat_humid) = column                  # intermediate variable, calculated in thermodynamics!

    # precompute scaling constants to minimize divisions (used to convert between humidity [kg/kg] and precipitation [m])
    pₛ = pres[end]                          # surface pressure
    (; Δt_sec) = time_stepping
    Δσ = geometry.σ_levels_thick            # layer thickness in sigma coordinates
    pₛ_gρ = pₛ/(planet.gravity * atmosphere.water_density)

    # thermodynamics
    Lᵥ = clausius_clapeyron.latent_heat_condensation    # latent heat of vaporization
    Lᵢ = clausius_clapeyron.latent_heat_fusion          # latent heat of fusion
    cₚ = clausius_clapeyron.heat_capacity               # heat capacity
    Rᵥ = clausius_clapeyron.R_vapour                    # gas constant for water vapour
    Lᵢ_cₚ = Lᵢ / cₚ

    (; time_scale, relative_humidity_threshold, freezing_threshold, melting_threshold) = condensation
    let_it_snow = condensation.snow                     # flag to switch snow on/off
    rain_flux_down = zero(eltype(column))               # start with zero rain/snow flux at top of atmosphere
    snow_flux_down = zero(eltype(column))               # both have units of [m] (precipitation)

    @inbounds for k in eachindex(column)

        # Condensation from humidity in this layer (for a negative humidity tendency)
        δq_cond = sat_humid[k] * relative_humidity_threshold - humid[k]
        
        # skip if no condensation has occurred yet in this layer or above
        if δq_cond < 0 || snow_flux_down > 0 || rain_flux_down > 0

            # 0. convert between humidity tendency [kg/kg/s] and precipitation amount [m]
            Δp_gρ = Δσ[k] * pₛ_gρ                           # pressure thickness of layer Δp times 1/g/ρ [m]
            ΔpΔt_gρ = Δp_gρ*Δt_sec                          # pressure thickness of layer Δp times Δt/g/ρ [ms]

            # 1. Melting of snow from layer above
            δT_melt = max(temp[k] - melting_threshold, 0)   # only if temperature above melting threshold
            E_avail = cₚ * δT_melt                          # available Energy [J / kg (of air)] for melting
            melt_depth = (E_avail / Lᵢ) * ΔpΔt_gρ           # "depth" of snow (in rainwater amount) [m] that can be melted
            melt_amount = min(snow_flux_down, melt_depth)   # cap to snow flux amount [m], don't melt more than available
            snow_flux_down -= melt_amount                   # move melted snow to rain
            rain_flux_down += melt_amount                   # this can evaporate now too
            δq_melt = melt_amount / Δp_gρ                   # convert back to humidity increase over timestep [kg/kg]

            # 2. Reevaporation and condensation
            r = condensation.reevaporation * max(δq_cond, 0)    # reevaporation efficiency [1], scale linearly with dryness
            rain_evaporated = min(r, 1) * rain_flux_down        # [m], min(r, 1) to not evaporate more than available
            rain_flux_down -= rain_evaporated                   # remove reevaporated rain
            δq_evap = rain_evaporated / Δp_gρ                   # convert to humidity tendency over timestep [kg/kg]
            δq = min(0, δq_cond) + δq_evap                      # sum with condensation and melting

            # effective latent heat L as weighted average from condensation/evaporation and melting (fusion)
            L = (Lᵥ * abs(δq) + Lᵢ * δq_melt) / (abs(δq) + δq_melt)
            δq += δq_melt                                       # only add now to allow distinction in line above

            # Solve for melting of snow, condensation, reevaporation (and possibly sublimation) implicitly in time
            # implicit correction, Frierson et al. 2006 eq. (21)
            # derivative of qsat wrt to temp, 1/cₚ included for fewer divisions
            dqsat_dT = sat_humid[k] * relative_humidity_threshold * L/(cₚ*Rᵥ*temp[k]^2)
            δq /= ((1 + L*dqsat_dT) * time_scale*Δt_sec)
            δT = -L/cₚ * δq                     # latent heat release for enthalpy conservation

            # If there is large-scale condensation at a level higher (i.e. smaller k) than
            # the cloud-top previously diagnosed due to convection, then increase the cloud-top
            # Fortran SPEEDY documentation Page 7 (last sentence)
            column.cloud_top = min(column.cloud_top, k)

            # 2. Precipitation (rain) due to large-scale condensation [kg/m²/s] /ρ for [m/s]
            δq = min(0, δq)                     # precipitation only for negative humidity tendency
            rain = ΔpΔt_gρ * -δq                # precipitation rain [m] on this layer k, Formula 4
            snow = zero(rain)                   # start with zero snow but potentially swap below

            # decide whether to turn precip into snow (all rain freezes to snow)
            rain, snow = let_it_snow && temp[k] < freezing_threshold ? (snow, rain) : (rain, snow)
            rain_flux_down += rain              # accumulate into downward fluxes [m] (used in layer below)
            snow_flux_down += snow              # accumulate into downward fluxes [m] (used in layer below)

            # latent heat release when freezing for enthalpy conservation
            δT = -(snow > 0) * Lᵢ_cₚ * δq        # snow for negative δq (condensation then freezing)

            # only accumulate into humid_tend now to allow humid_tend != 0 before this scheme is called
            humid_tend[k] += δq
            temp_tend[k] += δT
        end   
    end

    # precipitation from rain/snow whatever is fluxed out 
    column.precip_large_scale = rain_flux_down      # vertical integral [m]
	column.snow_large_scale   = snow_flux_down      # vertical integral [m]

    # convert to rain/snow fall rate [m/s]
    column.precip_rate_large_scale = rain_flux_down / Δt_sec
    column.snow_rate_large_scale = snow_flux_down / Δt_sec

    return nothing
end
