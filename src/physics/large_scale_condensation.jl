abstract type AbstractCondensation <: AbstractParameterization end

export ImplicitCondensation
"""Large-scale condensation with implicit latent heat release.
$(TYPEDFIELDS)"""
@kwdef struct ImplicitCondensation{NF} <: AbstractCondensation
    "[OPTION] Relative humidity threshold [1 = 100%] to trigger condensation"
    relative_humidity_threshold::NF = 0.95

    "[OPTION] Reevaporation efficiency [1/(kg/kg)], 0 for no reevaporation"
    reevaporation::NF = 30

    "[OPTION] Convert rain below freezing to snow?"
    snow::Bool = true

    "[OPTION] Freezing temperature for snow fall [K]"
    freezing_threshold::NF = 263

    "[OPTION] Melting temperature for snow fall [K]"
    melting_threshold::NF = 278

    "[OPTION] Time scale in multiples of time step Δt, the larger the less immediate"
    time_scale::NF = 3
end

Adapt.@adapt_structure ImplicitCondensation
ImplicitCondensation(SG::SpectralGrid; kwargs...) = ImplicitCondensation{SG.NF}(; kwargs...)
initialize!(::ImplicitCondensation, ::PrimitiveEquation) = nothing

# function barrier
@propagate_inbounds function parameterization!(ij, diagn, progn, lsc::ImplicitCondensation, model)
    (; clausius_clapeyron, geometry, planet, atmosphere, time_stepping) = model
    large_scale_condensation!(ij, diagn, lsc, clausius_clapeyron, geometry, planet, atmosphere, time_stepping)
end

"""
$(TYPEDSIGNATURES)
Large-scale condensation for a `column` by relaxation back to 100%
relative humidity. Calculates the tendencies for specific humidity
and temperature from latent heat release and integrates the
large-scale precipitation vertically for output."""
@propagate_inbounds function large_scale_condensation!(
    ij,
    diagn::DiagnosticVariables,
    condensation::ImplicitCondensation,
    clausius_clapeyron::AbstractClausiusClapeyron,
    geometry::Geometry,
    planet::AbstractPlanet,
    atmosphere::AbstractAtmosphere,
    time_stepping,
)
    # use previous time step for more stable Euler forward step of the parameterizations
    temp = diagn.grid.temp_grid_prev                    # temperature [K]
    humid = diagn.grid.humid_grid_prev                  # specific humidity [kg/kg]
    temp_tend = diagn.tendencies.temp_tend_grid         # temperature tendency [K/s]
    humid_tend = diagn.tendencies.humid_tend_grid       # specific humidity tendency [kg/kg/s]
    nlayers = size(temp, 2)

    # precompute scaling constants to minimize divisions (used to convert between humidity [kg/kg] and precipitation [m])
    pₛ = diagn.grid.pres_grid_prev[ij]                  # surface pressure [Pa]
    (; Δt_sec) = time_stepping                          # time step [s]
    σ = geometry.σ_levels_full                          # vertical sigma coordinate [1]
    Δσ = geometry.σ_levels_thick                        # layer thickness in sigma coordinates
    ρ = atmosphere.water_density                        # air density [kg/m³]
    pₛ_gρ = pₛ/(planet.gravity * ρ)                     # [Pa / (m/s² * kg/m³)] = [Pa m² * s² / kg] = [m]

    # thermodynamics
    Lᵥ = clausius_clapeyron.latent_heat_condensation    # latent heat of vaporization
    Lᵢ = clausius_clapeyron.latent_heat_fusion          # latent heat of fusion
    cₚ = clausius_clapeyron.heat_capacity               # heat capacity
    Rᵥ = clausius_clapeyron.R_vapour                    # gas constant for water vapour
    Lᵥ_cₚ = Lᵥ / cₚ
    Lᵢ_cₚ = Lᵢ / cₚ
    cₚ_Lᵢ = cₚ / Lᵢ

    (; time_scale, relative_humidity_threshold, freezing_threshold, melting_threshold) = condensation
    let_it_snow = condensation.snow                     # flag to switch snow on/off
    rain_flux_down = zero(eltype(temp))                 # start with zero rain/snow flux at top of atmosphere
    snow_flux_down = zero(eltype(temp))                 # both have units of [m/s] (precipitation)

    for k in 1:nlayers

        # Condensation from humidity in this layer (for a negative humidity tendency)
        # relative to threshold that can be <100%, e.g. 95%
        sat_humid_k = saturation_humidity(temp[ij, k], σ[k]*pₛ, clausius_clapeyron)
        δq_cond = sat_humid_k * relative_humidity_threshold - humid[ij, k]
        
        # skip if no condensation has occurred yet in this layer or above
        if δq_cond < 0 || snow_flux_down > 0 || rain_flux_down > 0

            # 0. convert between humidity tendency [kg/kg/s] and precipitation amount [m] or rate [m/s]
            Δp_gρ = Δσ[k] * pₛ_gρ                           # pressure thickness of layer Δp times 1/g/ρ [m]
            Δp_Δtgρ = Δp_gρ/Δt_sec                          # pressure thickness of layer Δp times 1/Δt/g/ρ [m/s]
            Δtgρ_Δp = inv(Δp_Δtgρ)                          # [s/m]

            # 1. Melting of snow from layer above
            δT_melt = max(temp[ij, k] - melting_threshold, 0)   # only if temperature above melting threshold
            E_avail = cₚ_Lᵢ * δT_melt                       # available energy [1] in multiples of latent heat [J/kg of air] for melting
            melt_rate_max = E_avail * Δp_Δtgρ               # max rate of snow melt [m/s]
            melt_rate = min(snow_flux_down, melt_rate_max)  # cap to snow flux [m/s], don't melt more than available
            snow_flux_down -= melt_rate                     # move melted snow to rain
            rain_flux_down += melt_rate                     # this can evaporate now too
            δq_melt = melt_rate * Δtgρ_Δp                   # convert back to humidity increase over timestep [kg/kg]
            δT = -Lᵢ_cₚ * δq_melt                           # just to calculate the latent heat required (explicit in time)
                                                            # don't use δq_melt as humidity tendency as it's only snow->rain

            # 2. Reevaporation and condensation
            dq = sat_humid_k - humid[ij, k]                 # reevaporation relative to 100% humidity
            r = condensation.reevaporation * max(dq, 0)     # reevaporation efficiency [1], scale linearly with dryness
            rain_evaporated = min(r, 1) * rain_flux_down    # [m/s], min(r, 1) to not evaporate more than available
            rain_flux_down -= rain_evaporated               # remove reevaporated rain
            δq_evap = rain_evaporated * Δtgρ_Δp             # convert to humidity tendency over timestep [kg/kg]
            δq = min(0, δq_cond) + δq_evap                  # [kg/kg] sum with condensation (negative) and evaporation (positive)
                                                            # division by timestep for tendency below

            # Solve for melting of snow, condensation, reevaporation (and possibly sublimation) implicitly in time
            # implicit correction, Frierson et al. 2006 eq. (21)
            # derivative of qsat wrt to temp
            T = temp[ij, k] + Δt_sec*δT                     # use temperature minus latent heat for melting in gradient
            dqsat_dT = sat_humid_k * relative_humidity_threshold * Lᵥ_cₚ/(Rᵥ*T^2)
            δq /= ((1 + Lᵥ_cₚ*dqsat_dT) * time_scale*Δt_sec)
            δT = -Lᵥ_cₚ * δq                                # latent heat release for enthalpy conservation

            # If there is large-scale condensation at a level higher (i.e. smaller k) than
            # the cloud-top previously diagnosed due to convection, then increase the cloud-top
            # Fortran SPEEDY documentation Page 7 (last sentence)
            diagn.physics.cloud_top[ij] = min(diagn.physics.cloud_top[ij], k)

            # 2. Precipitation (rain) due to large-scale condensation [kg/m²/s] /ρ for [m/s]
            δq_rain = min(0, δq)                # precipitation only for negative humidity tendency
            rain = Δp_gρ * -δq_rain             # precipitation rain [m/s] on this layer k, Formula 4
            snow = zero(rain)                   # start with zero snow but potentially swap below

            # decide whether to turn precip into snow (all rain freezes to snow)
            rain, snow = let_it_snow && temp[ij, k] < freezing_threshold ? (snow, rain) : (rain, snow)
            rain_flux_down += rain              # accumulate into downward fluxes [m/s] (used in layer below)
            snow_flux_down += snow              # accumulate into downward fluxes [m/s] (used in layer below)

            # latent heat release when freezing for enthalpy conservation
            δT -= (snow > 0) * Lᵢ_cₚ * δq_rain  # snow for negative δq (condensation then freezing)

            # only accumulate into humid_tend now to allow humid_tend != 0 before this scheme is called
            humid_tend[ij, k] += δq             # [kg/kg/s]
            temp_tend[ij, k] += δT              # [K/s]
        end
    end

    # avoid negative precipitation from rounding errors 
    rain_flux_down = max(rain_flux_down, 0)
    snow_flux_down = max(snow_flux_down, 0)

    # precipitation from rain/snow [m] during time step whatever is fluxed out 
    diagn.physics.rain_large_scale[ij] += Δt_sec*rain_flux_down
	diagn.physics.snow_large_scale[ij] += Δt_sec*snow_flux_down

    # TODO use [kg/m²/s] units? then * ρ here
    # and the rain/snow fall rate [m/s]
    diagn.physics.rain_rate_large_scale[ij] = rain_flux_down
    diagn.physics.snow_rate_large_scale[ij] = snow_flux_down

    return nothing
end

function variables(::ImplicitCondensation)
    return (
        DiagnosticVariable(name=:rain_large_scale, dims=Grid2D(), desc="Large-scale precipitation (rain, accumulated)", units="m"),
        DiagnosticVariable(name=:snow_large_scale, dims=Grid2D(), desc="Large-scale precipitation (snow, accumulated)", units="m"),
        DiagnosticVariable(name=:rain_rate_large_scale, dims=Grid2D(), desc="Large-scale precipitation (rain)", units="m/s"),
        DiagnosticVariable(name=:snow_rate_large_scale, dims=Grid2D(), desc="Large-scale precipitation (snow)", units="m/s"),
        DiagnosticVariable(name=:cloud_top, dims=Grid2D(), desc="Cloud top layer index", units="1"),
    )
end
