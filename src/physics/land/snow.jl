abstract type AbstractSnow <: AbstractParameterization end

export SnowModel    # maybe change for a more concise name later

"""
    SnowModel(; melting_threshold=275, runoff_time_scale=Year(1))

Single-column snow bucket model in equivalent liquid water depth. Snow accumulates
from the diagnosed precipitation, melts once the top soil layer exceeds
`melting_threshold`, and relaxes back to zero on the `runoff_time_scale`.
$(TYPEDFIELDS)"""
@kwdef mutable struct SnowModel{NF} <: AbstractSnow
    melting_threshold::NF = 275
    runoff_time_scale::Second = Year(1) 
end

# generator function
SnowModel(SG::SpectralGrid; kwargs...) = SnowModel{SG.NF}(; kwargs...)

# initialize component
initialize!(snow::SnowModel, model::PrimitiveEquation) = nothing

# set initial conditions for snow depth in initial conditions
function initialize!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    snow::SnowModel,
    model::PrimitiveEquation,
)
    set!(progn, model.geometry, snow_depth=0)
end

function timestep!(
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    snow::SnowModel,
    model::PrimitiveEquation,
)
    Δt = model.time_stepping.Δt_sec
    (; snow_depth) = progn.land                             # in equivalent liquid water height [m]
    (; soil_temperature) = progn.land
    
    (; mask) = model.land_sea_mask
    
    # Some thermodynamics needed by snow
	ρ_soil = model.land.thermodynamics.soil_density			# soil density [kg/m³]
	ρ_water = model.atmosphere.water_density				# water density [kg/m³]
    Lᵢ = model.clausius_clapeyron.latent_heat_fusion      	# latent heat of fusion
    cₛ = model.land.thermodynamics.heat_capacity_dry_soil
    z₁ = model.land.geometry.layer_thickness[1]
    (; melting_threshold) = snow
    r⁻¹ = inv(Second(snow.runoff_time_scale).value)

    # Snowfall rate in [kg/m²/s]
    snow_fall_rate = diagn.physics.snow_rate
    snow_melt_rate = diagn.physics.land.snow_melt_rate
    snow_runoff_rate = diagn.physics.land.snow_runoff_rate

    params = (;melting_threshold,cₛ,ρ_soil,z₁,Δt,ρ_water,Lᵢ,r⁻¹)

    launch!(architecture(snow_depth), LinearWorkOrder, size(snow_depth), land_snow_kernel!,
        snow_depth, soil_temperature, snow_melt_rate, snow_runoff_rate, snow_fall_rate, mask,
        params,
    )
	synchronize(architecture(snow_depth))
end

@kernel inbounds=true function land_snow_kernel!(
    snow_depth, soil_temperature, snow_melt_rate, snow_runoff_rate, snow_fall_rate, mask,
    params,
    )
    ij = @index(Global, Linear)             # every grid point ij

    if mask[ij] > 0                         # at least partially land
		
		(;melting_threshold, cₛ, ρ_soil, z₁, Δt, ρ_water, Lᵢ, r⁻¹) = params
        # check for melting of snow if temperature above melting threshold
        δT_melt = max(soil_temperature[ij, 1] - melting_threshold, 0)
	
        # energy available from soil warming above melting threshold [J/m²/s]
        E_avail = cₛ * δT_melt * z₁ / Δt  # [J/(m³ K)] * [K] * [kg/m³] * [m] / [s] = [J/m²/s] => we remove ρ_soil?

		# max melt rate allowed by available energy [m/s]
        melt_rate_max = E_avail / (ρ_water * Lᵢ)

	    # actual melt rate [m/s], limited by snow available
        melt_rate = min(snow_depth[ij] / Δt, melt_rate_max)
        snow_melt_rate[ij] = max(melt_rate * ρ_water, 0)     # store to pass to soil moisture [kg/m²/s]

        # runoff rate [m/s], limited by available snow and melt + snowfall rate
        runoff_rate = r⁻¹ * snow_depth[ij]
        # maximum sink rate [m/s]
        max_sink = max(snow_depth[ij] / Δt + snow_fall_rate[ij] / ρ_water, 0)
        # limit runoff to not exceed available snow after melt and snowfall
        runoff_rate = min(runoff_rate, max(max_sink - melt_rate, 0))
        
        # store to pass to soil moisture [kg/m²/s]
        snow_runoff_rate[ij] = runoff_rate * ρ_water

        # change snow depth by falling snow minus melting and runoff [m/s]
        dsnow = snow_fall_rate[ij] / ρ_water - melt_rate  - runoff_rate

        # Euler forward time step but cap at 0 depth to not melt more snow than available
        snow_depth[ij] = max(snow_depth[ij] + Δt * dsnow, 0)
    end
end