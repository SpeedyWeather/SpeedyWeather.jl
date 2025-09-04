abstract type AbstractSnow <: AbstractParameterization end

export SnowModel    # maybe change for a more concise name later
@kwdef mutable struct SnowModel{NF} <: AbstractSnow
    parameter1::NF = 0
    parameter2::NF = 0
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
    (; snow_depth) = progn.land
    (; mask) = model.land_sea_mask
    #p1 = snow.parameter1  # ??? PLV I do not know what to make of this
    # Some thermodynamics needed by snow
	cᵢ = thermodynamics.heat_capacity_snow
	cₛ = thermodynamics.heat_capacity_soil
	Lᵢ = clausius_clapeyron.latent_heat_fusion            # latent heat of fusion
    drag_snow::NF = 3.0e-5 # [m]
    soil_density::NF = 1.3e4 #[kg/m2]
    snow_depth_scale::NF = 0.1 #[m]                       #Tunable parameter
    # Some atmospheric variables needed to compute snow sublimation
    P = diagn.physics.total_precipitation_rate      # precipitation (rain+snow) in [m/s]
	S = diagn.physics.large_scale_snow_rate +  diagn.physics.convective_snow_rate    # Snowfall rate in [m/s]
    E = diagn.physics.land.surface_humidity_flux    # [kg/s/m²], divide by density for [m/s]
    R = diagn.physics.land.river_runoff             # diagnosed [m/s]
    ρ = column.surface_air_density
    V₀ = column.surface_wind_speed
	Sₐ = progn.land.snow_depth                             # snow accumulation at the surface (m of water)
	(; surface_pressure) = column
    (; soil_temperature, soil_moisture) = progn.land
	(; land_fraction, albedo_land) = column
	snow_albedo_old=model.land.thermodynamics.snow_albedo_old
	snow_albedo_fresh=model.land.thermodynamics.snow_albedo_fresh
    Δt = model.time_stepping.Δt_sec

    # PLV's code here updating snow_depth, e.g.
    ## First, lose snow by sublimation and lose liquid precip by runoff over snow
	infiltration_flux[ij] = (P[ij] - E[ij] - S[ij])/ρ   # Need to remove snowfall from total precipitation, to use just liquid part
	if Sₐ[ij] > snow_threshold
		Rₛ[ij] = infiltration_flux[ij]             # put all precipitation into runoff (no infiltration)
		infiltration_flux[ij] = 0.                 # Assume no infiltration if snow on surface        - R[ij]
		sat_humid_land[ij] = saturation_humidity(soil_temperature[ij,1], surface_pressure[ij], model.clausius_clapeyron)
		sublimation_flux_snow[ij] = isfinite(soil_temperature[ij,1]) && isfinite(Sₐ[ij]) ?
	                ρ[ij]*drag_snow*V₀[ij]*(sat_humid_land[ij]  - surface_humid)[ij] : zero(NF)
	else
		F[ij] = infiltration_flux[ij]                       # - R[ij]
	end
	## now check for melting of snow
	δT_melt[ij] = max(soil_temperature[ij,1] - melting_threshold, 0)   # only if temperature above melting threshold
	if δT_melt[ij] > 0
	    # energy available from soil warming above melting threshold [J/m²/s]
	    E_avail[ij] = cₛ * δT_melt[ij] * ρ_soil * z₁ / Δt  

	    # max melt rate allowed by available energy [m/s]
	    melt_rate_max[ij] = E_avail[ij] / (ρ_w * Lᵢ)  

	    # actual melt rate, limited by snow available
	    melt_rate[ij] = min(Sₐ[ij] / Δt, melt_rate_max)  # [m/s]

	    # snow water equivalent melted this timestep [m]
	    melt_depth[ij] = melt_rate[ij] * Δt

	    # convert to infiltration flux
	    melt_flux[ij] = melt_depth[ij]                   # [m], can add to infiltration

	    # energy used for melting [J/m²]
	    Q_melt[ij] = ρ_w * Lᵢ * melt_depth[ij]

	    # soil heat capacity [J/m²/K]
	    C_soil = ρ_soil * cₛ * z₁

	    # temperature tendency
	    δT[ij] = -Q_melt[ij] / C_soil
	    soil_temperature[ij,1] += δT[ij]
	end
	## Now compute the snow area cover fraction based on snow depth
	#σₛ = Sₐ / (10. * drag_snow + Sₐ) #JULES, from Betts et al.
	σₛ[ij] = min(1.0, Sₐ[ij] / snow_depth_scale)   # e.g. snow_depth_scale ≈ 0.1 m
	## Now compute the change in albedo
	τˢ = 10.0 * day   # tunable parameter
	albedo_snow[ij] = snow_albedo_old + (snow_albedo_fresh - snow_albedo_old) * exp(-Δt / τˢ) # time needs to be expresssed in days. 
	albedo_land[ij] = (1-σₛ[ij]) * albedo_land + σₛ[ij] * albedo_snow[ij]

	F[ij] = infiltration_flux[ij] + melt_flux[ij]  # Need to pass this back to bucket model
    R[ij] = + Rₛ[ij]   #??? PLV this needs to be returned to the bucket model, it may be enough to do it here

    launch!(architecture(snow_depth), LinearWorkOrder, size(snow_depth),
        land_snow_kernel!, snow_depth, mask, p1)
end

@kernel inbounds=true function land_snow_kernel!(
    snow_depth, mask, @Const(p1),
)
    ij = @index(Global, Linear)             # every grid point ij

    if mask[ij] > 0                         # at least partially land
        snow_depth[ij] = 0                  # dummy operation, replace with real logic
    end
end
