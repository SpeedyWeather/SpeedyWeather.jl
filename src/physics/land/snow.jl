abstract type AbstractSnow <: AbstractParameterization end

export SnowModel    # maybe change for a more concise name later
@kwdef mutable struct SnowModel{NF} <: AbstractSnow
	melting_threshold::NF = 275
	leakage_time_scale::Second = Year(1)
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
	(; snow_depth) = progn.land								# in equivalent iquid water height [m]
	(; soil_temperature) = progn.land
    
	(; mask) = model.land_sea_mask

    # Some thermodynamics needed by snow
	ρ_soil = model.land.thermodynamics.soil_density			# soil density [kg/m³]
	ρ_water = model.atmosphere.water_density				# water density [kg/m³]
    Lᵢ = model.clausius_clapeyron.latent_heat_fusion      	# latent heat of fusion
	cₛ = model.land.thermodynamics.heat_capacity_dry_soil
	z₁ = model.land.geometry.layer_thickness[1]
	(; melting_threshold) = snow
	r⁻¹ = 1 / Second(snow.leakage_time_scale).value

	# Snowfall rate in [kg/m²/s]
	snow_fall_rate = diagn.physics.snow_rate
	snow_melt_rate = diagn.physics.land.snow_melt_rate

    launch!(architecture(snow_depth), LinearWorkOrder, size(snow_depth), land_snow_kernel!,
		snow_depth, soil_temperature, snow_melt_rate, snow_fall_rate, mask,
		melting_threshold, cₛ, ρ_soil, z₁, Δt,
		ρ_water, Lᵢ, r⁻¹
		)
	synchronize(architecture(snow_depth))
end

@kernel inbounds=true function land_snow_kernel!(
    snow_depth, soil_temperature, snow_melt_rate, snow_fall_rate, mask,
	@Const(melting_threshold), @Const(cₛ), @Const(ρ_soil), @Const(z₁), @Const(Δt),
	@Const(ρ_water), @Const(Lᵢ), @Const(r⁻¹),
)
    ij = @index(Global, Linear)             # every grid point ij

    if mask[ij] > 0                         # at least partially land

		# check for melting of snow if temperature above melting threshold
		δT_melt = max(soil_temperature[ij, 1] - melting_threshold, 0)
	
		# energy available from soil warming above melting threshold [J/m²/s]
		E_avail = cₛ * δT_melt * ρ_soil * z₁ / Δt  

		# max melt rate allowed by available energy [m/s]
		melt_rate_max = E_avail / (ρ_water * Lᵢ) + r⁻¹*snow_depth[ij]   # plus leakage term

	    # actual melt rate [m/s], limited by snow available
	    melt_rate = min(snow_depth[ij] / Δt, melt_rate_max)
		snow_melt_rate[ij] = melt_rate*ρ_water		# store to pass to soil moisture [kg/m²/s]

		# change snow depth by falling snow minus melting [m/s]
		dsnow = snow_fall_rate[ij]/ρ_water - melt_rate 

		# Euler forward time step but cap at 0 depth to not melt more snow than available
		snow_depth[ij] = max(snow_depth[ij] + Δt * dsnow, 0)
    end
end
