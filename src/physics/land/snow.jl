abstract type AbstractSnow <: AbstractParameterization end

export SnowModel    # maybe change for a more concise name later
@kwdef mutable struct SnowModel{NF} <: AbstractSnow
	melting_threshold::NF = 275
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
	(; soil_temperature, soil_moisture) = progn.land
    
	(; mask) = model.land_sea_mask

    # Some thermodynamics needed by snow

	soil_density=snow.soil_density

	cᵢ = model.land.thermodynamics.heat_capacity_snow
	cₛ = model.land.thermodynamics.heat_capacity_dry_soil
	Lᵢ = model.clausius_clapeyron.latent_heat_fusion                  # latent heat of fusion

	# 
	S = diagn.physics.snow_rate   								# Snowfall rate in [m/s]



    launch!(architecture(snow_depth), LinearWorkOrder, size(snow_depth),
        land_snow_kernel!, snow_depth, mask, p1)
end

@kernel inbounds=true function land_snow_kernel!(
    snow_depth, mask, @Const(p1),
)
    ij = @index(Global, Linear)             # every grid point ij

    if mask[ij] > 0                         # at least partially land

		# check for melting of snow if temperature above melting threshold
		δT_melt = max(soil_temperature[ij, 1] - melting_threshold, 0)
	
		# energy available from soil warming above melting threshold [J/m²/s]
		E_avail = cₛ * δT_melt * ρ_soil * z₁ / Δt  

		# max melt rate allowed by available energy [m/s]
		melt_rate = E_avail / (ρ_w * Lᵢ)  

		# add the melt tendency to falling snow
		dsnow = melt_rate + S[ij]

		# Euler forward time step but cap at 0 depth to not melt more snow than available
		snow_depth[ij] = max(snow_depth[ij] + Δt * dsnow, 0)
    end
end
