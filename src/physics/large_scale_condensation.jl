abstract type AbstractCondensation <: AbstractParameterization end

export NoCondensation
struct NoCondensation <: AbstractCondensation end
NoCondesantion(::SpectralGrid) = NoCondensation()
initialize!(::NoCondensation, ::PrimitiveEquation) = nothing
large_scale_condensation!(::ColumnVariables, ::NoCondensation, ::PrimitiveEquation) = nothing

export ImplicitCondensation
"""
Large scale condensation as with implicit precipitation.
$(TYPEDFIELDS)"""
@kwdef struct ImplicitCondensation{NF<:AbstractFloat} <: AbstractCondensation
    "Relative humidity threshold [1 = 100%] to trigger condensation"
    relative_humidity_threshold::NF = 1
    freezing_threshold::NF = 263.0

    "Time scale in multiples of time step Δt, the larger the less immediate"
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
    #TODO not needed for NoCondensation
    saturation_humidity!(column, model.clausius_clapeyron)
    large_scale_condensation!(column, model.large_scale_condensation, model)
end

# function barrier for ImplicitCondensation to unpack model
function large_scale_condensation!( 
    column::ColumnVariables,
    scheme::ImplicitCondensation,
    model::PrimitiveWet,
)
    large_scale_condensation!(column, scheme,
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
    scheme::ImplicitCondensation,
    clausius_clapeyron::AbstractClausiusClapeyron,
    geometry::Geometry,
    planet::AbstractPlanet,
    atmosphere::AbstractAtmosphere,
    time_stepping::AbstractTimeStepper,
)

    (; pres, temp, humid) = column          # prognostic vars (from previous time step for numerical stability)
    (; temp_tend, humid_tend) = column      # tendencies to write into
    (; sat_humid) = column                  # intermediate variable, calculated in thermodynamics!
    
    # precompute scaling constant for precipitation output
    pₛ = pres[end]                          # surface pressure
    (; Δt_sec) = time_stepping
    Δσ = geometry.σ_levels_thick
    pₛΔt_gρ = (pₛ * Δt_sec)/(planet.gravity * atmosphere.water_density)

    (; Lᵥ, cₚ, Lᵥ_Rᵥ) = clausius_clapeyron
    Lᵥ_cₚ = Lᵥ/cₚ                           # latent heat of vaporization over heat capacity
    (; time_scale, relative_humidity_threshold) = scheme

    @inbounds for k in eachindex(column)
        if humid[k] > sat_humid[k]*relative_humidity_threshold

            # tendency for Implicit humid = sat_humid, divide by leapfrog time step below
            δq = sat_humid[k] * relative_humidity_threshold - humid[k]

            # implicit correction, Frierson et al. 2006 eq. (21)
            dqsat_dT = sat_humid[k] * relative_humidity_threshold * Lᵥ_Rᵥ/temp[k]^2       # derivative of qsat wrt to temp
            δq /= ((1 + Lᵥ_cₚ*dqsat_dT) * time_scale*Δt_sec) 

            # latent heat release for enthalpy conservation
            δT = -Lᵥ_cₚ * δq

            # If there is large-scale condensation at a level higher (i.e. smaller k) than
            # the cloud-top previously diagnosed due to convection, then increase the cloud-top
            # Fortran SPEEDY documentation Page 7 (last sentence)
            column.cloud_top = min(column.cloud_top, k)
    
            # 2. Precipitation due to large-scale condensation [kg/m²/s] /ρ for [m/s]
            # += for vertical integral
            precip = Δσ[k] * pₛΔt_gρ * -δq          # precipitation [m] on layer k, Formula 4
            column.precip_large_scale += precip     # integrate vertically, Formula 25, unit [m]

            # only accumulate into humid_tend now to allow humid_tend != 0 before this scheme is called
            humid_tend[k] += δq
            temp_tend[k] += δT
        end
    end

    # convert to rain rate [m/s]
    column.precip_rate_large_scale = column.precip_large_scale / Δt_sec
    if temp[:] < freezing_threshold
        column.snow_rate_large_scale = column.precip_rate_large_scale
        column.precip_rate_large_scale = 0.
    return nothing
end
