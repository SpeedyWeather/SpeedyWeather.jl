# function barrier for all AbstractCondensation
function large_scale_condensation!( 
    column::ColumnVariables,
    model::PrimitiveWet,
)
    saturation_humidity!(column, model.clausius_clapeyron)
    large_scale_condensation!(column,model.large_scale_condensation,model)
end

"""
Large scale condensation as with implicit precipitation.
$(TYPEDFIELDS)"""
Base.@kwdef struct ImplicitCondensation{NF<:AbstractFloat} <: AbstractCondensation{NF}
    "Flux limiter for latent heat release [K] per timestep"
    max_heating::NF = 0.2

    "Time scale in multiples of time step Δt"
    time_scale::NF = 7
end

ImplicitCondensation(SG::SpectralGrid;kwargs...) = ImplicitCondensation{SG.NF}(;kwargs...)

initialize!(scheme::ImplicitCondensation,model::PrimitiveEquation) = nothing

# function barrier for ImplicitCondensation to unpack model
function large_scale_condensation!( 
    column::ColumnVariables,
    scheme::ImplicitCondensation,
    model::PrimitiveWet,
)
    large_scale_condensation!(column,scheme,
        model.clausius_clapeyron,model.geometry,model.constants,model.time_stepping)
end

"""
$(TYPEDSIGNATURES)
Large-scale condensation for a `column` by relaxation back to 100%
relative humidity. Calculates the tendencies for specific humidity
and temperature from latent heat release and integrates the
large-scale precipitation vertically for output."""
function large_scale_condensation!(
    column::ColumnVariables{NF},
    scheme::ImplicitCondensation,
    clausius_clapeyron::AbstractClausiusClapeyron,
    geometry::Geometry,
    constants::DynamicsConstants,
    time_stepping::TimeStepper,
) where NF

    (;temp, humid, pres) = column           # prognostic vars: specific humidity, pressure
    (;temp_tend, humid_tend) = column       # tendencies to write into
    (;sat_humid) = column                   # intermediate variable, calculated in thermodynamics!
    
    # precompute scaling constant for precipitation output
    pₛ = pres[end]                          # surface pressure
    (;gravity, water_density) = constants
    (;Δt_sec) = time_stepping
    Δσ = geometry.σ_levels_thick
    pₛΔt_gρ = pₛ*Δt_sec/(gravity*water_density)

    (;Lᵥ, cₚ, Lᵥ_Rᵥ) = clausius_clapeyron
    Lᵥ_cₚ = Lᵥ/cₚ                           # latent heat of vaporization over heat capacity
    max_heating = scheme.max_heating/Δt_sec
    time_scale = scheme.time_scale

    @inbounds for k in eachindex(column)
        if humid[k] > sat_humid[k]

            # tendency for Implicit humid = sat_humid, divide by leapfrog time step below
            δq = sat_humid[k] - humid[k]

            # implicit correction, Frierson et al. 2006 eq. (21)
            dqsat_dT = sat_humid[k] * Lᵥ_Rᵥ/temp[k]^2       # derivative of qsat wrt to temp
            δq /= ((1 + Lᵥ_cₚ*dqsat_dT) * time_scale*Δt_sec) 

            # latent heat release with maximum heating limiter for stability
            δT = min(max_heating, -Lᵥ_cₚ * δq)
            δq = -δT/Lᵥ_cₚ                                  # also limit drying for enthalpy conservation

            # If there is large-scale condensation at a level higher (i.e. smaller k) than
            # the cloud-top previously diagnosed due to convection, then increase the cloud-top
            column.cloud_top = min(column.cloud_top, k)     # Page 7 (last sentence)
    
            # 2. Precipitation due to large-scale condensation [kg/m²/s] /ρ for [m/s]
            # += for vertical integral
            ΔpₖΔt_gρ = Δσ[k] * pₛΔt_gρ                      # Formula 4 *Δt for [m] of rain during Δt
            column.precip_large_scale -= ΔpₖΔt_gρ * δq      # Formula 25, unit [m]

            # only accumulate into humid_tend now to allow humid_tend != 0 before this scheme is called
            humid_tend[k] += δq
            temp_tend[k] += δT
        end
    end
end