"""
Large scale condensation as in Fortran SPEEDY with default values from therein.
$(TYPEDFIELDS)"""
Base.@kwdef struct SpeedyCondensation{NF<:AbstractFloat} <: AbstractCondensation{NF}

    "number of vertical levels"
    nlev::Int

    "Relative humidity threshold for boundary layer"
    threshold_boundary_layer::NF = 0.95

    "Vertical range of relative humidity threshold"
    threshold_range::NF = 0.1

    "Maximum relative humidity threshold [1]"
    threshold_max::NF = 0.9

    "Relaxation time for humidity [hrs]"
    time_scale::Second = Hour(4)

    "Flux limiter for humidity tendency"
    max_flux::NF = 0.1

    "Flux limiter for condensation heating"
    max_heating::NF = 0.001

    # precomputed arrays
    n_stratosphere_levels::Base.RefValue{Int} = Ref(0)
    humid_tend_max::Vector{NF} = zeros(NF,nlev)
    relative_threshold::Vector{NF} = zeros(NF,nlev)
    latent_heat_cₚ::Base.RefValue{NF} = Ref(zero(NF))
end

SpeedyCondensation(SG::SpectralGrid;kwargs...) = SpeedyCondensation{SG.NF}(nlev=SG.nlev;kwargs...)

"""
$(TYPEDSIGNATURES)
Initialize the SpeedyCondensation scheme."""
function initialize!(scheme::SpeedyCondensation,model::PrimitiveEquation)

    (;nlev,threshold_boundary_layer,threshold_range,threshold_max) = scheme
    (;humid_tend_max, max_flux, relative_threshold, time_scale) = scheme
    (;σ_levels_full) = model.geometry
    (;σ_tropopause,σ_boundary_layer) = model.atmosphere
   
    n = findlast(σ->σ<=σ_tropopause,σ_levels_full)
    scheme.n_stratosphere_levels[] = isnothing(n) ? 0 : n

    for k in 1:nlev
        # the relative humidity threshold above which condensation occurs per layer
        σₖ² = σ_levels_full[k]^2
        relative_threshold[k] = threshold_max + threshold_range * (σₖ² - 1)
        if σ_levels_full[k] >= σ_boundary_layer
            relative_threshold[k] = max(relative_threshold[k], threshold_boundary_layer)
        end

        # Impose a maximum heating rate to avoid grid-point storm instability
        # This formula does not appear in the original SPEEDY documentation
        # there's a (pres[end]/pres_ref)^2 which we set to 1 for simplicity
        humid_tend_max[k] = max_flux*σₖ² / time_scale.value
    end

    # to convert humidity tendency to temperature tendency
    # Fortran SPEEDY documentation eq. (23)
    scheme.latent_heat_cₚ[] = model.atmosphere.latent_heat_condensation/    # [J/kg]
                                model.atmosphere.cₚ     # [J/K/kg] specific heat capacity, const pres
end

"""
$(TYPEDSIGNATURES)
No condensation in a PrimitiveDry model."""
function large_scale_condensation!( 
    column::ColumnVariables,
    model::PrimitiveDry,
)
    return nothing
end

# function barrier for all AbstractCondensation
function large_scale_condensation!( 
    column::ColumnVariables,
    model::PrimitiveWet,
)
    saturation_humidity!(column, model.clausis_clapeyron)
    large_scale_condensation!(column,model.large_scale_condensation,model)
end

# function barrier for SpeedyCondensation to unpack model
function large_scale_condensation!( 
    column::ColumnVariables,
    scheme::SpeedyCondensation,
    model::PrimitiveWet,
)
    large_scale_condensation!(column,scheme,
        model.geometry,model.constants,model.time_stepping)
end

"""
$(TYPEDSIGNATURES)
Large-scale condensation for a `column` by relaxation back to a reference
relative humidity if larger than that. Calculates the tendencies for
specific humidity and temperature and integrates the large-scale
precipitation vertically for output."""
function large_scale_condensation!(
    column::ColumnVariables{NF},
    scheme::SpeedyCondensation,
    geometry::Geometry,
    constants::DynamicsConstants,
    time_stepping::TimeStepper,
    ) where NF

    (;relative_threshold, humid_tend_max) = scheme
    time_scale⁻¹ = inv(convert(NF,scheme.time_scale.value))
    n_stratosphere_levels = scheme.n_stratosphere_levels[]

    (;humid, pres) = column                 # prognostic vars: specific humidity, pressure
    (;temp_tend, humid_tend) = column       # tendencies to write into
    (;sat_humid) = column                   # intermediate variable, calculated in thermodynamics!
    (;nlev) = column
    pₛ = pres[end]                          # surface pressure
    pₛ_norm² = (pₛ/constants.pres_ref)^2
    max_heating = -scheme.max_heating*time_scale⁻¹

    # precompute scaling constant, undo radius scaling here as directly used for output
    (;gravity, water_density) = constants
    (;Δt_sec) = time_stepping
    (;σ_levels_thick) = geometry
    pₛΔt_gρ = pₛ*Δt_sec/gravity/water_density 

    # 1. Tendencies of humidity and temperature due to large-scale condensation
    @inbounds for k in n_stratosphere_levels+1:nlev   # top to bottom, skip stratospheric levels

        # Specific humidity threshold for condensation
        humid_threshold = relative_threshold[k] * sat_humid[k]                           

        if humid[k] > humid_threshold

            # accumulate in tendencies (nothing is added if humidity not above threshold)
            humid_tend_k = -(humid[k] - humid_threshold) * time_scale⁻¹                   # eq. 22

            # with flux limiter, use max and - as humid_tend_k is always negative
            humid_tend_k = max(humid_tend_k, -pₛ_norm²*humid_tend_max[k])

            # condensation heating, eq. 23 
            temp_tend[k] -= max(max_heating, scheme.latent_heat_cₚ[] * humid_tend_k)

            # If there is large-scale condensation at a level higher (i.e. smaller k) than
            # the cloud-top previously diagnosed due to convection, then increase the cloud-top
            column.cloud_top = min(column.cloud_top, k)             # Page 7 (last sentence)
    
            # 2. Precipitation due to large-scale condensation [kg/m²/s] /ρ for [m/s]
            # += for vertical integral
            ΔpₖΔt_gρ = σ_levels_thick[k]*pₛΔt_gρ                    # Formula 4 *Δt for [m] of rain during Δt
            column.precip_large_scale += -ΔpₖΔt_gρ * humid_tend_k   # Formula 25, unit [m]

            # only accumulate into humid_tend now to allow humid_tend != 0 before this scheme is called
            humid_tend[k] += humid_tend_k
        end
    end
end

"""
Large scale condensation as with immediate precipitation.
$(TYPEDFIELDS)"""
Base.@kwdef struct ImmediateCondensation{NF<:AbstractFloat} <: AbstractCondensation{NF}
    "Flux limiter for latent heat release [K] per timestep"
    max_heating::NF = 0.2
end

ImmediateCondensation(SG::SpectralGrid;kwargs...) = ImmediateCondensation{SG.NF}(;kwargs...)

initialize!(scheme::ImmediateCondensation,model::PrimitiveEquation) = nothing

# function barrier for ImmediateCondensation to unpack model
function large_scale_condensation!( 
    column::ColumnVariables,
    scheme::ImmediateCondensation,
    model::PrimitiveWet,
)
    large_scale_condensation!(column,scheme,
        model.clausis_clapeyron,model.geometry,model.constants,model.time_stepping)
end

"""
$(TYPEDSIGNATURES)
Large-scale condensation for a `column` by relaxation back to a reference
relative humidity if larger than that. Calculates the tendencies for
specific humidity and temperature and integrates the large-scale
precipitation vertically for output."""
function large_scale_condensation!(
    column::ColumnVariables{NF},
    scheme::ImmediateCondensation,
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
    (;σ_levels_thick) = geometry
    pₛΔt_gρ = pₛ*Δt_sec/gravity/water_density

    (;Lᵥ, cₚ, Lᵥ_Rᵥ) = clausius_clapeyron
    Lᵥ_cₚ = Lᵥ/cₚ                           # latent heat of vaporization over heat capacity
    max_heating = scheme.max_heating/Δt_sec

    @inbounds for k in eachindex(column)
        if humid[k] > sat_humid[k]

            # tendency for immediate humid = sat_humid
            humid_tend_k = (sat_humid[k] - humid[k])/(2Δt_sec)

            # implicit correction, Frierson et al. 2006 eq. (21)
            # dqsat_dT = grad_saturation_humidity(clausius_clapeyron,temp[k],pres[k])
            dqsat_dT = sat_humid[k] * Lᵥ_Rᵥ/temp[k]^2
            humid_tend_k /= (1 + Lᵥ_cₚ*dqsat_dT)

            # latent heat release with maximum heating limiter for stability
            # note that this violates enthalpy conservation
            temp_tend[k] += min(max_heating, -Lᵥ_cₚ * humid_tend_k)

            # If there is large-scale condensation at a level higher (i.e. smaller k) than
            # the cloud-top previously diagnosed due to convection, then increase the cloud-top
            column.cloud_top = min(column.cloud_top, k)             # Page 7 (last sentence)
    
            # 2. Precipitation due to large-scale condensation [kg/m²/s] /ρ for [m/s]
            # += for vertical integral
            ΔpₖΔt_gρ = σ_levels_thick[k]*pₛΔt_gρ                    # Formula 4 *Δt for [m] of rain during Δt
            column.precip_large_scale -= ΔpₖΔt_gρ * humid_tend_k    # Formula 25, unit [m]

            # only accumulate into humid_tend now to allow humid_tend != 0 before this scheme is called
            humid_tend[k] += humid_tend_k
        end
    end
end