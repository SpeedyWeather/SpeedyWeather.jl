"""
Large scale condensation as in Fortran SPEEDY with default values from therein.
$(TYPEDFIELDS)"""
Base.@kwdef struct SpeedyCondensation{NF<:AbstractFloat} <: AbstractCondensation{NF}

    "number of vertical levels"
    nlev::Int

    "Relative humidity threshold"
    relative_threshold::Float64 = 0.95

    "Relative humidity baseline for relaxation"
    relative_baseline::Float64 = 0.9

    "Relaxation time for humidity [hrs]"
    time_scale::Second = Hour(4)

    # precomputed arrays
    n_stratosphere_levels::Base.RefValue{Int} = Ref(0)
end

SpeedyCondensation(SG::SpectralGrid;kwargs...) = SpeedyCondensation{SG.NF}(nlev=SG.nlev;kwargs...)

"""
$(TYPEDSIGNATURES)
Initialize the SpeedyCondensation scheme."""
function initialize!(scheme::SpeedyCondensation,model::PrimitiveEquation)

    (;σ_levels_full) = model.geometry
    (;σ_tropopause) = model.atmosphere

    scheme.n_stratosphere_levels[] = findlast(σ->σ<=σ_tropopause,σ_levels_full)
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

"""Function barrier only."""
function large_scale_condensation!( 
    column::ColumnVariables,
    model::PrimitiveWet,
)
    large_scale_condensation!(column,model.large_scale_condensation,
        model.geometry,model.constants,model.atmosphere,model.time_stepping)
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
    atmosphere::AbstractAtmosphere,
    time_stepping::TimeStepper,
    ) where NF

    (;relative_threshold,relative_baseline) = scheme
    time_scale⁻¹ = scheme.time_scale.value
    n_stratosphere_levels = scheme.n_stratosphere_levels[]

    (;humid, pres) = column               # prognostic variables: specific humidity, surface pressure
    (;temp_tend, humid_tend) = column     # tendencies to write into
    (;sat_humid) = column                 # intermediate variable, calculate in thermodynamics!
    (;nlev) = column
    pₛ = pres[end]                         # surface pressure

    (;gravity, water_density) = constants
    (;Δt_sec) = time_stepping
    pₛΔt_gρ = pₛ*Δt_sec/gravity/water_density   # precompute constant

    (;σ_levels_thick) = geometry
    latent_heat = convert(NF, atmosphere.latent_heat_condensation)
    
    # 1. Tendencies of humidity and temperature due to large-scale condensation
    @inbounds for k in n_stratosphere_levels+1:nlev   # top to bottom, skip stratospheric levels

        # Specific humidity threshold for condensation
        humid_threshold = relative_threshold * sat_humid[k]               
        humid_baseline = relative_baseline * sat_humid[k]               

        if humid[k] > humid_threshold

            # accumulate in tendencies (nothing is added if humidity not above threshold)
            humid_tend_k = -(humid[k] - humid_baseline) * time_scale⁻¹                  # Formula 22
            # temp_tend[k] += -latent_heat * min(humid_tend_k, humid_tend_max[k]*pres[k]) # Formula 23
            temp_tend[k] += -latent_heat * humid_tend_k             # Formula 23, without limiter

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