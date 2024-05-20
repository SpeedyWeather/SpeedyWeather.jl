abstract type AbstractTemperatureRelaxation <: AbstractParameterization end

# function barrier to unpack model.temperature_relaxation
function temperature_relaxation!(column::ColumnVariables, model::PrimitiveEquation)
    temperature_relaxation!(column, model.temperature_relaxation, model)
end

export NoTemperatureRelaxation
struct NoTemperatureRelaxation <: AbstractTemperatureRelaxation end
NoTemperatureRelaxation(::SpectralGrid) = NoTemperatureRelaxation()
initialize!(::NoTemperatureRelaxation, ::PrimitiveEquation) = nothing
temperature_relaxation!(::ColumnVariables, ::NoTemperatureRelaxation, ::PrimitiveEquation) = nothing

export HeldSuarez

"""
Temperature relaxation from Held and Suarez, 1996 BAMS
$(TYPEDFIELDS)"""
Base.@kwdef struct HeldSuarez{NF<:AbstractFloat} <: AbstractTemperatureRelaxation
    # DIMENSIONS
    "number of latitude rings"
    nlat::Int

    "number of vertical levels"
    nlev::Int
    
    # OPTIONS
    "sigma coordinate below which faster surface relaxation is applied"
    σb::NF = 0.7

    "time scale for slow global relaxation"
    relax_time_slow::Second = Day(40)

    "time scale for faster tropical surface relaxation"
    relax_time_fast::Second = Day(4)

    "minimum equilibrium temperature [K]"
    Tmin::NF = 200    

    "maximum equilibrium temperature [K]"
    Tmax::NF = 315    

    "meridional temperature gradient [K]"
    ΔTy::NF = 60
    
    "vertical temperature gradient [K]"
    Δθz::NF = 10

    # precomputed constants, allocate here, fill in initialize!
    κ::Base.RefValue{NF} = Ref(zero(NF))
    p₀::Base.RefValue{NF} = Ref(zero(NF))

    temp_relax_freq::Matrix{NF} = zeros(NF, nlev, nlat) # (inverse) relax time scale per layer and lat
    temp_equil_a::Vector{NF} = zeros(NF, nlat)          # terms to calc equilibrium temper func
    temp_equil_b::Vector{NF} = zeros(NF, nlat)          # of latitude and pressure
end

"""
$(TYPEDSIGNATURES)
create a HeldSuarez temperature relaxation with arrays allocated given `spectral_grid`"""
function HeldSuarez(SG::SpectralGrid; kwargs...)
    (; NF, nlat, nlev) = SG
    return HeldSuarez{NF}(; nlev, nlat, kwargs...)
end

"""$(TYPEDSIGNATURES)
initialize the HeldSuarez temperature relaxation by precomputing terms for the
equilibrium temperature Teq."""
function initialize!(   scheme::HeldSuarez,
                        model::PrimitiveEquation)

    (; σ_levels_full, coslat, sinlat) = model.geometry
    (; σb, ΔTy, Δθz, relax_time_slow, relax_time_fast, Tmax) = scheme
    (; temp_relax_freq, temp_equil_a, temp_equil_b) = scheme
                           
    (; pres_ref) = model.atmosphere
    scheme.p₀[] = pres_ref              # surface reference pressure [Pa]
    scheme.κ[] = model.atmosphere.κ     # thermodynamic kappa R_dry/cₚ

    # slow relaxation everywhere, fast in the tropics
    kₐ = 1/relax_time_slow.value
    kₛ = 1/relax_time_fast.value

    for (j, (cosϕ, sinϕ)) = enumerate(zip(coslat, sinlat))     # use ϕ for latitude here
        for (k, σ) in enumerate(σ_levels_full)
            # Held and Suarez equation 4
            temp_relax_freq[k, j] =  kₐ + (kₛ - kₐ)*max(0, (σ-σb)/(1-σb))*cosϕ^4
        end

        # Held and Suarez equation 3, split into max(Tmin, (a - b*ln(p))*(p/p₀)^κ)
        # precompute a, b to simplify online calculation
        temp_equil_a[j] = Tmax - ΔTy*sinϕ^2 + Δθz*log(pres_ref)*cosϕ^2
        temp_equil_b[j] = -Δθz*cosϕ^2
    end
end

# function barrier
function temperature_relaxation!(
    column::ColumnVariables,
    scheme::HeldSuarez,
    model::PrimitiveEquation
)
    temperature_relaxation!(column, scheme, model.atmosphere)
end

"""$(TYPEDSIGNATURES)
Apply temperature relaxation following Held and Suarez 1996, BAMS."""
function temperature_relaxation!(   
    column::ColumnVariables,
    scheme::HeldSuarez,
    atmosphere::AbstractAtmosphere,
)
    (; temp, temp_tend, pres, ln_pres) = column
    j = column.jring[]                      # latitude ring index j
    (; Tmin, temp_relax_freq, temp_equil_a, temp_equil_b) = scheme

    # surface reference pressure [Pa] and thermodynamic kappa R_dry/cₚ
    (; pres_ref, κ) = atmosphere

    @inbounds for k in eachlayer(column)
        lnp = ln_pres[k]                    # logarithm of pressure at level k
        kₜ = temp_relax_freq[k, j]           # (inverse) relaxation time scale

        # Held and Suarez 1996, equation 3 with precomputed a, b during initilisation
        Teq = max(Tmin, (temp_equil_a[j] + temp_equil_b[j]*lnp)*(pres[k]/pres_ref)^κ)
        temp_tend[k] -= kₜ*(temp[k] - Teq)  # Held and Suarez 1996, equation 2
    end
end

export JablonowskiRelaxation

"""$(TYPEDSIGNATURES)
HeldSuarez-like temperature relaxation, but towards the Jablonowski temperature
profile with increasing temperatures in the stratosphere."""
Base.@kwdef mutable struct JablonowskiRelaxation{NF<:AbstractFloat} <: AbstractTemperatureRelaxation
    
    # DIMENSIONS
    nlat::Int
    nlev::Int

    # OPTIONS
    "sigma coordinate below which relax_time_fast is applied [1]"
    σb::NF = 0.7

    "sigma coordinate for tropopause temperature inversion"
    σ_tropopause::NF = 0.2

    "conversion from σ to Jablonowski's ηᵥ-coordinates"
    η₀::NF = 0.252

    "max amplitude of zonal wind [m/s]"
    u₀::NF = 35

    "temperature difference used for stratospheric lapse rate [K]"
    ΔT::NF = 4.8e5

    "Dry-adiabatic lapse rate [K/m]"
    lapse_rate::NF = 5/1000

    "time scale for slow global relaxation"
    relax_time_slow::Second = Day(40)
    
    "time scale for faster tropical surface relaxation"
    relax_time_fast::Second = Day(4)

    # precomputed constants, allocate here, fill in initialize!
    temp_relax_freq::Matrix{NF} = zeros(NF, nlev, nlat)   # (inverse) relax time scale per layer and lat
    temp_equil::Matrix{NF} = zeros(NF, nlev, nlat)        # terms to calc equilibrium temperature as func
end

"""
$(TYPEDSIGNATURES)
create a JablonowskiRelaxation temperature relaxation with arrays allocated given `spectral_grid`"""
function JablonowskiRelaxation(SG::SpectralGrid; kwargs...) 
    (; NF, nlat, nlev) = SG
    return JablonowskiRelaxation{NF}(; nlev, nlat, kwargs...)
end

"""$(TYPEDSIGNATURES)
initialize the JablonowskiRelaxation temperature relaxation by precomputing terms for the
equilibrium temperature Teq and the frequency (strength of relaxation)."""
function initialize!(   scheme::JablonowskiRelaxation,
                        model::PrimitiveEquation)

    (; σ_levels_full, coslat, sinlat) = model.geometry
    (; σb, relax_time_slow, relax_time_fast, η₀, u₀, ΔT) = scheme
    (; temp_relax_freq, temp_equil, σ_tropopause, lapse_rate) = scheme
    (; gravity) = model.planet
    (; R_dry, temp_ref) = model.atmosphere
    Ω = model.planet.rotation

    # slow relaxation [1/s] everywhere, fast in the tropics
    kₐ = 1/relax_time_slow.value
    kₛ = 1/relax_time_fast.value

    for (j, (cosϕ, sinϕ)) = enumerate(zip(coslat, sinlat))     # use ϕ for latitude here
        for (k, σ) in enumerate(σ_levels_full)
            # Held and Suarez equation 4
            temp_relax_freq[k, j] =  kₐ + (kₛ - kₐ)*max(0, (σ-σb)/(1-σb))*cosϕ^4
        
            # vertical profile
            Tη = temp_ref*σ^(R_dry*lapse_rate/gravity)      # Jablonowski and Williamson eq. 4

            if σ < σ_tropopause
                Tη += ΔT*(σ_tropopause-σ)^5      # Jablonowski and Williamson eq. 5
            end

            η = σ                   # Jablonowski and Williamson use η for σ coordinates
            ηᵥ = (η - η₀)*π/2       # auxiliary variable for vertical coordinate

            # amplitudes with height
            A1 = 3/4*η*π*u₀/R_dry*sin(ηᵥ)*sqrt(cos(ηᵥ))
            A2 = 2u₀*cos(ηᵥ)^(3/2)

            # Jablonowski and Williamson, eq. (6) 
            temp_equil[k, j] = Tη + A1*((-2sinϕ^6*(cosϕ^2 + 1/3) + 10/63)*A2 +
                                            (8/5*cosϕ^3*(sinϕ^2 + 2/3) - π/4)*Ω)
        end
    end
end

# function barrier
function temperature_relaxation!(
    column::ColumnVariables,
    scheme::JablonowskiRelaxation,
    model::PrimitiveEquation,
)
    temperature_relaxation!(column, scheme)
end


"""$(TYPEDSIGNATURES)
Apply HeldSuarez-like temperature relaxation to the Jablonowski and Williamson
vertical profile."""
function temperature_relaxation!(   column::ColumnVariables,
                                    scheme::JablonowskiRelaxation)

    (; temp, temp_tend) = column
    j = column.jring[]                      # latitude ring index j
    (; temp_relax_freq, temp_equil) = scheme

    @inbounds for k in eachlayer(column)
        kₜ = temp_relax_freq[k, j]           # (inverse) relaxation time scale

        # Held and Suarez 1996, equation 2, but using temp_equil from
        # Jablonowski and Williamson 2006, equation 6
        temp_tend[k] -= kₜ*(temp[k] - temp_equil[k, j])
    end
end
