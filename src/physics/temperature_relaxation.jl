struct NoTemperatureRelaxation{NF} <: TemperatureRelaxation{NF} end

"""NoBoundaryLayer scheme just passes."""
function temperature_relaxation!(   column::ColumnVariables,
                                    scheme::NoTemperatureRelaxation,
                                    model::PrimitiveEquation)
    return nothing
end

"""NoBoundaryLayer scheme does not need any initialisation."""
function initialize_temperature_relaxation!(K::ParameterizationConstants,
                                            scheme::NoTemperatureRelaxation,
                                            P::Parameters,
                                            G::Geometry)
    return nothing
end

@kwdef struct HeldSuarez{NF<:Real} <: TemperatureRelaxation{NF}
    σb::NF = 0.7                    # sigma coordinate below which linear drag is applied

    relax_time_slow::NF = 40.0*24   # [hours] time scale for slow global relaxation
    relax_time_fast::NF = 4.0*24    # [hours] time scale for faster tropical surface relaxation

    Tmin::NF = 200.0    # minimum temperature [K] in equilibrium temperature
    Tmax::NF = 315.0    # maximum temperature [K] in equilibrium temperature

    ΔTy::NF = 60.0      # meridional temperature gradient [K]
    Δθz::NF = 10.0      # vertical temperature gradient [K]
end

# generator so that HeldSuarez(Tmin=200::Int) is still possible → Float64
HeldSuarez(;kwargs...) = HeldSuarez{Float64}(;kwargs...)

function temperature_relaxation!(   column::ColumnVariables{NF},
                                    scheme::HeldSuarez,
                                    model::PrimitiveEquation) where NF

    (;temp,temp_tend,pres,ln_pres) = column
    j = column.jring[]                      # latitude ring index j
    (;temp_relax_freq,temp_equil_a,temp_equil_b) = model.parameterization_constants
    Tmin = convert(NF,scheme.Tmin)
    p₀ = convert(NF,model.parameters.pres_ref*100)  # [hPa] → [Pa]
    (;κ) = model.constants                  # R/cₚ = 2/7 

    @inbounds for k in eachlayer(column)
        lnp = ln_pres[k]                    # logarithm of pressure at level k
        kₜ = temp_relax_freq[k,j]           # (inverse) relaxation time scale

        # Held and Suarez 1996, equation 3 with precomputed a,b during initilisation
        Teq = max(Tmin,(temp_equil_a[j] + temp_equil_b[j]*lnp)*(pres[k]/p₀)^κ)
        temp_tend[k] -= kₜ*(temp[k] - Teq)  # Held and Suarez 1996, equation 2
    end
end

function initialize_temperature_relaxation!(K::ParameterizationConstants,
                                            scheme::HeldSuarez,
                                            P::Parameters,
                                            G::Geometry)

    (;σ_levels_full,radius,coslat,sinlat) = G
    (;σb,ΔTy,Δθz,relax_time_slow,relax_time_fast,Tmax) = scheme
    (;temp_relax_freq,temp_equil_a,temp_equil_b) = K
    p₀ = P.pres_ref*100                         # [hPa] → [Pa]

    # slow relaxation everywhere, fast in the tropics
    kₐ = radius/(relax_time_slow*3600)    # scale with radius as ∂ₜT is; hrs -> sec
    kₛ = radius/(relax_time_fast*3600)

    for (j,(cosϕ,sinϕ)) = enumerate(zip(coslat,sinlat))     # use ϕ for latitude here
        for (k,σ) in enumerate(σ_levels_full)
            # Held and Suarez equation 4
            temp_relax_freq[k,j] =  kₐ + (kₛ - kₐ)*max(0,(σ-σb)/(1-σb))*cosϕ^4
        end

        # Held and Suarez equation 3, split into max(Tmin,(a - b*ln(p))*(p/p₀)^κ)
        # precompute a,b to simplify online calculation
        temp_equil_a[j] = Tmax - ΔTy*sinϕ^2 + Δθz*log(p₀)*cosϕ^2
        temp_equil_b[j] = -Δθz*cosϕ^2
    end
end

"""
    J = JablonowskiRelaxation{NF}()

HeldSuarez-like temperature relaxation, but towards the Jablonowski temperature
profile with increasing temperatures in the stratosphere."""
@kwdef struct JablonowskiRelaxation{NF<:Real} <: TemperatureRelaxation{NF}
    σb::NF = 0.7        # sigma coordinate below which relax_time_fast is applied

    η₀::NF = 0.252      # conversion from σ to Jablonowski's ηᵥ-coordinates
    u₀::NF = 35.0       # max amplitude of zonal wind [m/s]
    ΔT::NF = 4.8e5      # temperature difference used for stratospheric lapse rate [K]

    relax_time_slow::NF = 40.0*24   # [hours] time scale for slow global relaxation
    relax_time_fast::NF = 4.0*24    # [hours] time scale for faster tropical surface relaxation
end

# generator so that arguments are converted to Float64
JablonowskiRelaxation(;kwargs...) = JablonowskiRelaxation{Float64}(;kwargs...)

function temperature_relaxation!(   column::ColumnVariables{NF},
                                    scheme::JablonowskiRelaxation,
                                    model::PrimitiveEquation) where NF

    (;temp,temp_tend) = column
    j = column.jring[]                      # latitude ring index j
    (;temp_relax_freq,temp_equil) = model.parameterization_constants

    @inbounds for k in eachlayer(column)
        kₜ = temp_relax_freq[k,j]           # (inverse) relaxation time scale

        # Held and Suarez 1996, equation 2, but using temp_equil from
        # Jablonowski and Williamson 2006, equation 6
        temp_tend[k] -= kₜ*(temp[k] - temp_equil[k,j])  
    end
end

function initialize_temperature_relaxation!(K::ParameterizationConstants,
                                            scheme::JablonowskiRelaxation,
                                            P::Parameters,
                                            G::Geometry)

    (;σ_levels_full,radius,coslat,sinlat) = G
    (;σb,relax_time_slow,relax_time_fast,η₀,u₀,ΔT) = scheme
    (;temp_relax_freq,temp_equil) = K
    (;gravity,radius,rotation) = P.planet
    (;lapse_rate,R_dry,σ_tropopause,temp_ref) = P

    Γ = lapse_rate/1000                   # from [K/km] to [K/m]
    aΩ = radius*rotation

    # slow relaxation everywhere, fast in the tropics
    kₐ = radius/(relax_time_slow*3600)    # scale with radius as ∂ₜT is; hrs -> sec
    kₛ = radius/(relax_time_fast*3600)

    for (j,(cosϕ,sinϕ)) = enumerate(zip(coslat,sinlat))     # use ϕ for latitude here
        for (k,σ) in enumerate(σ_levels_full)
            # Held and Suarez equation 4
            temp_relax_freq[k,j] =  kₐ + (kₛ - kₐ)*max(0,(σ-σb)/(1-σb))*cosϕ^4
        
            # vertical profile
            Tη = temp_ref*σ^(R_dry*Γ/gravity)    # Jablonowski and Williamson eq. 4

            if σ < σ_tropopause
                Tη += ΔT*(σ_tropopause-σ)^5      # Jablonowski and Williamson eq. 5
            end

            η = σ                   # Jablonowski and Williamson use η for σ coordinates
            ηᵥ = (η - η₀)*π/2       # auxiliary variable for vertical coordinate

            # amplitudes with height
            A1 = 3/4*η*π*u₀/R_dry*sin(ηᵥ)*sqrt(cos(ηᵥ))
            A2 = 2u₀*cos(ηᵥ)^(3/2)

            # Jablonowski and Williamson, eq. (6) 
            temp_equil[k,j] = Tη + A1*((-2sinϕ^6*(cosϕ^2 + 1/3) + 10/63)*A2 + (8/5*cosϕ^3*(sinϕ^2 + 2/3) - π/4)*aΩ)
        end
    end
end 