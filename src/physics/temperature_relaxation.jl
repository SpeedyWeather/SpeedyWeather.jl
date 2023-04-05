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