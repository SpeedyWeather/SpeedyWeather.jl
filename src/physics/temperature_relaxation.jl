struct NoTemperatureRelaxation{NF} <: TemperatureRelaxation{NF} end
NoTemperatureRelaxation(SG::SpectralGrid) = NoTemperatureRelaxation{SG.NF}()

"""$(TYPEDSIGNATURES) just passes."""
function temperature_relaxation!(   column::ColumnVariables,
                                    scheme::NoTemperatureRelaxation)
    return nothing
end

"""$(TYPEDSIGNATURES) just passes, does not need any initialization."""
function initialize!(   scheme::NoTemperatureRelaxation,
                        model::PrimitiveEquation)
    return nothing
end

"""
Struct that defines the temperature relaxation from Held and Suarez, 1996 BAMS
$(TYPEDFIELDS)"""
Base.@kwdef struct HeldSuarez{NF<:AbstractFloat} <: TemperatureRelaxation{NF}
    # DIMENSIONS
    "number of latitude rings"
    nlat::Int

    "number of vertical levels"
    nlev::Int
    
    # OPTIONS
    "sigma coordinate below which faster surface relaxation is applied"
    σb::Float64 = 0.7

    "time scale [hrs] for slow global relaxation"
    relax_time_slow::Float64 = 40*24

    "time scale [hrs] for faster tropical surface relaxation"
    relax_time_fast::Float64 = 4*24

    "minimum equilibrium temperature [K]"
    Tmin::Float64 = 200    

    "maximum equilibrium temperature [K]"
    Tmax::Float64 = 315    

    "meridional temperature gradient [K]"
    ΔTy::Float64 = 60
    
    "vertical temperature gradient [K]"
    Δθz::Float64 = 10

    # precomputed constants, allocate here, fill in initialize!
    κ::Base.RefValue{NF} = Ref(zero(NF))
    p₀::Base.RefValue{NF} = Ref(zero(NF))

    temp_relax_freq::Matrix{NF} = zeros(NF,nlev,nlat)   # (inverse) relax time scale per layer and lat
    temp_equil_a::Vector{NF} = zeros(NF,nlat)           # terms to calc equilibrium temper func
    temp_equil_b::Vector{NF} = zeros(NF,nlat)           # of latitude and pressure
end

"""
$(TYPEDSIGNATURES)
create a HeldSuarez temperature relaxation with arrays allocated given `spectral_grid`"""
function HeldSuarez(SG::SpectralGrid;kwargs...) 
    (;NF, Grid, nlat_half, nlev) = SG
    nlat = RingGrids.get_nlat(Grid,nlat_half)
    return HeldSuarez{NF}(;nlev,nlat,kwargs...)
end

"""$(TYPEDSIGNATURES)
initialize the HeldSuarez temperature relaxation by precomputing terms for the
equilibrium temperature Teq."""
function initialize!(   scheme::HeldSuarez,
                        model::PrimitiveEquation)

    (;σ_levels_full, radius, coslat, sinlat) = model.geometry
    (;σb, ΔTy, Δθz, relax_time_slow, relax_time_fast, Tmax) = scheme
    (;temp_relax_freq, temp_equil_a, temp_equil_b) = scheme
    
    p₀ = model.atmosphere.pres_ref*100                      # [hPa] → [Pa]
    scheme.p₀[] = p₀
    scheme.κ[] = model.constants.κ                          # thermodynamic kappa

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

"""$(TYPEDSIGNATURES)
Apply temperature relaxation following Held and Suarez 1996, BAMS."""
function temperature_relaxation!(   column::ColumnVariables{NF},
                                    scheme::HeldSuarez) where NF

    (;temp, temp_tend, pres, ln_pres) = column
    j = column.jring[]                      # latitude ring index j

    (;temp_relax_freq, temp_equil_a, temp_equil_b) = scheme
    Tmin = convert(NF,scheme.Tmin)
    
    p₀ = scheme.p₀[]                        # reference surface pressure
    κ = scheme.κ[]                          # thermodynamic kappa

    @inbounds for k in eachlayer(column)
        lnp = ln_pres[k]                    # logarithm of pressure at level k
        kₜ = temp_relax_freq[k,j]           # (inverse) relaxation time scale

        # Held and Suarez 1996, equation 3 with precomputed a,b during initilisation
        Teq = max(Tmin,(temp_equil_a[j] + temp_equil_b[j]*lnp)*(pres[k]/p₀)^κ)
        temp_tend[k] -= kₜ*(temp[k] - Teq)  # Held and Suarez 1996, equation 2
    end
end

"""
Struct that defines the temperature relaxation from Held and Suarez, 1996 BAMS
$(TYPEDFIELDS)"""
Base.@kwdef struct MyTemperatureRelaxation{NF<:AbstractFloat} <: TemperatureRelaxation{NF}
    # DIMENSIONS
    "number of latitude rings"
    nlat::Int

    "number of vertical levels"
    nlev::Int
    
    relax_time::Float64 = 40*24

    temp_relax_freq::Base.RefValue{NF} = Ref(zero(NF))   # (inverse) relax time scale per layer and lat
    temp_equil_a::Vector{NF} = zeros(NF,nlat)   # terms to calc equilibrium temper func
    temp_equil_s::Vector{NF} = zeros(NF,nlat)      # of latitude and pressure
end

"""
$(TYPEDSIGNATURES)
create a HeldSuarez temperature relaxation with arrays allocated given `spectral_grid`"""
function MyTemperatureRelaxation(SG::SpectralGrid;kwargs...) 
    (;NF, Grid, nlat_half, nlev) = SG
    nlat = RingGrids.get_nlat(Grid,nlat_half)
    return HeldSuarez{NF}(;nlev,nlat,kwargs...)
end

"""$(TYPEDSIGNATURES)
initialize the HeldSuarez temperature relaxation by precomputing terms for the
equilibrium temperature Teq."""
function initialize!(   scheme::MyTemperatureRelaxation,
                        model::PrimitiveEquation)
                        
    scheme.temp_relax_freq[] = radius/(3600 * scheme.relax_time)
    
    return nothing
end

"""$(TYPEDSIGNATURES)
Apply temperature relaxation following Held and Suarez 1996, BAMS."""
function temperature_relaxation!(   column::ColumnVariables{NF},
                                    scheme::MyTemperatureRelaxation) where NF

    (;temp, temp_tend) = column
    (;temp_relax_freq, temp_equil_a, temp_equil_s) = scheme
    
    @inbounds for k in eachlayer(column)
        kₜ = temp_relax_freq          # (inverse) relaxation time scale
        temp_tend[k] -= kₜ*(temp[k] - temp_equil_a[j])  # Held and Suarez 1996, equation 2
    end
    temp_tend[end] += kₜ*(temp[end] - temp_equil_a[j])  # Held and Suarez 1996, equation 2
    temp_tend[end] -= kₜ*(temp[end] - temp_equil_s[j])  # Held and Suarez 1996, equation 2
end

"""$(TYPEDSIGNATURES)
HeldSuarez-like temperature relaxation, but towards the Jablonowski temperature
profile with increasing temperatures in the stratosphere."""
Base.@kwdef struct JablonowskiRelaxation{NF<:AbstractFloat} <: TemperatureRelaxation{NF}
    # DIMENSIONS
    nlat::Int
    nlev::Int

    # OPTIONS
    "sigma coordinate below which relax_time_fast is applied"
    σb::Float64= 0.7

    "conversion from σ to Jablonowski's ηᵥ-coordinates"
    η₀::Float64 = 0.252

    "max amplitude of zonal wind [m/s]"
    u₀::Float64 = 35

    "temperature difference used for stratospheric lapse rate [K]"
    ΔT::Float64 = 4.8e5

    "[hours] time scale for slow global relaxation"
    relax_time_slow::NF = 40*24
    
    "[hours] time scale for fast aster tropical surface relaxation"
    relax_time_fast::NF = 4*24

    # precomputed constants, allocate here, fill in initialize!
    temp_relax_freq::Matrix{NF} = zeros(NF,nlev,nlat)   # (inverse) relax time scale per layer and lat
    temp_equil::Matrix{NF} = zeros(NF,nlev,nlat)        # terms to calc equilibrium temperature as func
end

"""
$(TYPEDSIGNATURES)
create a JablonowskiRelaxation temperature relaxation with arrays allocated given `spectral_grid`"""
function JablonowskiRelaxation(SG::SpectralGrid;kwargs...) 
    (;NF, Grid, nlat_half, nlev) = SG
    nlat = RingGrids.get_nlat(Grid,nlat_half)
    return JablonowskiRelaxation{NF}(;nlev,nlat,kwargs...)
end

"""$(TYPEDSIGNATURES)
initialize the JablonowskiRelaxation temperature relaxation by precomputing terms for the
equilibrium temperature Teq and the frequency (strength of relaxation)."""
function initialize!(   scheme::JablonowskiRelaxation,
                        model::PrimitiveEquation)

    (;σ_levels_full, radius, coslat, sinlat) = model.geometry
    (;σb, relax_time_slow, relax_time_fast, η₀, u₀, ΔT) = scheme
    (;temp_relax_freq, temp_equil) = scheme
    (;gravity, rotation) = model.planet
    (;lapse_rate, R_dry, σ_tropopause, temp_ref) = model.atmosphere

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
            temp_equil[k,j] = Tη + A1*((-2sinϕ^6*(cosϕ^2 + 1/3) + 10/63)*A2 +
                                            (8/5*cosϕ^3*(sinϕ^2 + 2/3) - π/4)*aΩ)
        end
    end
end 

"""$(TYPEDSIGNATURES)
Apply HeldSuarez-like temperature relaxation to the Jablonowski and Williamson
vertical profile."""
function temperature_relaxation!(   column::ColumnVariables,
                                    scheme::JablonowskiRelaxation)

    (;temp, temp_tend) = column
    j = column.jring[]                      # latitude ring index j
    (;temp_relax_freq, temp_equil) = scheme

    @inbounds for k in eachlayer(column)
        kₜ = temp_relax_freq[k,j]           # (inverse) relaxation time scale

        # Held and Suarez 1996, equation 2, but using temp_equil from
        # Jablonowski and Williamson 2006, equation 6
        temp_tend[k] -= kₜ*(temp[k] - temp_equil[k,j])  
    end
end
