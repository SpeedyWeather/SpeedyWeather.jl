"""
Struct holding the parameters needed at runtime in number format NF.
"""
Base.@kwdef struct DynamicsConstants{NF<:AbstractFloat} <: AbstractDynamicsConstants{NF}
    # PHYSICAL CONSTANTS
    radius::NF              # Radius of Planet
    rotation::NF            # Angular frequency of Planet's rotation
    gravity::NF             # Gravitational acceleration
    R_dry::NF               # specific gas constant for dry air [J/kg/K]
    layer_thickness::NF     # shallow water layer thickness [m]
    μ_virt_temp::NF         # used for virt temp calculation
    κ::NF                   # = R_dry/cₚ, gas const for air over heat capacity

    # CORIOLIS FREQUENCY (scaled by radius as is vorticity) = 2Ω*sin(lat)*radius
    f_coriolis::Vector{NF}

    
end

"""
Generator function for a DynamicsConstants struct.
"""
function DynamicsConstants( spectral_grid::SpectralGrid{NF},
                            planet::Planet,
                            atmosphere::AbstractAtmosphere,
                            geometry::Geometry) where NF

    # PHYSICAL CONSTANTS
    (;R_dry, R_vapour, cₚ) = atmosphere
    (;radius) = spectral_grid
    (;rotation, gravity) = planet
    layer_thickness = planet.layer_thickness*1000   # ShallowWater: convert from [km]s to [m]
    ξ = R_dry/R_vapour              # Ratio of gas constants: dry air / water vapour [1]
    μ_virt_temp = (1-ξ)/ξ           # used in Tv = T(1+μq), for conversion from humidity q
                                    # and temperature T to virtual temperature Tv
    κ = R_dry/cₚ                    # = 2/7ish for diatomic gas

    # CORIOLIS FREQUENCY (scaled by radius as is vorticity)
    f_coriolis = 2rotation*geometry.sinlat*radius

    # This implies conversion to NF
    return DynamicsConstants{NF}(   radius,rotation,gravity,R_dry,
                                    layer_thickness,μ_virt_temp,κ,
                                    f_coriolis)
end

    # # VERTICAL REFERENCE TEMPERATURE PROFILE
    # # integrate hydrostatic equation from pₛ to p, use ideal gas law p = ρRT and linear
    # # temperature decrease with height: T = Tₛ - ΔzΓ with lapse rate Γ
    # # for stratosphere (σ < σ_tropopause) increase temperature (Jablonowski & Williamson. 2006, eq. 5)
    # RΓg⁻¹ = R_dry*lapse_rate/(1000*gravity)     # convert lapse rate from [K/km] to [K/m]
    # temp_ref_profile = [temp_ref*σ^RΓg⁻¹ for σ in σ_levels_full]
    # temp_ref_profile .+= [σ < σ_tropopause ? ΔT_stratosphere*(σ_tropopause-σ)^5 : 0 for σ in σ_levels_full]

        # # GEOPOTENTIAL, coefficients to calculate geopotential
        # Δp_geopot_half, Δp_geopot_full = initialize_geopotential(σ_levels_full,σ_levels_half,R_dry)

        # # LAPSE RATE correction
        # lapserate_corr = lapserate_correction(σ_levels_full,σ_levels_half,Δp_geopot_full)
        # TROPOPAUSE/STRATOSPHERIC LEVELS
        n_stratosphere_levels = sum(σ_levels_full .< σ_tropopause)  # of levels above σ_tropopause

        # ADIABATIC TERM, Simmons and Burridge, 1981, eq. 3.12
    # precompute ln(σ_k+1/2) - ln(σ_k-1/2) but swap sign, include 1/Δσₖ
    σ_lnp_A = log.(σ_levels_half[1:end-1]./σ_levels_half[2:end]) ./ σ_levels_thick
    σ_lnp_A[1] = 0  # the corresponding sum is 1:k-1 so 0 to replace log(0) from above
    
    # precompute the αₖ = 1 - p_k-1/2/Δpₖ*log(p_k+1/2/p_k-1/2) term in σ coordinates
    σ_lnp_B = 1 .- σ_levels_half[1:end-1]./σ_levels_thick .*
                    log.(σ_levels_half[2:end]./σ_levels_half[1:end-1])
    σ_lnp_B[1] = σ_levels_half[1] <= 0 ? log(2) : σ_lnp_B[1]    # set α₁ = log(2), eq. 3.19
    σ_lnp_B .*= -1  # absorb sign from -1/Δσₖ only, eq. 3.12
