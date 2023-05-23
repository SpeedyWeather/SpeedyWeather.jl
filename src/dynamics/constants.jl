"""
Struct holding constants needed at runtime for the dynamical core in number format NF.
$(TYPEDFIELDS)"""
@kwdef struct DynamicsConstants{NF<:AbstractFloat} <: AbstractDynamicsConstants{NF}
    # PHYSICAL CONSTANTS
    radius::NF              # Radius of Planet
    rotation::NF            # Angular frequency of Planet's rotation
    gravity::NF             # Gravitational acceleration
    layer_thickness::NF     # shallow water layer thickness [m]
    
    # THERMODYNAMICS
    R_dry::NF               # specific gas constant for dry air [J/kg/K]
    R_vapour::NF            # specific gas constant for water vapour [J/kg/K]
    μ_virt_temp::NF         # used for virt temp calculation
    cₚ::NF                  # specific heat at constant pressure [J/K/kg]
    κ::NF                   # = R_dry/cₚ, gas const for air over heat capacity

    # CORIOLIS FREQUENCY (scaled by radius as is vorticity) = 2Ω*sin(lat)*radius
    f_coriolis::Vector{NF}

    # ADIABATIC TERM
    σ_lnp_A::Vector{NF}             # σ-related factors needed for adiabatic terms 1st
    σ_lnp_B::Vector{NF}             # and 2nd

    # GEOPOTENTIAL INTEGRATION (on half/full levels)
    Δp_geopot_half::Vector{NF}      # = R*(ln(p_k+1) - ln(p_k+1/2)), for half level geopotential
    Δp_geopot_full::Vector{NF}      # = R*(ln(p_k+1/2) - ln(p_k)), for full level geopotential

    # REFERENCE TEMPERATURE PROFILE for implicit
    temp_ref_profile::Vector{NF}    # reference temperature profile
end

"""
$(TYPEDSIGNATURES)
Generator function for a DynamicsConstants struct.
"""
function DynamicsConstants( spectral_grid::SpectralGrid,
                            planet::AbstractPlanet,
                            atmosphere::AbstractAtmosphere,
                            geometry::Geometry)

    # PHYSICAL CONSTANTS
    (;R_dry, R_vapour, lapse_rate, cₚ) = atmosphere
    (;ΔT_stratosphere, σ_tropopause, temp_ref) = atmosphere
    (;NF, radius) = spectral_grid
    (;rotation, gravity) = planet
    layer_thickness = atmosphere.layer_thickness*1000   # ShallowWater: convert from [km] to [m]
    ξ = R_dry/R_vapour              # Ratio of gas constants: dry air / water vapour [1]
    μ_virt_temp = (1-ξ)/ξ           # used in Tv = T(1+μq), for conversion from humidity q
                                    # and temperature T to virtual temperature Tv
    κ = R_dry/cₚ                    # = 2/7ish for diatomic gas

    # CORIOLIS FREQUENCY (scaled by radius as is vorticity)
    (;sinlat) = geometry
    f_coriolis = 2rotation*sinlat*radius

    # ADIABATIC TERM, Simmons and Burridge, 1981, eq. 3.12
    (;σ_levels_half,σ_levels_full,σ_levels_thick) = geometry
    # precompute ln(σ_k+1/2) - ln(σ_k-1/2) but swap sign, include 1/Δσₖ
    σ_lnp_A = log.(σ_levels_half[1:end-1]./σ_levels_half[2:end]) ./ σ_levels_thick
    σ_lnp_A[1] = 0  # the corresponding sum is 1:k-1 so 0 to replace log(0) from above
    
    # precompute the αₖ = 1 - p_k-1/2/Δpₖ*log(p_k+1/2/p_k-1/2) term in σ coordinates
    σ_lnp_B = 1 .- σ_levels_half[1:end-1]./σ_levels_thick .*
                    log.(σ_levels_half[2:end]./σ_levels_half[1:end-1])
    σ_lnp_B[1] = σ_levels_half[1] <= 0 ? log(2) : σ_lnp_B[1]    # set α₁ = log(2), eq. 3.19
    σ_lnp_B .*= -1  # absorb sign from -1/Δσₖ only, eq. 3.12

    # GEOPOTENTIAL
    Δp_geopot_half, Δp_geopot_full = initialize_geopotential(σ_levels_full,σ_levels_half,R_dry)

    # VERTICAL REFERENCE TEMPERATURE PROFILE (TODO: don't initialize here but in initalize! ?)
    # integrate hydrostatic equation from pₛ to p, use ideal gas law p = ρRT and linear
    # temperature decrease with height: T = Tₛ - ΔzΓ with lapse rate Γ
    # for stratosphere (σ < σ_tropopause) increase temperature (Jablonowski & Williamson. 2006, eq. 5)
    RΓg⁻¹ = R_dry*lapse_rate/(1000*gravity)     # convert lapse rate from [K/km] to [K/m]
    temp_ref_profile = [temp_ref*σ^RΓg⁻¹ for σ in σ_levels_full]
    temp_ref_profile .+= [σ < σ_tropopause ? ΔT_stratosphere*(σ_tropopause-σ)^5 : 0 for σ in σ_levels_full]

    # This implies conversion to NF
    return DynamicsConstants{NF}(;  radius,rotation,gravity,layer_thickness,
                                    R_dry,R_vapour,μ_virt_temp,cₚ,κ,
                                    f_coriolis,
                                    σ_lnp_A,σ_lnp_B,
                                    Δp_geopot_half, Δp_geopot_full,
                                    temp_ref_profile)
end
