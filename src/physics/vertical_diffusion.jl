abstract type AbstractVerticalDiffusion <: AbstractParameterization end

## FUNCTION BARRIERS for all AbstractVerticalDiffusion
function vertical_diffusion!(   column::ColumnVariables,
                                model::PrimitiveEquation)
    vertical_diffusion!(column, model.vertical_diffusion, model)
end

## NO VERTICAL DIFFUSION
export NoVerticalDiffusion
struct NoVerticalDiffusion <: AbstractVerticalDiffusion end
NoVerticalDiffusion(SG::SpectralGrid) = NoVerticalDiffusion()

# define dummy functions
initialize!(::NoVerticalDiffusion,::PrimitiveEquation) = nothing
vertical_diffusion!(::ColumnVariables, ::NoVerticalDiffusion, ::PrimitiveEquation) = nothing


Base.@kwdef struct BulkRichardsonDiffusion{NF} <: AbstractVerticalDiffusion
    nlev::Int

    "[OPTION] von Kármán constant [1]"
    κ::NF = 0.4

    "[OPTION] roughness length [m]"
    z₀::NF = 3.21e-5

    "[OPTION] Critical Richardson number for stable mixing cutoff [1]"
    Ri_c::NF = 1

    "[OPTION] Fraction of surface boundary layer"
    fb::NF = 0.1

    "[OPTION] diffuse static energy?"
    diffuse_static_energy::Bool = true

    "[OPTION] diffuse momentum?"
    diffuse_momentum::Bool = true

    "[OPTION] diffuse humidity? Ignored for PrimitiveDryModels"
    diffuse_humidity::Bool = true

    # precomputed constants
    sqrtC_max::Base.RefValue{NF} = Ref(zero(NF))

    # precomputed operators
    ∇²_above::Vector{NF} = zeros(NF, nlev)
    ∇²_below::Vector{NF} = zeros(NF, nlev)
end

BulkRichardsonDiffusion(SG::SpectralGrid; kwargs...) = BulkRichardsonDiffusion{SG.NF}(; nlev=SG.nlev, kwargs...)

function initialize!(scheme::BulkRichardsonDiffusion, model::PrimitiveEquation)
    
    (; nlev) = model.geometry
    nlev == 1 && return nothing     # no diffusion for 1-layer model

    # ∇² operator on σ levels like 1/Δσ² but for variable Δσ
    # also includes a 1/2 so that the diffusion coefficients on full levels can be added
    # which is equivalent to interpolating them on half levels for a ∂σ (K ∂σ) formulation
    # with σ-dependent diffusion coefficient K
    σ = model.geometry.σ_levels_full
    σ_half = model.geometry.σ_levels_half

    for k in 1:nlev
        σ₋ = k <= 1    ? -Inf : σ[k-1]   # sets the gradient across surface and top to 0
        σ₊ = k >= nlev ?  Inf : σ[k+1]   # = no flux boundary conditions
        scheme.∇²_above[k] = inv(2*(σ[k] - σ₋) * (σ_half[k+1] - σ_half[k]))
        scheme.∇²_below[k] = inv(2*(σ₊ - σ[k]) * (σ_half[k+1] - σ_half[k]))
    end

    # Typical height Z of lowermost layer from geopotential of reference surface temperature
    # minus surface geopotential (orography * gravity), simplification compared to
    # Frierson to reduce the number of expensive log calls given that z ≈ Z for most
    # surface temperature variations
    (; temp_ref) = model.atmosphere
    (; gravity) = model.planet
    (; Δp_geopot_full) = model.geopotential
    Z = temp_ref * Δp_geopot_full[end] / gravity
    
    # maximum drag Cmax from that height, stable conditions would decrease Cmax towards 0
    # Frierson 2006, eq (12)
    (; κ, z₀) = scheme
    scheme.sqrtC_max[] = κ/log(Z/z₀)

    return nothing
end

function vertical_diffusion!(
    column::ColumnVariables,
    scheme::BulkRichardsonDiffusion,
    model::PrimitiveEquation)

    (; diffuse_momentum, diffuse_static_energy, diffuse_humidity) = scheme

    # escape immediately if all diffusions disabled
    any((diffuse_momentum, diffuse_static_energy, diffuse_humidity)) || return nothing
    
    K, kₕ = get_diffusion_coefficients!(column, scheme, model.atmosphere, model.planet)

    (; u, u_tend, v, v_tend, dry_static_energy, temp_tend, humid, humid_tend) = column
    diffuse_momentum                            && vertical_diffusion!(u_tend, u, K, kₕ, scheme)
    diffuse_momentum                            && vertical_diffusion!(v_tend, v, K, kₕ, scheme)
    model isa PrimitiveWet && diffuse_humidity  && vertical_diffusion!(humid_tend, humid, K, kₕ, scheme)
    K ./= model.atmosphere.heat_capacity        # put dry static energy => temperature conversion into K
    diffuse_static_energy                       && vertical_diffusion!(temp_tend, dry_static_energy, K, kₕ, scheme)
    return nothing
end

function get_diffusion_coefficients!(
    column::ColumnVariables,
    scheme::BulkRichardsonDiffusion,
    atmosphere::AbstractAtmosphere,
    planet::AbstractPlanet,
)
    K = column.b        # reuse work array for diffusion coefficients
    (; Ri_c, κ, z₀, fb) = scheme
    (; nlev, u, v, geopot, orography) = column
    gravity⁻¹ = inv(planet.gravity)

    # Boundary layer depth is highest layer for which Ri < Ri_c (the "critical" threshold)
    # as well as all layers below
    Ri = bulk_richardson!(column, atmosphere)
    kₕ::Int = nlev
    while kₕ > 0 && Ri[kₕ] < Ri_c
        kₕ -= 1
    end
    kₕ += 1  # uppermost layer where Ri < Ri_c
    column.boundary_layer_depth = kₕ

    if kₕ <= nlev   # boundary layer depth is at least 1 layer thick (calculate diffusion)

        # Calculate diffusion coefficients following Frierson 2006, eq. 16-20
        h = max(geopot[kₕ]*gravity⁻¹ - orography, 0)# always positive to avoid error in log 
        Ri_a = Ri[nlev]                             # surface bulk Richardson number
        Ri_a = clamp(Ri_a, 0, Ri_c)                 # cases of eq. 12-14
        sqrtC = scheme.sqrtC_max[]*(1-Ri_a/Ri_c)    # sqrt of eq. 12-14
        surface_speed = sqrt(u[nlev]^2 + v[nlev]^2)
        K0 = κ * surface_speed * sqrtC              # height-independent K eq. 19, 20

        K[1:kₕ-1] .= 0                              # diffusion above boundary layer 0
        for k in kₕ:nlev
            z = max(geopot[k]*gravity⁻¹ - orography, z₀)    # height [m] above surface
            zmin = min(z, fb*h)         # height [m] to evaluate Kb(z) at
            K_k = K0 * zmin             # = κuₐ√Cz in eq. (19, 20)

            # multiply with z-dependent factor in eq. (18) ?
            K_k *= z < fb*h ? 1 : zfac(z, h, fb)

            # multiply with Ri-dependent factor in eq. (20) ?
            K_k *= Ri[kₕ] <= 0 ? 1 : Rifac(Ri[kₕ], Ri_c, zmin, z₀)
            K[k] = K_k                  # write diffusion coefficient into array
        end
    else
        fill!(K,0)
    end

    # return diffusion coefficients and height index of boundary layer
    return K, kₕ
end

# z-dependent factor in Frierson, 2006 eq (18)
@inline zfac(z,h,fb) = z/(fb*h)*(1 - (z - fb*h) / ((1-fb) * h))^2

# Ri-dependent factor in Frierson, 2006 eq (20)
@inline function Rifac(Ri, Ri_c, z, z₀)
    Ri_Ri_c = Ri/Ri_c
    inv(1 + Ri_Ri_c * log(z/z₀) / (1 - Ri_Ri_c))
end

function vertical_diffusion!(   
    tend::AbstractVector,       # tendency to accumulate diffusion into
    var::AbstractVector,        # variable calculate diffusion from
    K::AbstractVector,          # diffusion coefficients
    kₕ::Int,                    # uppermost layer that's still within the boundary layer
    scheme::BulkRichardsonDiffusion,
)
    (; ∇²_above, ∇²_below) = scheme
    nlev = length(tend)
    nlev == 1 && return nothing     # escape immediately for single-layer (no diffusion)

    @boundscheck nlev == length(var) == length(K) || throw(BoundsError)

    for k in nlev:-1:kₕ          # diffusion only in surface boundary layer of thickness h
        
        # sets the gradient across surface and top to 0 = no flux boundary conditions
        k₋ = max(k, 1)      # index above (- in σ direction which is 0 at top and 1 at surface)
        k₊ = min(k, nlev)   # index below (+ in σ direction which is 0 at top and 1 at surface)

        K_∂var_below = (var[k₊] - var[k]) * (K[k₊] + K[k])  # average diffusion coefficient K here
        K_∂var_above = (var[k] - var[k₋]) * (K[k] + K[k₋])  # but 1/2 is already baked into the ∇² operators

        tend[k] += ∇²_below[k]*K_∂var_below - ∇²_above[k]*K_∂var_above
    end
end

"""
$(TYPEDSIGNATURES)
Calculate the bulk richardson number following Frierson, 2007.
For vertical stability in the boundary layer."""
function bulk_richardson!(
    column::ColumnVariables,
    atmosphere::AbstractAtmosphere,
)
    cₚ = atmosphere.heat_capacity
    (; u, v, geopot, temp_virt, nlev) = column
    bulk_richardson = column.a      # reuse work array

    # surface layer
    V² = u[nlev]^2 + v[nlev]^2
    Θ₀ = cₚ*temp_virt[nlev]
    Θ₁ = Θ₀ + geopot[nlev]
    bulk_richardson[nlev] = geopot[nlev]*(Θ₁ - Θ₀) / (Θ₀*V²)

    @inbounds for k in 1:nlev-1
        V² = u[k]^2 + v[k]^2
        virtual_dry_static_energy = cₚ * temp_virt[k] + geopot[k]
        bulk_richardson[k] = geopot[k]*(virtual_dry_static_energy - Θ₁) / (Θ₁*V²)
    end

    return bulk_richardson
end