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

    "[OPTION] Surface layer fraction"
    fb::NF = 0.1

    "[OPTION] diffuse static energy?"
    diffuse_static_energy::Bool = true

    "[OPTION] diffuse momentum?"
    diffuse_momentum::Bool = true

    "[OPTION] diffuse humidity? Ignored for PrimitiveDryModels"
    diffuse_humidity::Bool = true

    # precomputed constants
    κZsqrtC_g::Base.RefValue{NF} = Ref(zero(NF))

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
        k₋ = max(k, 1)      # sets the gradient across surface and top to 0
        k₊ = min(k, nlev)   # = no flux boundary conditions
        scheme.∇²_above[k] = inv(2*(σ[k] - σ[k₋]) * (σ_half[k+1] - σ_half[k]))
        scheme.∇²_below[k] = inv(2*(σ[k₊] - σ[k]) * (σ_half[k+1] - σ_half[k]))
    end

    # Typical height Z of lowermost layer from geopotential of reference surface temperature
    # minus surface geopotential (orography * gravity)
    (; temp_ref) = model.atmosphere
    (; gravity) = model.planet
    (; Δp_geopot_full) = model.geopotential
    Z = temp_ref * Δp_geopot_full[end] / gravity
    
    # maximum drag Cmax from that height, stable conditions would decrease Cmax towards 0
    # Frierson 2006, eq (12)
    (; κ, z₀) = scheme
    scheme.κZsqrtC_g[] = Z*κ^2/(log(Z/z₀)*gravity)

    return nothing
end

function vertical_diffusion!(
    column::ColumnVariables,
    scheme::BulkRichardsonDiffusion,
    model::PrimitiveEquation)

    (; diffuse_momentum, diffuse_static_energy, diffuse_humidity) = scheme

    # escape immediately if all diffusions disabled
    any((diffuse_momentum, diffuse_static_energy, diffuse_humidity)) || return nothing
    
    K, h = get_diffusion_coefficients!(column, scheme, model.atmosphere)

    (; u, u_tend, v, v_tend, dry_static_energy, temp_tend, humid, humid_tend) = column
    diffuse_momentum                            && vertical_diffusion!(u_tend, u, K, h, scheme)
    diffuse_momentum                            && vertical_diffusion!(v_tend, v, K, h, scheme)
    diffuse_static_energy                       && vertical_diffusion!(temp_tend, dry_static_energy, K, h, scheme)
    model isa PrimitiveWet && diffuse_humidity  && vertical_diffusion!(humid_tend, humid, K, h, scheme)
    return nothing
end

function get_diffusion_coefficients!(
    column::ColumnVariables,
    scheme::BulkRichardsonDiffusion,
    atmosphere::AbstractAtmosphere,
)
    Ri = bulk_richardson!(column, atmosphere)

    K = column.b        # reuse work array for diffusion coefficients
    (; Ri_c, fb, κsqrtC_g) = scheme
    (; nlev, geopot) = column
    Ri_a = Ri[nlev]     # surface bulk Richardson number

    # Boundary layer depth is highest layer for which Ri < Ri_c (the "critical" threshold)
    # as well as all layers below
    h::Int = nlev
    while Ri[h] < Ri_c && h > 0
        h -= 1
    end
    h += 1  # uppermost layer where Ri < Ri_c
    column.boundary_layer_depth = h

    # Calculate diffusion coefficients
    surface_speed = sqrt(u[nlev]^2 + v[nlev]^2)
    Kb = κZsqrtC_g * surface_speed

    for k in 1:nlev
        K[k] = k > h ? 0 : Kb
    end

    return K, h
end

function vertical_diffusion!(   
    tend::AbstractVector,       # tendency to accumulate diffusion into
    var::AbstractVector,        # variable calculate diffusion from
    K::AbstractVector,          # diffusion coefficients
    h::Int,                     # uppermost layer that's still within the boundary layer
    scheme::BulkRichardsonDiffusion,
)
    (; ∇²_above, ∇²_below) = scheme
    nlev = length(tend)
    nlev == 1 && return nothing     # escape immediately for single-layer (no diffusion)

    @boundscheck nlev == length(var) == length(K) || throw(BoundsError)

    for k in nlev:-1:h          # diffusion only in surface boundary layer of thickness h
        
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