abstract type AbstractVerticalDiffusion <: AbstractParameterization end

export BulkRichardsonDiffusion
@kwdef struct BulkRichardsonDiffusion{NF, RefNF <: Ref{NF}, VectorType} <: AbstractVerticalDiffusion
    "[OPTION] von Kármán constant [1]"
    κ::NF = 0.4

    "[OPTION] roughness length [m]"
    z₀::NF = 3.21e-5

    "[OPTION] Critical Richardson number for stable mixing cutoff [1]"
    Ri_c::NF = 10

    "[OPTION] Fraction of surface boundary layer"
    fb::NF = 0.1

    "[OPTION] diffuse static energy?"
    diffuse_static_energy::Bool = true

    "[OPTION] diffuse momentum?"
    diffuse_momentum::Bool = true

    "[OPTION] diffuse humidity? Ignored for PrimitiveDryModels"
    diffuse_humidity::Bool = true

    "[DERIVED] Logarithm of typical lowermost layer height Z over roughness length z₀"
    logZ_z₀::RefNF = Ref(zero(NF))

    "[DERIVED] Sqrt of maximum diffusion coefficient C"
    sqrtC_max::RefNF = Ref(zero(NF))

    "[DERIVED] Vertical Laplace operator, operator for cells above"
    ∇²_above::VectorType

    "[DERIVED] Vertical Laplace operator, operator for cells below"
    ∇²_below::VectorType
end

Adapt.@adapt_structure BulkRichardsonDiffusion

# generator function
function BulkRichardsonDiffusion(SG::SpectralGrid; kwargs...)
    arch = SG.architecture
    ∇²_above = on_architecture(arch, zeros(SG.NF, SG.nlayers))
    ∇²_below = on_architecture(arch, zeros(SG.NF, SG.nlayers))
    return BulkRichardsonDiffusion{SG.NF, Ref{SG.NF}, SG.VectorType}(; ∇²_above, ∇²_below, kwargs...)
end

variables(::BulkRichardsonDiffusion) = (
    # TODO change to height in meters or index?
    DiagnosticVariable(name=:boundary_layer_height, dims=Grid2D(), desc="Boundary layer height", units="1"),
)

# function barrier to unpack only model components needed
function initialize!(scheme::BulkRichardsonDiffusion, model::PrimitiveEquation)
    model_parameters = get_model_parameters(model)
    initialize!(scheme, model_parameters)
end

function initialize!(scheme::BulkRichardsonDiffusion, model)

    (; nlayers) = model.geometry
    nlayers == 1 && return nothing     # no diffusion for 1-layer model

    # ∇² operator on σ levels like 1/Δσ² but for variable Δσ
    # also includes a 1/2 so that the diffusion coefficients on full levels can be added
    # which is equivalent to interpolating them on half levels for a ∂σ (K ∂σ) formulation
    # with σ-dependent diffusion coefficient K
    σ = on_architecture(CPU(), model.geometry.σ_levels_full)
    σ_half = on_architecture(CPU(), model.geometry.σ_levels_half)
    ∇²_above = on_architecture(CPU(), scheme.∇²_above)
    ∇²_below = on_architecture(CPU(), scheme.∇²_below)

    for k in 1:nlayers
        σ₋ = k <= 1       ? -Inf : σ[k-1]   # sets the gradient across surface and top to 0
        σ₊ = k >= nlayers ?  Inf : σ[k+1]   # = no flux boundary conditions
        ∇²_above[k] = inv(2*(σ[k] - σ₋) * (σ_half[k+1] - σ_half[k]))
        ∇²_below[k] = inv(2*(σ₊ - σ[k]) * (σ_half[k+1] - σ_half[k]))
    end

    arch = model.architecture
    scheme.∇²_above .= on_architecture(arch, ∇²_above)
    scheme.∇²_below .= on_architecture(arch, ∇²_below)

    # Typical height Z of lowermost layer from geopotential of reference surface temperature
    # minus surface geopotential (orography * gravity), simplification compared to
    # Frierson to reduce the number of expensive log calls given that z ≈ Z for most
    # surface temperature variations
    (; temp_ref) = model.atmosphere
    (; gravity) = model.planet
    Δp_geopot_full = on_architecture(CPU(), model.geopotential.Δp_geopot_full)
    Z = GPUArrays.@allowscalar temp_ref * Δp_geopot_full[end] / gravity
    
    # maximum drag Cmax from that height, stable conditions would decrease Cmax towards 0
    # Frierson 2006, eq (12)
    (; κ, z₀) = scheme
    scheme.logZ_z₀[] = log(Z/z₀)        # ≈ log(z/z₀) avoid log during integration
    scheme.sqrtC_max[] = κ/log(Z/z₀)

    return nothing
end

# function barrier
parameterization!(ij, diagn, progn, diffusion::BulkRichardsonDiffusion, model) =
    vertical_diffusion!(ij, diagn, diffusion, model.atmosphere, model.planet, model.orography, model.class)

function vertical_diffusion!(
    ij,
    diagn,
    diffusion::BulkRichardsonDiffusion,
    atmosphere,
    planet,
    orography,
    Model,
)
    (; diffuse_momentum, diffuse_static_energy, diffuse_humidity) = diffusion

    # escape immediately if all diffusions disabled
    any((diffuse_momentum, diffuse_static_energy, diffuse_humidity)) || return nothing
    
    K, kₕ = get_diffusion_coefficients!(ij, diagn, diffusion, atmosphere, planet, orography)

    u_tend = diagn.tendencies.u_tend_grid
    v_tend = diagn.tendencies.v_tend_grid
    temp_tend = diagn.tendencies.temp_tend_grid
    humid_tend = diagn.tendencies.humid_tend_grid

    # TODO previous time step?
    u = diagn.grid.u_grid
    v = diagn.grid.v_grid
    humid = diagn.grid.humid_grid

    diffuse_momentum                            && vertical_diffusion!(ij, u_tend, u, K, kₕ, diffusion)
    diffuse_momentum                            && vertical_diffusion!(ij, v_tend, v, K, kₕ, diffusion)
    Model <: PrimitiveWet && diffuse_humidity   && vertical_diffusion!(ij, humid_tend, humid, K, kₕ, diffusion)

    if diffuse_static_energy
        # compute dry static energy on the fly
        dry_static_energy = diagn.dynamics.a_grid
        cₚ = atmosphere.heat_capacity
        T = diagn.grid.temp_grid
        Φ = diagn.grid.geopotential

        for k in 1:size(T, 2)
            dry_static_energy[ij, k] = cₚ * T[ij, k] + Φ[ij, k]
            K[ij, k] /= cₚ        # put temperature => dry static energy conversion into K
        end

        vertical_diffusion!(ij, temp_tend, dry_static_energy, K, kₕ, diffusion)
    end
    return nothing
end

function get_diffusion_coefficients!(
    ij,
    diagn,
    scheme::BulkRichardsonDiffusion,
    atmosphere::AbstractAtmosphere,
    planet::AbstractPlanet,
    orography,
)
    # reuse scratch array for diffusion coefficients
    K = diagn.dynamics.b_grid
    nlayers = size(K, 2)
    (; Ri_c, κ, z₀, fb) = scheme
    logZ_z₀ = scheme.logZ_z₀[]
    gravity⁻¹ = inv(planet.gravity)

    u = diagn.grid.u_grid
    v = diagn.grid.v_grid
    geopotential = diagn.grid.geopotential

    # Boundary layer depth is highest layer for which Ri < Ri_c (the "critical" threshold)
    # as well as all layers below
    Ri = bulk_richardson!(ij, diagn, atmosphere)
    kₕ::Int = nlayers
    while kₕ > 0 && Ri[ij, kₕ] < Ri_c
        kₕ -= 1
    end
    kₕ += 1  # uppermost layer where Ri < Ri_c

    # for output, TODO as layer index or height?
    diagn.physics.boundary_layer_height[ij] = kₕ

    # diffusion above boundary layer is 0
    for k in 1:nlayers  # set for all layers 
        K[ij, k] = 0                
    end

    if kₕ <= nlayers    # boundary layer depth is at least 1 layer thick (calculate diffusion)

        # Calculate diffusion coefficients following Frierson 2006, eq. 16-20
        # h always non-negative to avoid error in log 
        h = max(geopotential[ij, kₕ]*gravity⁻¹ - orography[ij], 0)  
        Ri_N = Ri[ij, nlayers]                      # surface bulk Richardson number
        Ri_N = clamp(Ri_N, 0, Ri_c)                 # cases of eq. 12-14
        sqrtC = scheme.sqrtC_max[]*(1-Ri_N/Ri_c)    # sqrt of eq. 12-14
        surface_speed = sqrt(u[ij, nlayers]^2 + v[ij, nlayers]^2)
        K0 = κ * surface_speed * sqrtC              # height-independent K eq. 19, 20                           

        for k in kₕ:nlayers
            # height [m] above surface
            z = max(geopotential[ij, k]*gravity⁻¹ - orography[ij], z₀)    
            zmin = min(z, fb*h)         # height [m] to evaluate Kb(z) at
            K_k = K0 * zmin             # = κ*u_N*√Cz in eq. (19, 20)
            
            # multiply with z-dependent factor in eq. (18) ?
            K_k *= z < fb*h ? one(K0) : zfac(z, h, fb)

            # multiply with Ri-dependent factor in eq. (20) ?
            # TODO use Ri[kₕ] or Ri_N here? 
            K_k *= Ri[ij, kₕ] <= 0 ? one(K0) : Rifac(Ri[ij, kₕ], Ri_c, logZ_z₀)
            K[ij, k] = K_k              # write diffusion coefficient into array
        end
    end

    # return diffusion coefficients and height index of boundary layer
    return K, kₕ
end

# z-dependent factor in Frierson, 2006 eq (18)
@inline zfac(z, h, fb) = z/(fb*h)*(1 - (z - fb*h) / ((1-fb) * h))^2

# Ri-dependent factor in Frierson, 2006 eq (20)
@inline function Rifac(Ri, Ri_c, z, z₀)
    Ri_Ri_c = Ri/Ri_c
    return inv(1 + Ri_Ri_c * log(z/z₀) / (1 - Ri_Ri_c))
end

# Approximate: Ri-dependent factor in Frierson, 2006 eq (20)
# because 1 / (1 + log(z/z₀)) is so weakly dependent on z for 10-10000m
@inline function Rifac(Ri, Ri_c, logz_z₀)
    Ri_Ri_c = Ri/Ri_c
    return inv(1 + Ri_Ri_c * logz_z₀ / (1 - Ri_Ri_c))
end

function vertical_diffusion!(   
    ij,                     # horizontal grid point ij
    tend::AbstractField,    # tendency to accumulate diffusion into
    var::AbstractField,     # variable calculate diffusion from
    K::AbstractField,       # diffusion coefficients
    kₕ::Int,                # uppermost layer that's still within the boundary layer
    diffusion::BulkRichardsonDiffusion,
)
    (; ∇²_above, ∇²_below) = diffusion
    nlayers = size(tend, 2)

    for k in kₕ:nlayers         # diffusion only in surface boundary layer of thickness h
        
        # sets the gradient across surface and top to 0 = no flux boundary conditions
        k₋ = max(k, 1)          # index above (- in σ direction which is 0 at top and 1 at surface)
        k₊ = min(k, nlayers)    # index below (+ in σ direction which is 0 at top and 1 at surface)

        K_∂var_below = (var[ij, k₊] - var[ij, k]) * (K[ij, k₊] + K[ij, k])  # average diffusion coefficient K here
        K_∂var_above = (var[ij, k] - var[ij, k₋]) * (K[ij, k] + K[ij, k₋])  # but 1/2 is already baked into the ∇² operators

        tend[ij, k] += ∇²_below[k]*K_∂var_below - ∇²_above[k]*K_∂var_above
    end
end

"""
$(TYPEDSIGNATURES)
Calculate the bulk Richardson number following Frierson, 2007.
For vertical stability in the boundary layer."""
function bulk_richardson!(
    ij,
    diagn,
    atmosphere::AbstractAtmosphere,
)
    # reuse work array
    Ri = diagn.dynamics.a_grid
    nlayers = size(Ri, 2)
    surface = nlayers       # surface index
    cₚ = atmosphere.heat_capacity

    # TODO previous time step?
    u = diagn.grid.u_grid
    v = diagn.grid.v_grid
    Φ = diagn.grid.geopotential
    Tᵥ = diagn.grid.temp_virt_grid

    # surface layer
    V² = u[ij, surface]^2 + v[ij, surface]^2
    Θ₀ = cₚ*Tᵥ[ij, surface]
    Θ₁ = Θ₀ + Φ[ij, surface]
    Ri[ij, surface] = Φ[ij, surface]*(Θ₁ - Θ₀) / (Θ₀*V²)

    @inbounds for k in 1:nlayers-1
        V² = u[ij, k]^2 + v[ij, k]^2
        virtual_dry_static_energy = cₚ * Tᵥ[ij, k] + Φ[ij, k]
        Ri[ij, k] = Φ[ij, k]*(virtual_dry_static_energy - Θ₁) / (Θ₁*V²)
    end

    return Ri
end