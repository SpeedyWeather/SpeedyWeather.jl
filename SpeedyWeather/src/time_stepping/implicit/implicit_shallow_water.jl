# dispatch over model.implicit and time stepping
implicit_correction!(vars::Variables, model::AbstractModel) =
    implicit_correction!(vars, model.implicit, model.time_stepping, model)

# model.implicit=nothing (for BarotropicModel)
initialize!(::Nothing, ::Variables, ::AbstractModel) = nothing
implicit_correction!(::Variables, ::Nothing, ::AbstractModel) = nothing

# SHALLOW WATER MODEL
export ImplicitShallowWater

"""Struct that holds various precomputed arrays for the semi-implicit correction to
prevent gravity waves from amplifying in the shallow water model. The implicit time step
between i-1 and i+1 is controlled by the centering parameter `α` in the range [0, 1]

    α = 0   means the gravity wave terms are evaluated at i-1 (forward)
    α = 0.5 evaluates at i+1 and i-1 (centered implicit)
    α = 1   evaluates at i+1 (backward implicit)
    α ∈ [0.5, 1] are also possible which controls the strength of the gravity wave dampening.
    α = 0.5 slows gravity waves and prevents them from amplifying
    α > 0.5 will dampen the gravity waves within days to a few timesteps (α=1)

Fields are
$(TYPEDFIELDS)"""
@kwdef struct ImplicitShallowWater{NF} <: AbstractImplicit
    "[OPTION] Centering coefficient for semi-implicit, 0.5 = Crank-Nicolson, 1 = Backwards Euler"
    centering::NF = 0.5
end

"""
$(TYPEDSIGNATURES)
Generator using the resolution from `spectral_grid`."""
ImplicitShallowWater(SG::SpectralGrid; kwargs...) = ImplicitShallowWater{SG.NF}(; kwargs...)
initialize!(::ImplicitShallowWater, ::Variables, ::AbstractModel) = nothing

"""
$(TYPEDSIGNATURES)
Apply correction to the tendencies in `diagn` to prevent the gravity waves from amplifying.
The correction is implicitly evaluated using the parameter `implicit.α` to switch between
forward, centered implicit or backward evaluation of the gravity wave terms."""
function implicit_correction!(
        vars::Variables,
        implicit::ImplicitShallowWater,
        time_stepping::AbstractLeapfrog,
        model::ShallowWater,
    )

    div_tend = get_tendency_step(vars.tendencies.divergence, time_stepping, implicit)
    η_tend = get_tendency_step(vars.tendencies.η, time_stepping, implicit)
    div_old, div_new = get_steps(vars.prognostic.divergence)    # divergence at t, t+dt
    η_old, η_new = get_steps(vars.prognostic.η)                 # η at t, t+dt

    H = model.atmosphere.layer_thickness        # layer thickness [m], undisturbed, no mountains
    g = model.planet.gravity                    # gravitational acceleration [m/s²]
    
    # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog)
    Δt = time_step(time_stepping, vars.prognostic.clock)       
    ξ = implicit.centering * Δt
    
    # Get precomputed l_indices from the spectrum
    l_indices = div_tend.spectrum.l_indices

    # GPU kernel launch
    arch = architecture(div_tend)
    launch!(
        arch, SpectralWorkOrder, size(div_tend), implicit_leapfrog_shallow_water_kernel!,
        div_tend, η_tend, div_old, div_new, η_old, η_new, l_indices, H, g, ξ
    )
    return nothing
end

@kernel inbounds = true function implicit_leapfrog_shallow_water_kernel!(
        div_tend, η_tend, div_old, div_new, η_old, η_new, l_indices,
        H, g, ξ,
    )
    I = @index(Global, Cartesian)
    lm = I[1]  # single index lm corresponding to harmonic l, m
    k = I[2]   # layer index

    # Use precomputed l index from spectrum
    l = l_indices[lm]

    # eigenvalue, with without 1/radius², 1-based -l*(l+1) → -l*(l-1)
    ∇² = -l * (l - 1)

    # Calculate the G = N(Vⁱ) + NI(Vⁱ⁻¹ - Vⁱ) term.
    # Vⁱ is a prognostic variable at time step i
    # N is the right hand side of ∂V\∂t = N(V)
    # NI is the part of N that's calculated semi-implicitily: N = NE + NI
    G_div = div_tend[lm, k] - g * ∇² * (η_old[lm] - η_new[lm])
    G_η = η_tend[lm] - H * (div_old[lm, k] - div_new[lm, k])

    # Using the Gs correct the tendencies for semi-implicit time stepping
    S⁻¹ = inv(1 - ξ^2 * H * g * ∇²)  # operator to invert
    div_tend[lm, k] = S⁻¹ * (G_div - ξ * g * ∇² * G_η)
    η_tend[lm] = G_η - ξ * H * div_tend[lm, k]
end

"""
$(TYPEDSIGNATURES)
Apply correction to the tendencies in `diagn` to prevent the gravity waves from amplifying.
The correction is implicitly evaluated using the parameter `implicit.α` to switch between
forward, centered implicit or backward evaluation of the gravity wave terms."""
function implicit_correction!(
        vars::Variables,
        implicit::ImplicitShallowWater,
        time_stepping::AbstractNCycleLorenz,
        model::ShallowWater
    )
    # implicit timestep ξ = α*Δt, depending on centering (0.5 Crank Nicolson, 1 backward Euler)
    (; Δt) = time_stepping    
    w = weight_coefficient(time_stepping, vars.prognostic.clock)
    ξ = w * implicit.centering * Δt

    # use Hotta et al. 2016 notation for current tendency F and previous (averaged) tendency G
    F_div = get_step(vars.tendencies.divergence, 1)    # no tuple (get_steps) here: a tuple of step
    G_div = get_step(vars.tendencies.divergence, 2)    # views breaks Enzyme on Julia >= 1.11
    F_η = get_step(vars.tendencies.η, 1)
    G_η = get_step(vars.tendencies.η, 2)

    H = model.atmosphere.layer_thickness        # layer thickness [m], undisturbed, no mountains
    g = model.planet.gravity                    # gravitational acceleration [m/s²]
        
    # Get precomputed l_indices from the spectrum
    (; l_indices) = F_div.spectrum

    # GPU kernel launch
    arch = architecture(F_div)
    launch!(
        arch, SpectralWorkOrder, size(F_div), implicit_ncycle_lorenz_shallow_water_kernel!,
        F_div, G_div, F_η, G_η, l_indices, H, g, ξ, w
    )
    return nothing
end

@kernel inbounds = true function implicit_ncycle_lorenz_shallow_water_kernel!(
        F_div, G_div, F_η, G_η, l_indices, H, g, ξ, w
    )
    I = @index(Global, Cartesian)
    lm = I[1]  # single index lm corresponding to harmonic l, m
    k = I[2]   # layer index

    # Use precomputed l index from spectrum
    l = l_indices[lm]

    # eigenvalue, with without 1/radius², 1-based -l*(l+1) → -l*(l-1)
    ∇² = -l * (l - 1)

    # N-Cycle Lorenz tendency average to form the RHS of the solve
    RHS_div = w*F_div[lm, k] + (1 - w)*G_div[lm, k]
    RHS_η = w*F_η[lm] + (1 - w)*G_η[lm]

    # Implicit solve
    S⁻¹ = inv(1 - ξ^2 * H * g * ∇²)
    δdiv = S⁻¹ * (RHS_div - ξ * g * ∇² * RHS_η)
    δη = RHS_η - ξ * H * δdiv

    # Reobtain F_explicit + L_implicit by undoing the Lorenz tendency average
    # and store in F tendencies so it can be plugged into Lorenz time stepping as normal
    F_div[lm, k] = (δdiv - (1 - w)*G_div[lm, k]) / w
    F_η[lm, k] = (δη - (1 - w)*G_η[lm]) / w
end