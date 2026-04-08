abstract type AbstractImplicit <: AbstractModelComponent end

# Barotropic model has no implicit correction
initialize!(::Nothing, ::Real, ::DiagnosticVariables, ::AbstractModel) = nothing
implicit_correction!(::DiagnosticVariables, ::PrognosticVariables, ::Nothing, ::AbstractModel, ::Real) = nothing
set_initialized!(::Nothing) = nothing

# =============================================================================
# SHALLOW WATER MODEL
# =============================================================================

export ImplicitShallowWater

"""Precomputed arrays for the semi-implicit correction in the shallow water model.
Gravity waves are evaluated implicitly with weight α ∈ [0.5, 1]:
  α = 0.5  centred implicit (Crank–Nicolson), second-order, no damping
  α = 1    fully backward implicit, strongly damps gravity waves

The effective implicit time step is ξ = α·Δt.
$(TYPEDFIELDS)"""
@kwdef mutable struct ImplicitShallowWater{NF} <: AbstractImplicit
    "[OPTION] Semi-implicit coefficient, 0.5 ≤ α ≤ 1"
    α::NF = 1

    "Implicit time step ξ = α·Δt [s/m, radius-scaled]"
    time_step::NF = 0
end

ImplicitShallowWater(SG::SpectralGrid; kwargs...) = ImplicitShallowWater{SG.NF}(; kwargs...)

"""$(TYPEDSIGNATURES)
Set the implicit time step ξ = α·dt. For the Lorenz N-cycle, dt = Δt (not 2Δt)."""
function initialize!(implicit::ImplicitShallowWater, dt::Real, args...)
    implicit.time_step = implicit.α * dt
end

set_initialized!(::ImplicitShallowWater) = nothing

"""$(TYPEDSIGNATURES)
Semi-implicit correction for the shallow water model (Lorenz N-cycle version).

For each spectral mode (l, m):

1. Extract the purely explicit tendencies by removing the implicit terms at x:
     F^E_div  = div_tend  - L^I_div(pres_x)   = div_tend  + g∇²·pres_x
     F^E_pres = pres_tend - L^I_pres(div_x)   = pres_tend + H·div_x

2. Blend with the stored G (step 2 of each prognostic variable):
     G_div  ← w·F^E_div  + (1-w)·G_div_stored
     G_pres ← w·F^E_pres + (1-w)·G_pres_stored
   Store new G back into step 2 for the next time step.

3. Form the implicit RHS by re-adding the implicit terms at x:
     rhs_D = G_div  - g∇²·pres_x   (= G_div  + g·l·(l-1)·pres_x)
     rhs_η = G_pres - H·div_x

4. Solve the coupled implicit system with ξ = α·Δt:
     S⁻¹ = inv(1 - ξ²·H·g·∇²)
     δD  = S⁻¹·(rhs_D - ξ·g·∇²·rhs_η)
     δη  = rhs_η - ξ·H·δD

The output tendencies δD and δη are written back to `div_tend` and `pres_tend`;
`lorenz_step!` then updates the state as x ← x + Δt·δ."""
function implicit_correction!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    implicit::ImplicitShallowWater,
    model::ShallowWater,
    w::Real,
)
    (; div_tend, pres_tend) = diagn.tendencies
    div_x, G_div   = get_steps(progn.div)    # step 1 = current state, step 2 = stored G
    pres_x, G_pres = get_steps(progn.pres)

    H  = model.atmosphere.layer_thickness
    g  = model.planet.gravity
    ξ  = implicit.time_step

    l_indices = div_tend.spectrum.l_indices
    arch = architecture(div_tend)

    launch!(arch, SpectralWorkOrder, size(div_tend),
            implicit_shallow_water_lorenz_kernel!,
            div_tend, pres_tend, div_x, G_div, pres_x, G_pres,
            l_indices, H, g, ξ, w)

    zero_last_degree!(div_tend)
    zero_last_degree!(pres_tend)
end

@kernel inbounds=true function implicit_shallow_water_lorenz_kernel!(
    div_tend, pres_tend, div_x, G_div, pres_x, G_pres, @Const(l_indices),
    @Const(H), @Const(g), @Const(ξ), @Const(w)
)
    I = @index(Global, Cartesian)
    lm = I[1]
    k  = I[2]

    l  = l_indices[lm]
    ∇² = -l*(l-1)     # eigenvalue of ∇², always ≤ 0; 1-based l so -l*(l-1) not -l*(l+1)

    # --- Step 1: extract explicit tendencies (remove implicit at current x) ---
    # div_tend = F^E_div + L^I_div(pres_x),  L^I_div(η) = -g·∇²·η = g·l·(l-1)·η
    # → F^E_div = div_tend - (-g·∇²·pres_x) = div_tend + g·∇²·pres_x
    f_e_div  = div_tend[lm, k] + g * ∇² * pres_x[lm]

    # pres_tend = F^E_pres + L^I_pres(div_x),  L^I_pres(D) = -H·D
    # → F^E_pres = pres_tend + H·div_x
    f_e_pres = pres_tend[lm] + H * div_x[lm, k]

    # --- Step 2: blend with stored G and update storage ---
    g_div  = w * f_e_div  + (1 - w) * G_div[lm, k]
    g_pres = w * f_e_pres + (1 - w) * G_pres[lm]
    G_div[lm, k] = g_div
    G_pres[lm]   = g_pres    # safe: SWE has nlayers=1, so k=1 always

    # --- Step 3: re-add implicit at x to form the solve RHS ---
    rhs_D = g_div  - g * ∇² * pres_x[lm]    # = G_div + g·l·(l-1)·pres_x
    rhs_η = g_pres - H  * div_x[lm, k]

    # --- Step 4: implicit solve ---
    S⁻¹ = inv(1 - ξ^2 * H * g * ∇²)
    div_tend[lm, k] = S⁻¹ * (rhs_D - ξ * g * ∇² * rhs_η)
    pres_tend[lm]   = rhs_η - ξ * H * div_tend[lm, k]
end
