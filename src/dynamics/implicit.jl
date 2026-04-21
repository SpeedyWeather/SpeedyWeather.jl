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

# =============================================================================
# PRIMITIVE EQUATION MODEL
# =============================================================================

export ImplicitPrimitiveEquation

"""Precomputed arrays for the semi-implicit correction in the primitive equation model.
$(TYPEDFIELDS)"""
@kwdef mutable struct ImplicitPrimitiveEquation{
    NF,
    VectorType,
    MatrixType,
    TensorType,
} <: AbstractImplicit

    "Spectral resolution"
    trunc::Int

    "Number of vertical layers"
    nlayers::Int

    "Time-step coefficient: 0=explicit, 0.5=centred implicit, 1=backward implicit"
    α::NF = 1

    "Reinitialize at restart when initialized=true"
    reinitialize::Bool = true

    "Flag automatically set to true when initialize! has been called"
    initialized::Bool = false

    "vertical temperature profile, obtained from diagn on first time step"
    temp_profile::VectorType = zeros(NF, nlayers)

    "time step 2α·Δt packed in RefValue for mutability"
    ξ::Base.RefValue{NF} = Ref{NF}(0)

    "divergence: operator for the geopotential calculation"
    R::MatrixType = zeros(NF, nlayers, nlayers)

    "divergence: the -RdTₖ∇² term excl the eigenvalues from ∇² for divergence"
    U::VectorType = zeros(NF, nlayers)

    "temperature: operator for the TₖD + κTₖDlnps/Dt term"
    L::MatrixType = zeros(NF, nlayers, nlayers)

    "pressure: vertical averaging of the -D̄ term in the log surface pres equation"
    W::VectorType = zeros(NF, nlayers)

    "components to construct L, 1/2Δσ"
    L0::VectorType = zeros(NF, nlayers)

    "vert advection term in the temperature equation (below+above)"
    L1::MatrixType = zeros(NF, nlayers, nlayers)

    "factor in front of the div_sum_above term"
    L2::VectorType = zeros(NF, nlayers)

    "_sum_above operator itself"
    L3::MatrixType = zeros(NF, nlayers, nlayers)

    "factor in front of div term in Dlnps/Dt"
    L4::VectorType = zeros(NF, nlayers)

    "for every l the matrix to be inverted"
    S::MatrixType = zeros(NF, nlayers, nlayers)

    "combined inverted operator: S = 1 - ξ²(RL + UW)"
    S⁻¹::TensorType = zeros(NF, trunc+2, nlayers, nlayers)
end

function ImplicitPrimitiveEquation(spectral_grid::SpectralGrid; kwargs...)
    (; NF, VectorType, MatrixType, TensorType, trunc, nlayers) = spectral_grid
    return ImplicitPrimitiveEquation{NF, VectorType, MatrixType, TensorType}(;
        trunc, nlayers, kwargs...)
end

function initialize!(
    I::ImplicitPrimitiveEquation,
    dt::Real,
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)
    model.dynamics || return nothing
    (; geometry, geopotential, atmosphere, adiabatic_conversion) = model
    initialize!(I, dt, diagn, geometry, geopotential, atmosphere, adiabatic_conversion)
end

function initialize!(
    implicit::ImplicitPrimitiveEquation,
    dt::Real,
    diagn::DiagnosticVariables,
    geometry::AbstractGeometry,
    geopotential::AbstractGeopotential,
    atmosphere::AbstractAtmosphere,
    adiabatic_conversion::AbstractAdiabaticConversion,
)
    NF = eltype(diagn)
    (implicit.initialized && !implicit.reinitialize) && return nothing

    (; trunc, nlayers, α) = implicit
    (; σ_levels_full, σ_levels_thick) = geometry
    (; R_dry, κ) = atmosphere
    (; Δp_geopot_half, Δp_geopot_full) = geopotential
    (; σ_lnp_A, σ_lnp_B) = adiabatic_conversion

    arch = architecture(implicit.S)

    temp_profile_cpu, S_cpu, S⁻¹_cpu, L_cpu, R_cpu, U_cpu, W_cpu, L0_cpu, L1_cpu, L2_cpu, L3_cpu, L4_cpu =
        on_architecture(CPU(), (implicit.temp_profile, implicit.S, implicit.S⁻¹, implicit.L,
                                implicit.R, implicit.U, implicit.W, implicit.L0, implicit.L1,
                                implicit.L2, implicit.L3, implicit.L4))

    σ_levels_full_cpu, σ_levels_thick_cpu, Δp_geopot_half_cpu, Δp_geopot_full_cpu, σ_lnp_A_cpu, σ_lnp_B_cpu, temp_average_cpu =
        on_architecture(CPU(), (σ_levels_full, σ_levels_thick, Δp_geopot_half, Δp_geopot_full,
                                σ_lnp_A, σ_lnp_B, diagn.temp_average))

    temp_profile_cpu .= temp_average_cpu
    all(isfinite.(temp_profile_cpu)) || return nothing

    ξ = α*dt
    implicit.ξ[] = ξ

    @inbounds for k in 1:nlayers
        R_cpu[1:k, k] .= -Δp_geopot_full_cpu[k]
        R_cpu[1:k-1, k] .+= -Δp_geopot_half_cpu[k]
    end
    U_cpu .= -R_dry*temp_profile_cpu

    L0_cpu .= 1 ./ 2σ_levels_thick_cpu
    L2_cpu .= κ*temp_profile_cpu.*σ_lnp_A_cpu
    L4_cpu .= κ*temp_profile_cpu.*σ_lnp_B_cpu

    @inbounds for k in 1:nlayers
        Tₖ = temp_profile_cpu[k]
        k_above = max(1, k-1)
        k_below = min(k+1, nlayers)
        ΔT_above = Tₖ - temp_profile_cpu[k_above]
        ΔT_below = temp_profile_cpu[k_below] - Tₖ
        σₖ = σ_levels_full_cpu[k]
        σₖ_above = σ_levels_full_cpu[k_above]

        for r in 1:nlayers
            L1_cpu[k, r]  = ΔT_below*σ_levels_thick_cpu[r]*σₖ
            L1_cpu[k, r] -= k>=r ? σ_levels_thick_cpu[r] : zero(NF)
            L1_cpu[k, r] += ΔT_above*σ_levels_thick_cpu[r]*σₖ_above
            L1_cpu[k, r] -= (k-1)>=r ? σ_levels_thick_cpu[r] : zero(NF)
        end

        L3_cpu[1:k, k]    .= 0
        L3_cpu[k+1:end, k] .= σ_levels_thick_cpu[k]
    end

    L_cpu .= Diagonal(L0_cpu)*L1_cpu .+ Diagonal(L2_cpu)*L3_cpu .+ Diagonal(L4_cpu)
    W_cpu .= -σ_levels_thick_cpu

    I_mat = LinearAlgebra.I(nlayers)
    @inbounds for l in 1:trunc+1
        eigenvalue = -l*(l-1)
        S_cpu .= I_mat .- ξ^2*eigenvalue*(R_cpu*L_cpu .+ U_cpu*W_cpu')
        luS = LinearAlgebra.lu!(S_cpu)
        Sinv = L1_cpu
        Sinv .= I_mat
        LinearAlgebra.ldiv!(luS, Sinv)
        S⁻¹_cpu[l, :, :] .= Sinv
    end

    implicit.temp_profile .= on_architecture(arch, temp_profile_cpu)
    implicit.S    .= on_architecture(arch, S_cpu)
    implicit.S⁻¹  .= on_architecture(arch, S⁻¹_cpu)
    implicit.L    .= on_architecture(arch, L_cpu)
    implicit.R    .= on_architecture(arch, R_cpu)
    implicit.U    .= on_architecture(arch, U_cpu)
    implicit.W    .= on_architecture(arch, W_cpu)
    implicit.L0   .= on_architecture(arch, L0_cpu)
    implicit.L1   .= on_architecture(arch, L1_cpu)
    implicit.L2   .= on_architecture(arch, L2_cpu)
    implicit.L3   .= on_architecture(arch, L3_cpu)
    implicit.L4   .= on_architecture(arch, L4_cpu)
end

set_initialized!(implicit::ImplicitPrimitiveEquation) = (implicit.initialized = true)

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
