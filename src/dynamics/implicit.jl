abstract type AbstractImplicit <: AbstractModelComponent end

# model.implicit = nothing (for BarotropicModel)
initialize!(::Nothing, dt::Real, ::Any...) = nothing
lorenz_implicit_correction!(::DiagnosticVariables, ::PrognosticVariables, ::Nothing, ::AbstractModel) = nothing
_strip_implicit!(::DiagnosticVariables, ::PrognosticVariables, ::Nothing, ::AbstractModel) = nothing

"""$(TYPEDSIGNATURES)
No-op strip for models without implicit terms (e.g. BarotropicModel)."""
lorenz_strip_implicit!(::DiagnosticVariables, ::PrognosticVariables, ::Nothing, ::AbstractModel) = nothing

# ============================================================================
# SHALLOW WATER MODEL
# ============================================================================
export ImplicitShallowWater

"""Struct that holds various precomputed arrays for the semi-implicit correction to
prevent gravity waves from amplifying in the shallow water model.

The centering parameter `α` controls the implicit evaluation:

    α = 0.5  Crank-Nicolson: 2nd-order accurate, slows but does not damp gravity waves
    α = 1.0  Backward Euler: 1st-order accurate, strongly damps gravity waves

Unlike the leapfrog version (where ξ = α·2Δt), for the Lorenz N-cycle ξ = α·Δt,
using only a single time level (Hotta et al. 2016, Eq. 20).

Fields are
$(TYPEDFIELDS)"""
@kwdef mutable struct ImplicitShallowWater{NF} <: AbstractImplicit
    "[OPTION] coefficient for semi-implicit computations to filter gravity waves, 0.5 <= α <= 1"
    α::NF = 1

    "Implicit time step ξ = α·Δt [s], set in initialize!"
    time_step::NF = 0
end

"""$(TYPEDSIGNATURES)
Generator using the resolution from `spectral_grid`."""
ImplicitShallowWater(SG::SpectralGrid; kwargs...) = ImplicitShallowWater{SG.NF}(; kwargs...)

"""$(TYPEDSIGNATURES)
Update the implicit time step ξ = α·Δt for the Lorenz N-cycle shallow water model.
`dt` is the radius-scaled Δt [s/m]. Extra args are accepted and ignored so that
the same call site works for both ShallowWater and PrimitiveEquation."""
function initialize!(implicit::ImplicitShallowWater, dt::Real, args...)
    implicit.time_step = implicit.α * dt    # ξ = α·Δt  (Hotta et al. 2016, Eq. 20)
end

# ImplicitShallowWater has no precomputed arrays requiring initialization tracking
set_initialized!(implicit::ImplicitShallowWater) = nothing
set_initialized!(implicit::Nothing) = nothing

# Strip L_I(x) from diagn.tendencies in-place before Lorenz accumulation so that
# accumulate_tendency! stores only F_E into G (Hotta et al. 2016, Eq. 19 uses F_E).
# For ShallowWater: L_I_div = -g·∇²·η, L_I_pres = -H·D (both already in diagn.tendencies
# from bernoulli_potential_swm! and volume_flux_divergence!).
function _strip_implicit!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    ::ImplicitShallowWater,
    model::ShallowWater,
)
    div_tend  = diagn.tendencies.div_tend
    pres_tend = diagn.tendencies.pres_tend
    div_x     = get_step(progn.div,  1)
    pres_x    = get_step(progn.pres, 1)
    H = model.atmosphere.layer_thickness
    g = model.planet.gravity
    l_indices = div_tend.spectrum.l_indices
    arch = architecture(div_tend)
    launch!(arch, SpectralWorkOrder, size(div_tend),
            _strip_implicit_swm_kernel!,
            div_tend, pres_tend, div_x, pres_x, l_indices, H, g)
    return nothing
end

# F_E_div  = F_div  - (-g·∇²·η) = F_div  + g·∇²·η  (∇²≤0: removes positive gravity wave term)
# F_E_pres = F_pres - (-H·D)    = F_pres + H·D
@kernel inbounds=true function _strip_implicit_swm_kernel!(
    div_tend, pres_tend,
    @Const(div_x), @Const(pres_x),
    l_indices,
    @Const(H), @Const(g),
)
    I  = @index(Global, Cartesian)
    lm = I[1]
    k  = I[2]
    l  = l_indices[lm]
    ∇² = -l*(l-1)
    div_tend[lm, k] += g * ∇² * pres_x[lm]
    pres_tend[lm]   += H * div_x[lm, k]
end

"""$(TYPEDSIGNATURES)
Apply the semi-implicit correction for the Lorenz N-cycle shallow water model,
implementing Hotta et al. (2016) Eq. 20:

    δx = (I - α·Δt·L_I)⁻¹ · (G_E + L_I·x)

Reads G_E from progn index 2 (the F_E accumulator, left unchanged) and the current
state x from progn index 1. Writes the corrected δx into diagn.tendencies.div_tend
and diagn.tendencies.pres_tend, which were pre-filled with G_E by the caller via a
G→tend copy before this function is called."""
function lorenz_implicit_correction!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    implicit::ImplicitShallowWater,
    model::ShallowWater,
)
    div_G    = get_step(progn.div,  2)          # F_E accumulator (read only)
    pres_G   = get_step(progn.pres, 2)          # F_E accumulator (read only)
    div_x    = get_step(progn.div,  1)          # current state
    pres_x   = get_step(progn.pres, 1)          # current state
    div_δx   = diagn.tendencies.div_tend        # output: corrected δD
    pres_δx  = diagn.tendencies.pres_tend       # output: corrected δη

    H = model.atmosphere.layer_thickness
    g = model.planet.gravity
    ξ = implicit.time_step

    l_indices = div_G.spectrum.l_indices
    arch = architecture(div_G)

    launch!(arch, SpectralWorkOrder, size(div_G),
            lorenz_implicit_shallow_water_kernel!,
            div_G, pres_G, div_x, pres_x, div_δx, pres_δx, l_indices, H, g, ξ)

    zero_last_degree!(div_δx)
    zero_last_degree!(pres_δx)

    return nothing
end

@kernel inbounds=true function lorenz_implicit_shallow_water_kernel!(
    @Const(div_G), @Const(pres_G),
    @Const(div_x), @Const(pres_x),
    div_δx, pres_δx,
    l_indices,
    @Const(H), @Const(g), @Const(ξ),
)
    I  = @index(Global, Cartesian)
    lm = I[1]
    k  = I[2]

    l  = l_indices[lm]
    ∇² = -l*(l-1)   # Laplacian eigenvalue, 1-based, always ≤ 0

    # Form exact Eq. 20 RHS = G_E + L_I(x_current):
    #   L_I_div  = -g·∇²·η  (≥0 since ∇²≤0)
    #   L_I_pres = -H·D
    G_div = div_G[lm, k] - g*∇²*pres_x[lm]
    G_η   = pres_G[lm]   - H*div_x[lm, k]

    # Scalar inversion: S⁻¹ = 1/(1 - ξ²·H·g·∇²). Since ∇²≤0, S⁻¹≤1.
    S⁻¹ = inv(1 - ξ^2*H*g*∇²)

    # Solve 2×2 system (I - ξ·L_I)·δx = G:
    #   δD + ξ·g·∇²·δη = G_div  =>  δD = S⁻¹·(G_div - ξ·g·∇²·G_η)
    #   δη + ξ·H·δD    = G_η    =>  δη = G_η - ξ·H·δD
    δD = S⁻¹*(G_div - ξ*g*∇²*G_η)
    div_δx[lm, k] = δD
    pres_δx[lm]   = G_η - ξ*H*δD
end

# ============================================================================
# PRIMITIVE EQUATION MODEL
# ============================================================================
export ImplicitPrimitiveEquation

"""Struct that holds various precomputed arrays for the semi-implicit correction to
prevent gravity waves from amplifying in the primitive equation model.

The operator S⁻¹ is precomputed for each spherical harmonic degree l and stored
as a (trunc+2) × nlayers × nlayers tensor. For the Lorenz N-cycle, ξ = α·Δt
(not α·2Δt as in leapfrog).

$(TYPEDFIELDS)"""
@kwdef mutable struct ImplicitPrimitiveEquation{
    NF,
    VectorType,
    MatrixType,
    TensorType,
} <: AbstractImplicit

    # DIMENSIONS
    "Spectral resolution"
    trunc::Int

    "Number of vertical layers"
    nlayers::Int

    # PARAMETERS
    "Time-step coefficient: 0=explicit, 0.5=centred implicit (Crank-Nicolson), 1=backward implicit"
    α::NF = 1

    "Reinitialize at restart when initialized=true"
    reinitialize::Bool = true

    "Flag automatically set to true when initialize! has been called"
    initialized::Bool = false

    # PRECOMPUTED ARRAYS, initialized with initialize!
    "Vertical temperature profile, obtained from diagn on first time step"
    temp_profile::VectorType = zeros(NF, nlayers)

    "Implicit time step ξ = α·Δt packed in RefValue for mutability"
    ξ::Base.RefValue{NF} = Ref{NF}(0)

    "Divergence: operator for the geopotential calculation"
    R::MatrixType = zeros(NF, nlayers, nlayers)

    "Divergence: the -RdTₖ∇² term excl the eigenvalues from ∇² for divergence"
    U::VectorType = zeros(NF, nlayers)

    "Temperature: operator for the TₖD + κTₖDlnps/Dt term"
    L::MatrixType = zeros(NF, nlayers, nlayers)

    "Pressure: vertical averaging of the -D̄ term in the log surface pres equation"
    W::VectorType = zeros(NF, nlayers)

    "Components to construct L, 1/2Δσ"
    L0::VectorType = zeros(NF, nlayers)

    "Vertical advection term in the temperature equation (below+above)"
    L1::MatrixType = zeros(NF, nlayers, nlayers)

    "Factor in front of the div_sum_above term"
    L2::VectorType = zeros(NF, nlayers)

    "_sum_above operator itself"
    L3::MatrixType = zeros(NF, nlayers, nlayers)

    "Factor in front of div term in Dlnps/Dt"
    L4::VectorType = zeros(NF, nlayers)

    "For every l the matrix to be inverted"
    S::MatrixType = zeros(NF, nlayers, nlayers)

    "Combined inverted operator: S⁻¹ = (I - ξ²(RL + UW))⁻¹"
    S⁻¹::TensorType = zeros(NF, trunc+2, nlayers, nlayers)
end

"""$(TYPEDSIGNATURES)
Generator using the resolution from SpectralGrid."""
function ImplicitPrimitiveEquation(spectral_grid::SpectralGrid; kwargs...)
    (; NF, VectorType, MatrixType, TensorType, trunc, nlayers) = spectral_grid
    return ImplicitPrimitiveEquation{NF, VectorType, MatrixType, TensorType}(;
        trunc, nlayers, kwargs...)
end

# Function barrier to unpack the model for primitive equation models
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

"""$(TYPEDSIGNATURES)
Initialize the implicit terms for the PrimitiveEquation models.
Precomputes S⁻¹[l, :, :] = (I - ξ²·(R·L + U·W'))⁻¹ for each spherical harmonic
degree l, where ξ = α·Δt for the Lorenz N-cycle (not α·2Δt as in leapfrog)."""
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

    # Transfer small arrays to CPU for precomputation
    temp_profile_cpu, S_cpu, S⁻¹_cpu, L_cpu, R_cpu, U_cpu, W_cpu,
    L0_cpu, L1_cpu, L2_cpu, L3_cpu, L4_cpu =
        on_architecture(CPU(), (implicit.temp_profile, implicit.S, implicit.S⁻¹,
                                implicit.L, implicit.R, implicit.U, implicit.W,
                                implicit.L0, implicit.L1, implicit.L2, implicit.L3, implicit.L4))

    σ_levels_full_cpu, σ_levels_thick_cpu, Δp_geopot_half_cpu, Δp_geopot_full_cpu,
    σ_lnp_A_cpu, σ_lnp_B_cpu, temp_average_cpu =
        on_architecture(CPU(), (σ_levels_full, σ_levels_thick, Δp_geopot_half, Δp_geopot_full,
                                σ_lnp_A, σ_lnp_B, diagn.temp_average))

    temp_profile_cpu .= temp_average_cpu
    all(isfinite.(temp_profile_cpu)) || return nothing  # bail if model has blown up

    # ξ = α·Δt  (Lorenz N-cycle: single Δt, not 2Δt as in leapfrog)
    ξ = α * dt
    implicit.ξ[] = ξ

    # DIVERGENCE OPERATORS (R, U)
    @inbounds for k in 1:nlayers
        R_cpu[1:k, k]   .= -Δp_geopot_full_cpu[k]
        R_cpu[1:k-1, k] .+= -Δp_geopot_half_cpu[k]
    end
    U_cpu .= -R_dry * temp_profile_cpu

    # TEMPERATURE OPERATOR (L)
    L0_cpu .= 1 ./ (2 .* σ_levels_thick_cpu)
    L2_cpu .= κ .* temp_profile_cpu .* σ_lnp_A_cpu
    L4_cpu .= κ .* temp_profile_cpu .* σ_lnp_B_cpu

    @inbounds for k in 1:nlayers
        Tₖ       = temp_profile_cpu[k]
        k_above  = max(1, k-1)
        k_below  = min(k+1, nlayers)
        ΔT_above = Tₖ - temp_profile_cpu[k_above]
        ΔT_below = temp_profile_cpu[k_below] - Tₖ
        σₖ       = σ_levels_full_cpu[k]
        σₖ_above = σ_levels_full_cpu[k_above]

        for r in 1:nlayers
            L1_cpu[k, r]  = ΔT_below * σ_levels_thick_cpu[r] * σₖ
            L1_cpu[k, r] -= k >= r ? σ_levels_thick_cpu[r] : zero(NF)
            L1_cpu[k, r] += ΔT_above * σ_levels_thick_cpu[r] * σₖ_above
            L1_cpu[k, r] -= (k-1) >= r ? σ_levels_thick_cpu[r] : zero(NF)
        end

        L3_cpu[1:k, k]     .= 0
        L3_cpu[k+1:end, k] .= σ_levels_thick_cpu[k]
    end

    L_cpu .= Diagonal(L0_cpu)*L1_cpu .+ Diagonal(L2_cpu)*L3_cpu .+ Diagonal(L4_cpu)

    # PRESSURE OPERATOR (W)
    W_cpu .= -σ_levels_thick_cpu

    # INVERT S = I - ξ²·(R·L + U·W') for each spherical harmonic degree l
    Id = LinearAlgebra.I(nlayers)
    @inbounds for l in 1:trunc+1
        eigenvalue = -l*(l-1)   # 1-based Laplacian eigenvalue: always ≤ 0
        S_cpu .= Id .- ξ^2 .* eigenvalue .* (R_cpu*L_cpu .+ U_cpu*W_cpu')

        luS  = LinearAlgebra.lu!(S_cpu)
        Sinv = L1_cpu       # reuse L1 as scratch space
        Sinv .= Id
        LinearAlgebra.ldiv!(luS, Sinv)
        S⁻¹_cpu[l, :, :] .= Sinv
    end

    # Transfer back to target architecture
    implicit.temp_profile .= on_architecture(arch, temp_profile_cpu)
    implicit.S            .= on_architecture(arch, S_cpu)
    implicit.S⁻¹          .= on_architecture(arch, S⁻¹_cpu)
    implicit.L            .= on_architecture(arch, L_cpu)
    implicit.R            .= on_architecture(arch, R_cpu)
    implicit.U            .= on_architecture(arch, U_cpu)
    implicit.W            .= on_architecture(arch, W_cpu)
    implicit.L0           .= on_architecture(arch, L0_cpu)
    implicit.L1           .= on_architecture(arch, L1_cpu)
    implicit.L2           .= on_architecture(arch, L2_cpu)
    implicit.L3           .= on_architecture(arch, L3_cpu)
    implicit.L4           .= on_architecture(arch, L4_cpu)

    return nothing
end

set_initialized!(implicit::ImplicitPrimitiveEquation) = (implicit.initialized = true)

"""$(TYPEDSIGNATURES)
Apply the semi-implicit correction for the Lorenz N-cycle primitive equation model,
implementing Hotta et al. (2016) Eq. 20:

    dx = (I - α·Δt·L_I)⁻¹ · (G + L_I·x)

The implicit operator L_I covers:
  - divergence D:           geopotential and surface pressure gradient
  - temperature T:          adiabatic expansion (L operator)
  - log surface pressure:   mass divergence (W operator)

Unlike the leapfrog version, L_I·x is evaluated at the current state x (index 1)
rather than xⁿ⁻¹, so there is no (div_old - div_new) differencing."""
function lorenz_implicit_correction!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    implicit::ImplicitPrimitiveEquation,
    model::PrimitiveEquation,
)
    implicit.α == 0 && return nothing   # skip if fully explicit

    (; nlayers) = implicit
    (; S⁻¹, R, U, L, W) = implicit
    ξ = implicit.ξ[]

    (; temp_tend, pres_tend, div_tend) = diagn.tendencies

    # Current state x at index 1 (no old time level)
    div_current = get_step(progn.div, 1)

    G            = diagn.dynamics.a     # work array: combined tendency G_D + ξRG_T + ξUG_lnps
    geopotential = diagn.dynamics.b     # work array: R·G_T term
    geopotential .= 0

    l_indices = temp_tend.spectrum.l_indices
    arch = architecture(temp_tend)

    lm_size = size(pres_tend, 1)
    launch!(arch, LinearWorkOrder, (lm_size,),
            lorenz_implicit_primitive_kernel!,
            temp_tend, pres_tend, div_tend, G, geopotential,
            div_current, S⁻¹, R, U, L, W, l_indices, ξ, nlayers)

    zero_last_degree!(div_tend)
    zero_last_degree!(pres_tend)
    zero_last_degree!(temp_tend)

    pres_tend.data[1:1] .= 0    # mass conservation

    return nothing
end

@kernel inbounds=true function lorenz_implicit_primitive_kernel!(
    temp_tend, pres_tend, div_tend, G, geopotential,
    div_current,
    @Const(S⁻¹), @Const(R), @Const(U), @Const(L), @Const(W), @Const(l_indices),
    @Const(ξ), @Const(nlayers),
)
    lm = @index(Global, Linear)

    l          = l_indices[lm]
    eigenvalue = -l*(l-1)   # Laplacian eigenvalue, 1-based, always ≤ 0

    # Step 1: Add L_I·x to temperature tendency
    # In leapfrog: L*(div_old - div_new) shifts implicit term to old level
    # Here: L*div_current evaluates directly at current state x
    for k in 1:nlayers
        temp_correction = zero(eltype(temp_tend))
        for r in 1:nlayers
            temp_correction += L[k, r] * div_current[lm, r]
        end
        temp_tend[lm, k] += temp_correction
    end

    # Step 2: Vertical integration of geopotential  R·G_T
    for k in 1:nlayers
        geopotential_val = zero(eltype(geopotential))
        for r in k:nlayers
            geopotential_val += R[k, r] * temp_tend[lm, r]
        end
        geopotential[lm, k] = geopotential_val
    end

    # Step 3: Form G = G_D + ξ·eigenvalue·(U·G_lnps + R·G_T)
    for k in 1:nlayers
        G[lm, k] = div_tend[lm, k] + ξ*eigenvalue*(U[k]*pres_tend[lm] + geopotential[lm, k])
    end

    # Step 4: Solve δD = S⁻¹·G
    for k in 1:nlayers
        div_val = zero(eltype(div_tend))
        for r in 1:nlayers
            div_val += S⁻¹[l, k, r] * G[lm, r]
        end
        div_tend[lm, k] = div_val
    end

    # Step 5a: Temperature correction  δT = G_T + ξ·L·δD
    for k in 1:nlayers
        temp_correction = zero(eltype(temp_tend))
        for r in 1:nlayers
            temp_correction += ξ * L[k, r] * div_tend[lm, r]
        end
        temp_tend[lm, k] += temp_correction
    end

    # Step 5b: Pressure correction  δlnpₛ = G_lnpₛ + ξ·W·δD
    pres_correction = zero(eltype(pres_tend))
    for k in 1:nlayers
        pres_correction += ξ * W[k] * div_tend[lm, k]
    end
    pres_tend[lm] += pres_correction
end