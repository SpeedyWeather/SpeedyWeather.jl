abstract type AbstractImplicit <: AbstractModelComponent end

# For barotrpoic model 
initialize!(::Nothing, dt::Real, α::Real, ::AbstractModel) = nothing
lorenz_implicit_correction!(::DiagnosticVariables, ::PrognosticVariables, ::Nothing, ::AbstractModel) = nothing
set_initialized!(::Nothing) = nothing

# Shallow water implicit 
export ImplicitShallowWater

"""
Struct that holds parameters for the semi-implicit correction in shallow water equations.

The implicit correction prevents gravity waves from amplifying. The parameter `α` controls
the centering of the implicit scheme:
    α = 0.5: Crank-Nicolson (2nd order, centered implicit)
    α = 1.0: Backward Euler (1st order)
    α ∈ [0.5, 1]: Controls gravity wave dampening strength

Fields are:
$(TYPEDFIELDS)
"""
@kwdef mutable struct ImplicitShallowWater{NF} <: AbstractImplicit
    "[OPTION] Coefficient for semi-implicit computations, 0.5 <= α <= 1"
    α::NF = 0.5
    
    "[DERIVED] Time step ξ = α*Δt (set during initialization)"
    time_step::NF = 0
end

"""$(TYPEDSIGNATURES)
Generator using the resolution from `spectral_grid`."""
ImplicitShallowWater(SG::SpectralGrid; kwargs...) = ImplicitShallowWater{SG.NF}(; kwargs...)

"""$(TYPEDSIGNATURES)
Initialize implicit solver with time step Δt and centering parameter α.

For Lorenz N-cycle: ξ = α*Δt (not 2α*Δt like Leapfrog)
"""
function initialize!(implicit::ImplicitShallowWater, dt::Real, α::Real, args...)
    implicit.α = α
    implicit.time_step = α * dt  # ξ = α*Δt for Lorenz N-cycle
end

set_initialized!(implicit::ImplicitShallowWater) = nothing


"""$(TYPEDSIGNATURES)
Apply implicit correction for Lorenz N-cycle (shallow water).

Algorithm (Hotta et al. 2016, Eq. 20):
1. Form complete tendencies: G_total = G_accumulated + L_I*x_current
2. Apply implicit operator inversion: (I - ξ²Hg∇²)^(-1)
3. Update G arrays with corrected tendencies
Apply semi-implicit correction to the accumulated tendencies for the Lorenz N-cycle
time stepping scheme.

This is the direct analog of the leapfrog implicit correction, but rewritten
for a single-time-level formulation following Hotta et al. (2016).

The governing equations are:

    ∂D/∂t = N_D + g ∇² η
    ∂η/∂t = N_η + H D

where the gravity-wave terms are treated semi-implicitly.

The Lorenz N-cycle implicit step solves:

    (I - ξ² g H ∇²) δD = G_D + g ∇² η
    δη = G_η + H δD

with ξ = α Δt.
"""
function lorenz_implicit_correction!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    implicit::ImplicitShallowWater,
    model::ShallowWater,
)

    # Accumulated tendencies G from the N-cycle
    (; div_tend, pres_tend) = diagn.tendencies

    # Current prognostic state (single time level!)
    div = get_step(progn.div, 1)
    pres = get_step(progn.pres, 1)

    # Physical parameters
    H = model.atmosphere.layer_thickness   # layer thickness [m]
    g = model.planet.gravity               # gravity [m/s²]

    # Implicit time step ξ = α Δt
    ξ = implicit.time_step

    # Spectral indices
    l_indices = div_tend.spectrum.l_indices

    # GPU / CPU kernel launch
    arch = architecture(div_tend)
    launch!(
        arch,
        SpectralWorkOrder,
        size(div_tend),
        lorenz_implicit_shallow_water_kernel!,
        div_tend,
        pres_tend,
        div,
        pres,
        l_indices,
        H,
        g,
        ξ,
    )

    zero_last_degree!(div_tend)
    zero_last_degree!(pres_tend)

    return nothing
end


@kernel inbounds=true function lorenz_implicit_shallow_water_kernel!(
    div_tend, pres_tend, div, pres, l_indices,
    @Const(H), @Const(g), @Const(ξ)
)
    I = @index(Global, Cartesian)
    lm = I[1]
    k  = I[2]

    # Spherical harmonic degree l
    l = l_indices[lm]

    # Laplacian eigenvalue (dimensionless, radius-scaled)
    ∇² = -l*(l-1)

    # Form the total tendencies G + L_I x 
    #   G_D = G_D + g∇² η
    #   G_η = G_η + H D
    G_div = div_tend[lm, k] + g * ∇² * pres[lm]
    G_η   = pres_tend[lm]  + H * div[lm, k]

    # Implicit solve:
    #   (1 - ξ² g H ∇²) δD = G_div - ξ g ∇² G_η
    #
    @assert isfinite(1 - ξ^2 * H * g * ∇²)

    S⁻¹ = inv(1 - ξ^2*H*g*∇²)

    δdiv = S⁻¹*(G_div-ξ*g*∇²*G_η)

    #Back-substitution for η
    δη = G_η-ξ*H*δdiv

    # Store corrected tendencies
    div_tend[lm, k] = δdiv
    pres_tend[lm]   = δη
end


# ============================================================================
# PRIMITIVE EQUATION IMPLICIT
# ============================================================================

export ImplicitPrimitiveEquation

"""
Struct that holds precomputed arrays for the semi-implicit correction in primitive equations.

The implicit correction prevents gravity waves from amplifying. Operators are precomputed
for efficiency and applied at each timestep.

$(TYPEDFIELDS)
"""
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
    "Time-step coefficient: 0.5=Crank-Nicolson (2nd order), 1.0=Backward Euler (1st order)"
    α::NF = 0.5

    "Reinitialize at restart when initialized=true"
    reinitialize::Bool = true

    "Flag automatically set to true when initialize! has been called"
    initialized::Bool = false

    # PRECOMPUTED ARRAYS (initialized in initialize!)
    "Vertical temperature profile, obtained from diagn on first time step"
    temp_profile::VectorType = zeros(NF, nlayers)

    "Time step ξ = α*Δt packed in RefValue for mutability"
    ξ::Base.RefValue{NF} = Ref{NF}(0)

    "Divergence: operator for geopotential calculation"
    R::MatrixType = zeros(NF, nlayers, nlayers)

    "Divergence: the -RdTₖ∇² term excl eigenvalues for divergence"
    U::VectorType = zeros(NF, nlayers)

    "Temperature: operator for TₖD + κTₖDlnps/Dt term"
    L::MatrixType = zeros(NF, nlayers, nlayers)

    "Pressure: vertical averaging of -D̄ term in log surface pressure equation"
    W::VectorType = zeros(NF, nlayers)

    # Components to construct L
    "1/(2Δσ)"
    L0::VectorType = zeros(NF, nlayers)

    "Vertical advection term (below+above)"
    L1::MatrixType = zeros(NF, nlayers, nlayers)

    "Factor in front of div_sum_above term"
    L2::VectorType = zeros(NF, nlayers)

    "Sum_above operator"
    L3::MatrixType = zeros(NF, nlayers, nlayers)

    "Factor in front of div term in Dlnps/Dt"
    L4::VectorType = zeros(NF, nlayers)

    # Inversion matrices
    "Temporary matrix for inversion"
    S::MatrixType = zeros(NF, nlayers, nlayers)

    "Combined inverted operator: S⁻¹ = (I - ξ²*∇²*(RL + UW'))^(-1)"
    S⁻¹::TensorType = zeros(NF, trunc+2, nlayers, nlayers)
end

"""$(TYPEDSIGNATURES)
Generator using the resolution from SpectralGrid."""
function ImplicitPrimitiveEquation(spectral_grid::SpectralGrid; kwargs...)
    (; NF, VectorType, MatrixType, TensorType, trunc, nlayers) = spectral_grid
    return ImplicitPrimitiveEquation{NF, VectorType, MatrixType, TensorType}(;
        trunc, nlayers, kwargs...)
end

# Function barrier to unpack the constants struct for primitive eq models
function initialize!(
    I::ImplicitPrimitiveEquation,
    dt::Real,
    α::Real,
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)
    model.dynamics || return nothing
    (; geometry, geopotential, atmosphere, adiabatic_conversion) = model
    initialize!(I, dt, α, diagn, geometry, geopotential, atmosphere, adiabatic_conversion)
end

"""$(TYPEDSIGNATURES)
Initialize the implicit terms for the PrimitiveEquation models."""
function initialize!(
    implicit::ImplicitPrimitiveEquation,
    dt::Real,
    α::Real,
    diagn::DiagnosticVariables,
    geometry::AbstractGeometry,
    geopotential::AbstractGeopotential,
    atmosphere::AbstractAtmosphere,
    adiabatic_conversion::AbstractAdiabaticConversion,
) 
    NF = eltype(diagn)
    # option to skip reinitialization at restart
    (implicit.initialized && !implicit.reinitialize) && return nothing

    (; trunc, nlayers) = implicit
    (; σ_levels_full, σ_levels_thick) = geometry
    (; R_dry, κ) = atmosphere
    (; Δp_geopot_half, Δp_geopot_full) = geopotential
    (; σ_lnp_A, σ_lnp_B) = adiabatic_conversion

    # Get the architecture to transfer back at the end
    arch = architecture(implicit.S)
    
    # Transfer all arrays that need to be computed to CPU
    # These are small (nlayers × nlayers) matrices, so CPU computation is more efficient
    temp_profile_cpu, S_cpu, S⁻¹_cpu, L_cpu, R_cpu, U_cpu, W_cpu, L0_cpu, L1_cpu, L2_cpu, L3_cpu, L4_cpu = 
        on_architecture(CPU(), (implicit.temp_profile, implicit.S, implicit.S⁻¹, implicit.L, 
                                implicit.R, implicit.U, implicit.W, implicit.L0, implicit.L1, 
                                implicit.L2, implicit.L3, implicit.L4))
    
    # Also transfer geometry and other arrays to CPU
    σ_levels_full_cpu, σ_levels_thick_cpu, Δp_geopot_half_cpu, Δp_geopot_full_cpu, σ_lnp_A_cpu, σ_lnp_B_cpu, temp_average_cpu = 
        on_architecture(CPU(), (σ_levels_full, σ_levels_thick, Δp_geopot_half, Δp_geopot_full, 
                                σ_lnp_A, σ_lnp_B, diagn.temp_average))

    # use current vertical temperature profile
    temp_profile_cpu .= temp_average_cpu

    # return immediately if temp_profile contains NaNs, model blew up in that case
    all(isfinite.(temp_profile_cpu)) || return nothing

    # set up R, U, L, W operators from
    # δD = G_D + ξ(RδT + Uδlnps)        divergence D correction
    # δT = G_T + ξLδD                   temperature T correction
    # δlnps = G_lnps + ξWδD             log surface pressure lnps correction
    #
    # For Lorenz N-cycle (Hotta et al. 2016, Eq. 20):
    # G_X is the accumulated weighted tendency: G = w*F_E(x) + (1-w)*G
    # The implicit terms L_I*x from the current state x are added to G
    # Then G is corrected: G ← (I - ξ*L_I)^(-1) * (G + L_I*x)
    #
    # R, U, L, W are linear operators that define the implicit coupling
    # S = I - ξ²*∇²*(RL + UW') is precomputed and inverted for each wavenumber

    # Update the implicit α and compute ξ
    implicit.α = α
    ξ = α * dt                      # ξ = α*Δt for Lorenz N-cycle (not 2α*Δt like leapfrog)
    implicit.ξ[] = ξ

    # DIVERGENCE OPERATORS (called g in Hoskins and Simmons 1975, eq 11 and Appendix 1)
    @inbounds for k in 1:nlayers                # vertical geopotential integration as matrix operator
        R_cpu[1:k, k] .= -Δp_geopot_full_cpu[k]         # otherwise equivalent to geopotential! with zero orography
        R_cpu[1:k-1, k] .+= -Δp_geopot_half_cpu[k]      # incl the minus but excluding the eigenvalues as with U
    end
    U_cpu .= -R_dry*temp_profile_cpu        # the R_d*Tₖ∇² term excl the eigenvalues from ∇² for divergence

    # TEMPERATURE OPERATOR (called τ in Hoskins and Simmons 1975, eq 9 and Appendix 1)
    L0_cpu .= 1 ./ (2 * σ_levels_thick_cpu)
    L2_cpu .= κ * temp_profile_cpu .* σ_lnp_A_cpu    # factor in front of the div_sum_above term
    L4_cpu .= κ * temp_profile_cpu .* σ_lnp_B_cpu    # factor in front of div term in Dlnps/Dt

    @inbounds for k in 1:nlayers
        Tₖ = temp_profile_cpu[k]                    # average temperature at k
        k_above = max(1, k-1)                       # layer index above
        k_below = min(k+1, nlayers)                 # layer index below
        ΔT_above = Tₖ - temp_profile_cpu[k_above]   # temperature difference to layer above
        ΔT_below = temp_profile_cpu[k_below] - Tₖ   # and to layer below
        σₖ = σ_levels_full_cpu[k]                   # should be Σ_r=1^k Δσᵣ for model top at >0hPa
        σₖ_above = σ_levels_full_cpu[k_above]

        for r in 1:nlayers
            L1_cpu[k, r] = ΔT_below * σ_levels_thick_cpu[r] * σₖ         # vert advection operator below
            L1_cpu[k, r] -= k >= r ? σ_levels_thick_cpu[r] : zero(NF)

            L1_cpu[k, r] += ΔT_above * σ_levels_thick_cpu[r] * σₖ_above   # vert advection operator above
            L1_cpu[k, r] -= (k-1) >= r ? σ_levels_thick_cpu[r] : zero(NF)
        end

        # _sum_above operator itself
        L3_cpu[1:k, k] .= 0                              # fill upper triangle + diagonal with zeros
        L3_cpu[k+1:end, k] .= σ_levels_thick_cpu[k]      # vert integration top to k-1
    end

    L_cpu .= Diagonal(L0_cpu) * L1_cpu .+ Diagonal(L2_cpu) * L3_cpu .+ Diagonal(L4_cpu)  # combine all operators into L

    # PRESSURE OPERATOR (called πᵣ in Hoskins and Simmons, 1975 Appendix 1)
    W_cpu .= -σ_levels_thick_cpu                # the -D̄ term in the log surface pres equation

    # solving the equations above for δD yields
    # δD = S⁻¹*G, with G = G_D + ξ*∇²*(R*G_T + U*G_lnps) and the operator S
    # S = I - ξ²*∇²*(RL + UW') that has to be inverted to obtain δD from the Gs
    I_matrix = LinearAlgebra.I(nlayers)
    @inbounds for l in 1:trunc+1
        eigenvalue = -l*(l-1)           # 1-based, -l*(l+1) → -l*(l-1)
        S_cpu .= I_matrix .- ξ^2 * eigenvalue * (R_cpu * L_cpu .+ U_cpu * W_cpu')

        # inv(S) but saving memory:
        luS = LinearAlgebra.lu!(S_cpu)      # in-place LU decomposition (overwriting S)
        Sinv = L1_cpu                       # reuse L1 matrix to store inv(S)
        Sinv .= I_matrix                    # use ldiv! so last arg needs to be unity matrix
        LinearAlgebra.ldiv!(luS, Sinv)      # now do S\I = S⁻¹ via LU decomposition
        S⁻¹_cpu[l, :, :] .= Sinv            # store in array
    end
    
    # Transfer computed results back to the original architecture
    implicit.temp_profile .= on_architecture(arch, temp_profile_cpu)
    implicit.S .= on_architecture(arch, S_cpu)
    implicit.S⁻¹ .= on_architecture(arch, S⁻¹_cpu)
    implicit.L .= on_architecture(arch, L_cpu)
    implicit.R .= on_architecture(arch, R_cpu)
    implicit.U .= on_architecture(arch, U_cpu)
    implicit.W .= on_architecture(arch, W_cpu)
    implicit.L0 .= on_architecture(arch, L0_cpu)
    implicit.L1 .= on_architecture(arch, L1_cpu)
    implicit.L2 .= on_architecture(arch, L2_cpu)
    implicit.L3 .= on_architecture(arch, L3_cpu)
    implicit.L4 .= on_architecture(arch, L4_cpu)
end

set_initialized!(implicit::ImplicitPrimitiveEquation) = (implicit.initialized = true)

"""$(TYPEDSIGNATURES)
Apply implicit correction for Lorenz N-cycle (primitive equations).

Following Hotta et al. (2016), Equation 20:
    dx = (I - ξ*L_I)^(-1) * (G + L_I*x)

where G is the accumulated weighted tendency and x is the current state.

Key difference from Leapfrog:
- Leapfrog: Corrects using (x^(i-1) - x^(i+1)), the difference between time steps
- Lorenz N-cycle: Adds implicit terms from current state x^(i) to G, then applies correction

Algorithm:
1. Compute G_total = G_accumulated + L_I*x_current
2. Solve for corrected divergence: δD = S⁻¹*G_total  
3. Back-substitute to get corrected temperature and pressure tendencies
"""
function lorenz_implicit_correction!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    implicit::ImplicitPrimitiveEquation,
    model::PrimitiveEquation,
)
    # escape immediately if explicit
    implicit.α == 0 && return nothing

    (; nlayers, trunc) = implicit
    (; S⁻¹, R, U, L, W) = implicit
    ξ = implicit.ξ[]

    (; temp_tend, pres_tend, div_tend) = diagn.tendencies
    div_x, _ = get_steps(progn.div)   # Current state (index 1)
    
    G = diagn.dynamics.a              # Combined tendency workspace
    geopot_from_Gtemp = diagn.dynamics.b   # R*G_temp workspace

    l_indices = temp_tend.spectrum.l_indices

    arch = architecture(temp_tend)
    lm_size = size(pres_tend, 1)
    
    # Launch kernel for all spectral modes
    launch!(arch, LinearWorkOrder, (lm_size,),
            lorenz_implicit_primitive_kernel!,
            temp_tend, pres_tend, div_tend, G, geopot_from_Gtemp,
            div_x, S⁻¹, R, U, L, W, l_indices, ξ, nlayers)

    zero_last_degree!(div_tend)
    zero_last_degree!(pres_tend)
    zero_last_degree!(temp_tend)
    pres_tend.data[1:1] .= 0    # Mass conservation
    
    return nothing
end

@kernel inbounds=true function lorenz_implicit_primitive_kernel!(
    temp_tend, pres_tend, div_tend, G, geopot_from_Gtemp,
    div_x, @Const(S⁻¹), @Const(R), @Const(U), @Const(L), @Const(W), 
    @Const(l_indices), @Const(ξ), @Const(nlayers)
)
    lm = @index(Global, Linear)
    
    l = l_indices[lm]
    eigenvalue = -l*(l-1)
    
    # ============================================================================
    # Step 1: Compute geopotential from ORIGINAL temperature tendency
    # ============================================================================
    # geopot = R*G_temp (using accumulated temp_tend, before adding implicit terms)
    # This is the same calculation as in leapfrog, but temp_tend is already the
    # accumulated weighted tendency G, not a difference between time steps
    for k in 1:nlayers
        geopot_val = zero(eltype(geopot_from_Gtemp))
        for r in k:nlayers  # Integration from k to surface
            geopot_val += R[k, r] * temp_tend[lm, r]
        end
        geopot_from_Gtemp[lm, k] = geopot_val
    end
    
    # ============================================================================
    # Step 2: Form combined divergence tendency
    # ============================================================================
    # For Lorenz N-cycle: G = G_div + ∇²*U*div_x + ξ*∇²*(U*G_pres + R*G_temp)
    # 
    # Contrast with Leapfrog: G_div + ∇²*U*(pres_old - pres_new) + ...
    # 
    # The key difference: we add the implicit term ∇²*U*div_x from the CURRENT state
    # rather than a difference (pres_old - pres_new)
    for k in 1:nlayers
        G[lm, k] = div_tend[lm, k] + 
                   eigenvalue * U[k] * div_x[lm, k] +  # L_I*x term from current state
                   ξ * eigenvalue * (U[k] * pres_tend[lm] + geopot_from_Gtemp[lm, k])
    end
    
    # ============================================================================
    # Step 3: Solve for corrected divergence: δD = S⁻¹*G
    # ============================================================================
    # This step is identical to leapfrog
    # Store temporarily in geopot_from_Gtemp (reuse workspace)
    for k in 1:nlayers
        corrected_div = zero(eltype(div_tend))
        for r in 1:nlayers
            corrected_div += S⁻¹[l, k, r] * G[lm, r]
        end
        geopot_from_Gtemp[lm, k] = corrected_div  # Temporary storage
    end
    
    # Copy corrected divergence back
    for k in 1:nlayers
        div_tend[lm, k] = geopot_from_Gtemp[lm, k]
    end
    
    # ============================================================================
    # Step 4: Back-substitute for temperature
    # ============================================================================
    # For Lorenz N-cycle: δT = G_temp + L*div_x + ξ*L*δD
    #
    # Contrast with Leapfrog: δT = G_temp + L*(div_old - div_new) + ξ*L*δD
    #
    # Again, we add L*div_x from the CURRENT state, not a difference
    for k in 1:nlayers
        temp_correction = zero(eltype(temp_tend))
        
        # Add L*div_x (implicit term from current state)
        for r in 1:nlayers
            temp_correction += L[k, r] * div_x[lm, r]
        end
        
        # Add ξ*L*δD (coupling to corrected divergence)
        for r in 1:nlayers
            temp_correction += ξ * L[k, r] * div_tend[lm, r]
        end
        
        temp_tend[lm, k] += temp_correction
    end
    
    # ============================================================================
    # Step 5: Back-substitute for pressure
    # ============================================================================
    # For Lorenz N-cycle: δln(pₛ) = G_pres + W*div_x + ξ*W*δD
    #
    # Contrast with Leapfrog: δln(pₛ) = G_pres + W*(div_old - div_new) + ξ*W*δD
    pres_correction = zero(eltype(pres_tend))
    
    for k in 1:nlayers
        # W*div_x (implicit term from current state)
        pres_correction += W[k] * div_x[lm, k]
        
        # ξ*W*δD (coupling to corrected divergence)
        pres_correction += ξ * W[k] * div_tend[lm, k]
    end
    
    pres_tend[lm] += pres_correction
end