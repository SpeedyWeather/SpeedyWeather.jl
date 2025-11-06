abstract type AbstractImplicit <: AbstractModelComponent end

# model.implicit=nothing (for BarotropicModel)
initialize!(::Nothing, dt::Real, ::DiagnosticVariables, ::AbstractModel) = nothing
implicit_correction!(::DiagnosticVariables, ::PrognosticVariables, ::Nothing, ::AbstractModel) = nothing

# SHALLOW WATER MODEL
export ImplicitShallowWater

"""Struct that holds various precomputed arrays for the semi-implicit correction to
prevent gravity waves from amplifying in the shallow water model. The implicit time step
between i-1 and i+1 is controlled by the parameter `α` in the range [0, 1]
    
    α = 0   means the gravity wave terms are evaluated at i-1 (forward)
    α = 0.5 evaluates at i+1 and i-1 (centered implicit)
    α = 1   evaluates at i+1 (backward implicit)
    α ∈ [0.5, 1] are also possible which controls the strength of the gravity wave dampening.
    α = 0.5 slows gravity waves and prevents them from amplifying
    α > 0.5 will dampen the gravity waves within days to a few timesteps (α=1)

Fields are
$(TYPEDFIELDS)"""
@kwdef mutable struct ImplicitShallowWater{NF} <: AbstractImplicit
    "[OPTION] coefficient for semi-implicit computations to filter gravity waves, 0.5 <= α <= 1"
    α::NF = 1
    
    "Time step [s], = αdt = 2αΔt (for leapfrog)"
    time_step::NF = 0
end

"""
$(TYPEDSIGNATURES)
Generator using the resolution from `spectral_grid`."""
ImplicitShallowWater(SG::SpectralGrid; kwargs...) = ImplicitShallowWater{SG.NF}(; kwargs...)

"""
$(TYPEDSIGNATURES)
Update the implicit terms in `implicit` for the shallow water model as they depend on the time step `dt`."""
function initialize!(implicit::ImplicitShallowWater, dt::Real, args...)
    implicit.time_step = implicit.α*dt  # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog) from input dt
end

# implicit shallow water has no precomputed arrays, so implicit.initialized is not defined
set_initialized!(implicit::ImplicitShallowWater) = nothing
set_initialized!(implicit::Nothing) = nothing

"""
$(TYPEDSIGNATURES)
Apply correction to the tendencies in `diagn` to prevent the gravity waves from amplifying.
The correction is implicitly evaluated using the parameter `implicit.α` to switch between
forward, centered implicit or backward evaluation of the gravity wave terms."""
function implicit_correction!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    implicit::ImplicitShallowWater,
    model::ShallowWater)

    (; div_tend, pres_tend) = diagn.tendencies  # tendency of divergence and pressure/η
    div_old, div_new   = get_steps(progn.div)   # divergence at t, t+dt
    pres_old, pres_new = get_steps(progn.pres)  # pressure/η at t, t+dt

    # unpack with [] as stored in a RefValue for mutation during initialization
    H = model.atmosphere.layer_thickness        # layer thickness [m], undisturbed, no mountains
    g = model.planet.gravity                    # gravitational acceleration [m/s²]
    ξ = implicit.time_step                      # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog)

    # Get precomputed l_indices from the spectrum
    l_indices = div_tend.spectrum.l_indices

    # GPU kernel launch
    arch = architecture(div_tend)
    launch!(arch, SpectralWorkOrder, size(div_tend), implicit_shallow_water_kernel!,
            div_tend, pres_tend, div_old, div_new, pres_old, pres_new, l_indices, H, g, ξ)
end

@kernel inbounds=true function implicit_shallow_water_kernel!(
    div_tend, pres_tend, div_old, div_new, pres_old, pres_new, l_indices,
    @Const(H), @Const(g), @Const(ξ)
)
    I = @index(Global, Cartesian)
    lm = I[1]  # single index lm corresponding to harmonic l, m
    k = I[2]   # layer index
    
    # Use precomputed l index from spectrum
    l = l_indices[lm]
    
    # eigenvalue, with without 1/radius², 1-based -l*(l+1) → -l*(l-1)
    ∇² = -l*(l-1)
    
    # Calculate the G = N(Vⁱ) + NI(Vⁱ⁻¹ - Vⁱ) term.
    # Vⁱ is a prognostic variable at time step i
    # N is the right hand side of ∂V\∂t = N(V)
    # NI is the part of N that's calculated semi-implicitily: N = NE + NI
    G_div = div_tend[lm, k] - g*∇²*(pres_old[lm] - pres_new[lm])
    G_η   = pres_tend[lm] - H*(div_old[lm, k] - div_new[lm, k])
    
    # Using the Gs correct the tendencies for semi-implicit time stepping
    S⁻¹ = inv(1 - ξ^2*H*g*∇²)  # operator to invert
    div_tend[lm, k] = S⁻¹*(G_div - ξ*g*∇²*G_η)
    pres_tend[lm] = G_η - ξ*H*div_tend[lm, k]
end

# CPU version (original implementation)
function implicit_correction_cpu!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    implicit::ImplicitShallowWater,
    model::ShallowWater)

    (; div_tend, pres_tend) = diagn.tendencies  # tendency of divergence and pressure/η
    div_old, div_new   = get_steps(progn.div)   # divergence at t, t+dt
    pres_old, pres_new = get_steps(progn.pres)  # pressure/η at t, t+dt

    # unpack with [] as stored in a RefValue for mutation during initialization
    H = model.atmosphere.layer_thickness        # layer thickness [m], undisturbed, no mountains
    g = model.planet.gravity                    # gravitational acceleration [m/s²]
    ξ = implicit.time_step                      # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog)

    lmax, mmax = size(div_tend, OneBased, as=Matrix)

    @inbounds for k in eachmatrix(div_tend)
        lm = 0
        for m in 1:mmax
            for l in m:lmax
                lm += 1                     # single index lm corresponding to harmonic l, m with a LowerTriangularMatrix
                ∇² = -l*(l-1)               # eigenvalue, with without 1/radius², 1-based -l*(l+1) → -l*(l-1)

                # calculate the G = N(Vⁱ) + NI(Vⁱ⁻¹ - Vⁱ) term.
                # Vⁱ is a prognostic variable at time step i
                # N is the right hand side of ∂V\∂t = N(V)
                # NI is the part of N that's calculated semi-implicitily: N = NE + NI
                G_div = div_tend[lm, k] - g*∇²*(pres_old[lm] - pres_new[lm])
                G_η   = pres_tend[lm] - H*(div_old[lm, k] - div_new[lm, k])

                # using the Gs correct the tendencies for semi-implicit time stepping
                S⁻¹ = inv(1 - ξ^2*H*g*∇²)                   # operator to invert
                div_tend[lm, k] = S⁻¹*(G_div - ξ*g*∇²*G_η)
                pres_tend[lm] = G_η - ξ*H*div_tend[lm, k]
            end
        end
    end
end

export ImplicitPrimitiveEquation

"""
Struct that holds various precomputed arrays for the semi-implicit correction to
prevent gravity waves from amplifying in the primitive equation model.
$(TYPEDFIELDS)"""
@kwdef mutable struct ImplicitPrimitiveEquation{
    NF,             # number format
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
    "Time-step coefficient: 0=explicit, 0.5=centred implicit, 1=backward implicit"
    α::NF = 1

    "Reinitialize at restart when initialized=true"
    reinitialize::Bool = true

    "Flag automatically set to true when initialize! has been called"
    initialized::Bool = false

    # PRECOMPUTED ARRAYS, to be initiliazed with initialize!
    "vertical temperature profile, obtained from diagn on first time step"
    temp_profile::VectorType = zeros(NF, nlayers)

    "time step 2α*Δt packed in RefValue for mutability"
    ξ::Base.RefValue{NF} = Ref{NF}(0)

    "divergence: operator for the geopotential calculation"
    R::MatrixType = zeros(NF, nlayers, nlayers)

    "divergence: the -RdTₖ∇² term excl the eigenvalues from ∇² for divergence"
    U::VectorType = zeros(NF, nlayers)

    "temperature: operator for the TₖD + κTₖDlnps/Dt term"
    L::MatrixType = zeros(NF, nlayers, nlayers)

    "pressure: vertical averaging of the -D̄ term in the log surface pres equation"
    W::VectorType = zeros(NF, nlayers)

    "components to construct L, 1/ 2Δσ"
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

"""$(TYPEDSIGNATURES)
Generator using the resolution from SpectralGrid."""
function ImplicitPrimitiveEquation(spectral_grid::SpectralGrid; kwargs...)
    (; NF, VectorType, MatrixType, TensorType, trunc, nlayers) = spectral_grid
    return ImplicitPrimitiveEquation{NF, VectorType, MatrixType, TensorType}(;
        trunc, nlayers, kwargs...)
end

# function barrier to unpack the constants struct for primitive eq models
function initialize!(
    I::ImplicitPrimitiveEquation,
    dt::Real,
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)
    (; geometry, geopotential, atmosphere, adiabatic_conversion) = model
    initialize!(I, dt, diagn, geometry, geopotential, atmosphere, adiabatic_conversion)
end

"""$(TYPEDSIGNATURES)
Initialize the implicit terms for the PrimitiveEquation models."""
function initialize!(
    implicit::ImplicitPrimitiveEquation,
    dt::Real,                                           # the scaled time step radius*dt
    diagn::DiagnosticVariables{NF},
    geometry::AbstractGeometry,
    geopotential::AbstractGeopotential,
    atmosphere::AbstractAtmosphere,
    adiabatic_conversion::AbstractAdiabaticConversion,
) where NF

    # option to skip reinitialization at restart
    (implicit.initialized && !implicit.reinitialize) && return nothing

    (; trunc, nlayers, α) = implicit
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

    # return immediately if temp_profile contains NaRs, model blew up in that case
    all(isfinite.(temp_profile_cpu)) || return nothing

    # set up R, U, L, W operators from
    # δD = G_D + ξ(RδT + Uδlnps)        divergence D correction
    # δT = G_T + ξLδD                   temperature T correction
    # δlnps = G_lnps + ξWδD             log surface pressure lnps correction
    #
    # G_X is the uncorrected explicit tendency calculated as RHS_expl(Xⁱ) + RHS_impl(Xⁱ⁻¹)
    # with RHS_expl being the nonlinear terms calculated from the centered time step i
    # and RHS_impl are the linear terms that are supposed to be calcualted semi-implicitly
    # however, they have sofar only been evaluated explicitly at time step i-1
    # and are subject to be corrected to δX following the equations above
    # R, U, L, W are linear operators that are therefore defined here and inverted
    # to obtain δD first, and then δT and δlnps through substitution

    ξ = α*dt                        # dt = 2Δt for leapfrog, but = Δt, Δ/2 in first_timesteps!
    implicit.ξ[] = ξ                # also store in Implicit struct

    # DIVERGENCE OPERATORS (called g in Hoskins and Simmons 1975, eq 11 and Appendix 1)
    @inbounds for k in 1:nlayers                # vertical geopotential integration as matrix operator
        R_cpu[1:k, k] .= -Δp_geopot_full_cpu[k]         # otherwise equivalent to geopotential! with zero orography
        R_cpu[1:k-1, k] .+= -Δp_geopot_half_cpu[k]      # incl the minus but excluding the eigenvalues as with U
    end
    U_cpu .= -R_dry*temp_profile_cpu        # the R_d*Tₖ∇² term excl the eigenvalues from ∇² for divergence

    # TEMPERATURE OPERATOR (called τ in Hoskins and Simmons 1975, eq 9 and Appendix 1)
    L0_cpu .= 1 ./ 2σ_levels_thick_cpu
    L2_cpu .= κ*temp_profile_cpu.*σ_lnp_A_cpu    # factor in front of the div_sum_above term
    L4_cpu .= κ*temp_profile_cpu.*σ_lnp_B_cpu    # factor in front of div term in Dlnps/Dt

    @inbounds for k in 1:nlayers
        Tₖ = temp_profile_cpu[k]                    # average temperature at k
        k_above = max(1, k-1)                       # layer index above
        k_below = min(k+1, nlayers)                 # layer index below
        ΔT_above = Tₖ - temp_profile_cpu[k_above]   # temperature difference to layer above
        ΔT_below = temp_profile_cpu[k_below] - Tₖ   # and to layer below
        σₖ = σ_levels_full_cpu[k]                   # should be Σ_r=1^k Δσᵣ for model top at >0hPa
        σₖ_above = σ_levels_full_cpu[k_above]

        for r in 1:nlayers
            L1_cpu[k, r] = ΔT_below*σ_levels_thick_cpu[r]*σₖ         # vert advection operator below
            L1_cpu[k, r] -= k>=r ? σ_levels_thick_cpu[r] : zero(NF)

            L1_cpu[k, r] += ΔT_above*σ_levels_thick_cpu[r]*σₖ_above   # vert advection operator above
            L1_cpu[k, r] -= (k-1)>=r ? σ_levels_thick_cpu[r] : zero(NF)
        end

        # _sum_above operator itself
        L3_cpu[1:k, k] .= 0                              # fill upper triangle + diagonal with zeros
        L3_cpu[k+1:end, k] .= σ_levels_thick_cpu[k]      # vert integration top to k-1
    end

    L_cpu .= Diagonal(L0_cpu)*L1_cpu .+ Diagonal(L2_cpu)*L3_cpu .+ Diagonal(L4_cpu)  # combine all operators into L

    # PRESSURE OPERATOR (called πᵣ in Hoskins and Simmons, 1975 Appendix 1)
    W_cpu .= -σ_levels_thick_cpu                # the -D̄ term in the log surface pres equation

    # solving the equations above for δD yields
    # δD = SG, with G = G_D + ξRG_T + ξUG_lnps and the operator S
    # S = 1 - ξ²(RL + UW) that has to be inverted to obtain δD from the Gs
    I = LinearAlgebra.I(nlayers)
    @inbounds for l in 1:trunc+1
        eigenvalue = -l*(l-1)           # 1-based, -l*(l+1) → -l*(l-1)
        S_cpu .= I .- ξ^2*eigenvalue*(R_cpu*L_cpu .+ U_cpu*W_cpu')

        # inv(S) but saving memory:
        luS = LinearAlgebra.lu!(S_cpu)      # in-place LU decomposition (overwriting S)
        Sinv = L1_cpu                       # reuse L1 matrix to store inv(S)
        Sinv .= I                           # use ldiv! so last arg needs to be unity matrix
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
Apply the implicit corrections to dampen gravity waves in the primitive equation models."""
function implicit_correction!(
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
    div_old, div_new = get_steps(progn.div)
    G = diagn.dynamics.a        # reuse work arrays, used for combined tendency G
    geopot = diagn.dynamics.b   # used for geopotential
    geopot .= 0

    # Get precomputed l_indices from the spectrum
    l_indices = temp_tend.spectrum.l_indices

    arch = architecture(temp_tend)
    
    # MOVE THE IMPLICIT TERMS OF THE TEMPERATURE EQUATION FROM TIME STEP i TO i-1
    # geopotential and linear pressure gradient (divergence equation) are already evaluated at i-1
    # so is the -D̄ term for surface pressure in tendencies!

    # SEMI IMPLICIT CORRECTIONS FOR DIVERGENCE
    # calculate the combined tendency G = G_D + ξRG_T + ξUG_lnps to solve for divergence δD

    # Kernel 1: Temperature update and geopotential calculation
    launch!(arch, SpectralWorkOrder, size(temp_tend), 
            implicit_temp_geopot_kernel!, 
            temp_tend, geopot, div_old, div_new, R, L, nlayers)
    
    # Kernel 2: Calculate G = G_D + ξRG_T + ξUG_lnps
    launch!(arch, SpectralWorkOrder, size(div_tend),
            implicit_calculate_G_kernel!,
            G, div_tend, geopot, pres_tend, U, l_indices, ξ, nlayers, trunc)
    
    # Kernel 3: Solve for divergence δD = S⁻¹G
    launch!(arch, SpectralWorkOrder, size(div_tend),
            implicit_solve_divergence_kernel!,
            div_tend, G, S⁻¹, l_indices, nlayers, trunc)
    
    # Kernel 4: Correct temperature
    launch!(arch, SpectralWorkOrder, size(temp_tend),
            implicit_correct_temp_kernel!,
            temp_tend, div_tend, L, l_indices, ξ, nlayers, trunc)
    
    # Kernel 5: Correct pressure (parallelize only over lm, not k)
    lm_size = size(pres_tend, 1)
    launch!(arch, SpectralWorkOrder, (lm_size,),
            implicit_correct_pres_kernel!,
            pres_tend, div_tend, W, l_indices, ξ, nlayers, trunc)

    zero_last_degree!(div_tend)
    zero_last_degree!(pres_tend)
    zero_last_degree!(temp_tend)

    return 
end

# Alternative single-kernel version (parallelize only over lm, not k)
# This version is simpler and avoids all race conditions by doing all work per spectral mode
function implicit_correction_lm_only!(
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
    div_old, div_new = get_steps(progn.div)
    G = diagn.dynamics.a        # reuse work arrays, used for combined tendency G
    geopot = diagn.dynamics.b   # used for geopotential
    geopot .= 0

    # Get precomputed l_indices from the spectrum
    l_indices = temp_tend.spectrum.l_indices

    arch = architecture(temp_tend)
    
    # Single kernel: All implicit correction steps for each spectral mode
    lm_size = size(pres_tend, 1)
    launch!(arch, LinearWorkOrder, (lm_size,),
            implicit_primitive_single_kernel!,
            temp_tend, pres_tend, div_tend, G, geopot,
            div_old, div_new, S⁻¹, R, U, L, W, l_indices, ξ, nlayers)

    zero_last_degree!(div_tend)
    zero_last_degree!(pres_tend)
    zero_last_degree!(temp_tend)

    return geopot
end

# Single kernel that does all steps for one spectral mode
@kernel inbounds=true function implicit_primitive_single_kernel!(
    temp_tend, pres_tend, div_tend, G, geopot,
    div_old, div_new, @Const(S⁻¹), @Const(R), @Const(U), @Const(L), @Const(W), @Const(l_indices),
    @Const(ξ), @Const(nlayers)
)
    lm = @index(Global, Linear)
    
    # Get degree l for this spectral mode
    l = l_indices[lm]
    
    # Step 1: Move implicit terms of temperature equation from time step i to i-1
    # RHS_expl(Vⁱ) + RHS_impl(Vⁱ⁻¹) = RHS(Vⁱ) + RHS_impl(Vⁱ⁻¹ - Vⁱ)
    for k in 1:nlayers
        for r in 1:nlayers
            temp_tend[lm, k] += L[k, r] * (div_old[lm, r] - div_new[lm, r])
        end
    end

    for k in 1:nlayers
        # skip 1:k-1 as integration is surface to k
        for r in k:nlayers
            geopot[lm, k] += R[k, r] * temp_tend[lm, r]
        end
    end
            
    eigenvalue = -l*(l-1)  # 1-based, -l*(l+1) → -l*(l-1)
        
    # Step 2: Calculate the ξ*R*G_T term, vertical integration of geopotential
    # (excl ξ, this is done in step 3)
        
    # Step 3: Calculate the G = G_D + ξRG_T + ξUG_lnps terms
    # ∇² not part of U so *eigenvalues here
    for k in 1:nlayers
        G[lm, k] = div_tend[lm, k] + ξ*eigenvalue*(U[k]*pres_tend[lm] + geopot[lm, k])
    end
        
    # Step 4: Now solve δD = S⁻¹G to correct divergence tendency
    # Zero div_tend first so it can be used as an accumulator
    for k in 1:nlayers
        div_tend[lm, k] = 0
        for r in 1:nlayers
            div_tend[lm, k] += S⁻¹[l, k, r] * G[lm, r]
        end
    end
        
    # Step 5: Semi implicit corrections for temperature and pressure
        
    # Step 5a: Temperature correction δT = G_T + ξLδD
    for k in 1:nlayers
        for r in 1:nlayers
            temp_tend[lm, k] += ξ * L[k, r] * div_tend[lm, r]
        end
    end
        
    # Step 5b: Pressure correction δlnpₛ = G_lnpₛ + ξWδD
    for k in 1:nlayers
        pres_tend[lm] += ξ * W[k] * div_tend[lm, k]
    end
end

# Kernel 1: Combined temperature update and geopotential calculation
@kernel inbounds=true function implicit_temp_geopot_kernel!(
    temp_tend, geopot, div_old, div_new, R, L, @Const(nlayers)
)
    I = @index(Global, Cartesian)
    lm = I[1]
    k = I[2]
    
    # Move implicit terms of temperature equation from time step i to i-1
    # RHS_expl(Vⁱ) + RHS_impl(Vⁱ⁻¹) = RHS(Vⁱ) + RHS_impl(Vⁱ⁻¹ - Vⁱ)
    for r in 1:nlayers
        temp_tend[I] += L[k, r] * (div_old[lm, r] - div_new[lm, r])
    end
    
    # Calculate the ξ*R*G_T term, vertical integration of geopotential
    # (excl ξ, this is done in kernel 2)
    # skip 1:k-1 as integration is surface to k
    sum_val = zero(eltype(geopot))
    for r in k:nlayers
        sum_val += R[k, r] * temp_tend[lm, r]
    end
    geopot[I] = sum_val
end


# Kernel 2: Calculate G (must complete before kernel 3)
@kernel inbounds=true function implicit_calculate_G_kernel!(
    G, div_tend, geopot, pres_tend, U, l_indices,
    @Const(ξ), @Const(nlayers), @Const(trunc)
)
    I = @index(Global, Cartesian)
    lm = I[1]
    k = I[2]
    
    # Use precomputed l index from spectrum
    l = l_indices[lm]
    
    # Skip out-of-bounds modes
    eigenvalue = -l*(l-1)  # 1-based, -l*(l+1) → -l*(l-1)
    
    # Calculate the G = G_D + ξRG_T + ξUG_lnps terms using geopot from kernel 1
    # ∇² not part of U so *eigenvalues here
    G[lm, k] = div_tend[lm, k] + ξ*eigenvalue*(U[k]*pres_tend[lm] + geopot[lm, k])
end

# Kernel 3: Solve for divergence δD = S⁻¹G
@kernel inbounds=true function implicit_solve_divergence_kernel!(
    div_tend, G, S⁻¹, l_indices,
    @Const(nlayers), @Const(trunc)
)
    I = @index(Global, Cartesian)
    lm = I[1]
    k = I[2]
    
    # Use precomputed l index from spectrum
    l = l_indices[lm]
    
    # Skip out-of-bounds modes
    # Now solve δD = S⁻¹G to correct divergence tendency
    sum_val = zero(eltype(div_tend))
    for r in 1:nlayers
        sum_val += S⁻¹[l, k, r] * G[lm, r]
    end
    div_tend[lm, k] = sum_val   
end

# Kernel 4: Correct temperature
@kernel inbounds=true function implicit_correct_temp_kernel!(
    temp_tend, div_tend, L, l_indices, @Const(ξ), @Const(nlayers), @Const(trunc)
)
    I = @index(Global, Cartesian)
    lm = I[1]
    k = I[2]
    
    # Get degree l
    l = l_indices[lm]
 
    # Semi implicit correction for temperature
    # δT = G_T + ξLδD
    for r in 1:nlayers
        temp_tend[lm, k] += ξ * L[k, r] * div_tend[lm, r]
    end
end

# Kernel 5: Correct pressure
@kernel inbounds=true function implicit_correct_pres_kernel!(
    pres_tend, div_tend, W, l_indices, @Const(ξ), @Const(nlayers), @Const(trunc)
)
    lm = @index(Global, Linear)
    
    # Get degree l
    l = l_indices[lm]
    
    # Skip out-of-bounds modes
    # Semi implicit correction for pressure
    # δlnpₛ = G_lnpₛ + ξWδD
    # Accumulate contributions from all layers for this spectral mode
    for k in 1:nlayers
        pres_tend[lm] += ξ * W[k] * div_tend[lm, k]
    end
end

# CPU version (original implementation)
function implicit_correction_cpu!(
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

    # MOVE THE IMPLICIT TERMS OF THE TEMPERATURE EQUATION FROM TIME STEP i TO i-1
    # geopotential and linear pressure gradient (divergence equation) are already evaluated at i-1
    # so is the -D̄ term for surface pressure in tendencies!
    (; temp_tend) = diagn.tendencies
    div_old, div_new = get_steps(progn.div)     # divergence at i-1 (old), i (new, i.e. current)

    for k in eachmatrix(temp_tend, div_old, div_new)
        for r in eachmatrix(temp_tend, div_old, div_new)
            for lm in eachharmonic(temp_tend, div_old, div_new)
                # RHS_expl(Vⁱ) + RHS_impl(Vⁱ⁻¹) = RHS(Vⁱ) + RHS_impl(Vⁱ⁻¹ - Vⁱ)
                # for temperature tendency do the latter as its cheaper.
                temp_tend[lm, k] += L[k, r] * (div_old[lm, r] - div_new[lm, r])
                # temp_tend[lm, k] += L[k, r] * div_old[lm, r]    # for the former
            end
        end
    end

    # SEMI IMPLICIT CORRECTIONS FOR DIVERGENCE
    # calculate the combined tendency G = G_D + ξRG_T + ξUG_lnps to solve for divergence δD
    (; pres_tend, div_tend) = diagn.tendencies
    G = diagn.dynamics.a        # reuse work arrays, used for combined tendency G
    geopot = diagn.dynamics.b   # used for geopotential
    geopot .= 0

    for k in 1:nlayers
        for r in k:nlayers      # skip 1:k-1 as integration is surface to k
            for lm in eachharmonic(temp_tend, div_old, div_new)
                # 1. the ξ*R*G_T term, vertical integration of geopotential (excl ξ, this is done in 2.)
                geopot[lm, k] += R[k, r]*temp_tend[lm, r]
            end
        end
    end 

    
    for k in 1:nlayers
        # 2. the G = G_D + ξRG_T + ξUG_lnps terms using geopot from above
        lm = 0
        for m in 1:trunc+1              # loops over all columns/order m
            for l in m:trunc+1          # but skips the lmax+2 degree (1-based)
                lm += 1                 # single index lm corresponding to harmonic l, m
                                        # ∇² not part of U so *eigenvalues here
                eigenvalue = -l*(l-1)   # 1-based, -l*(l+1) → -l*(l-1)
                G[lm, k] = div_tend[lm, k] + ξ*eigenvalue*(U[k]*pres_tend[lm] + geopot[lm, k])

                # div_tend is now in G, fill with zeros here so that it can be used as an accumulator
                # in the δD = S⁻¹G calculation below
                div_tend[lm, k] = 0
            end
            lm += 1         # skip last row, LowerTriangularArrays are of size lmax+2 x mmax+1
        end
    end

    # NOW SOLVE THE δD = S⁻¹G to correct divergence tendency
    for k in eachmatrix(div_tend, G)
        for r in eachmatrix(div_tend, G)
            lm = 0
            for m in 1:trunc+1      # loops over all columns/order m
                for l in m:trunc+1  # but skips the lmax+2 degree (1-based)
                    lm += 1         # single index lm corresponding to harmonic l, m
                    div_tend[lm, k] += S⁻¹[l, k, r]*G[lm, r]
                end
                lm += 1             # skip last row, LowerTriMatrices are of size lmax+2 x mmax+1
            end
        end
    end

    # SEMI IMPLICIT CORRECTIONS FOR PRESSURE AND TEMPERATURE, insert δD to get δT, δlnpₛ
    for k in eachmatrix(div_tend, temp_tend)
        for r in eachmatrix(div_tend, temp_tend)
            for lm in eachharmonic(div_tend, temp_tend)
                # δT = G_T + ξLδD
                #temp_tend[lm, k] += ξ*L[k, r]*div_tend[lm, r]
            end
        end

        for lm in eachharmonic(div_tend, temp_tend)
            # δlnpₛ = G_lnpₛ + ξWδD
            pres_tend[lm] += ξ*W[k]*div_tend[lm, k]
        end
    end

    zero_last_degree!(div_tend)
    zero_last_degree!(pres_tend)
    zero_last_degree!(temp_tend)

    return nothing
end