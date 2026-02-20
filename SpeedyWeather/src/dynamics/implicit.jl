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
    return implicit.time_step = implicit.α * dt  # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog) from input dt
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
        model::ShallowWater
    )

    (; div_tend, pres_tend) = diagn.tendencies  # tendency of divergence and pressure/η
    div_old, div_new = get_steps(progn.div)   # divergence at t, t+dt
    pres_old, pres_new = get_steps(progn.pres)  # pressure/η at t, t+dt

    # unpack with [] as stored in a RefValue for mutation during initialization
    H = model.atmosphere.layer_thickness        # layer thickness [m], undisturbed, no mountains
    g = model.planet.gravity                    # gravitational acceleration [m/s²]
    ξ = implicit.time_step                      # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog)

    # Get precomputed l_indices from the spectrum
    l_indices = div_tend.spectrum.l_indices

    # GPU kernel launch
    arch = architecture(div_tend)
    launch!(
        arch, SpectralWorkOrder, size(div_tend), implicit_shallow_water_kernel!,
        div_tend, pres_tend, div_old, div_new, pres_old, pres_new, l_indices, H, g, ξ
    )

    zero_last_degree!(div_tend)
    zero_last_degree!(pres_tend)
    return nothing
end

@kernel inbounds = true function implicit_shallow_water_kernel!(
        div_tend, pres_tend, div_old, div_new, pres_old, pres_new, l_indices,
        @Const(H), @Const(g), @Const(ξ)
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
    G_div = div_tend[lm, k] - g * ∇² * (pres_old[lm] - pres_new[lm])
    G_η = pres_tend[lm] - H * (div_old[lm, k] - div_new[lm, k])

    # Using the Gs correct the tendencies for semi-implicit time stepping
    S⁻¹ = inv(1 - ξ^2 * H * g * ∇²)  # operator to invert
    div_tend[lm, k] = S⁻¹ * (G_div - ξ * g * ∇² * G_η)
    pres_tend[lm] = G_η - ξ * H * div_tend[lm, k]
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
        IntType,
    } <: AbstractImplicit

    # DIMENSIONS
    "Spectral resolution"
    trunc::IntType

    "Number of vertical layers"
    nlayers::IntType

    # PARAMETERS
    "Time-step coefficient: 0=explicit, 0.5=centred implicit, 1=backward implicit"
    α::NF = 1

    "Reinitialize at restart when initialized=true"
    reinitialize::Bool = true

    "Flag automatically set to true when initialize! has been called"
    initialized::Bool = false

    # PRECOMPUTED ARRAYS, to be initialized with initialize!
    "vertical temperature profile, obtained from diagn on first time step"
    temp_profile::VectorType = zeros(NF, nlayers)

    "time step 2α*Δt packed in RefValue for mutability"
    ξ::Base.RefValue{NF} = Ref{NF}(0)

    "divergence: operator for the geopotential calculation"
    R::MatrixType = zeros(NF, nlayers, nlayers)

    "divergence: the -RdTₖ∇² term excl the eigenvalues from ∇² for divergence"
    U::VectorType = zeros(NF, nlayers)

    "temperature: operator for the TₖD + κTₖD(ln pₛ)/Dt term"
    L::MatrixType = zeros(NF, nlayers, nlayers)

    "pressure: vertical averaging of the -D̄ term in the log surface pres equation"
    W::VectorType = zeros(NF, nlayers)

    "components to construct L, 1/ 2Δσ"
    L0::VectorType = zeros(NF, nlayers)

    "vert advection term in the temperature equation (below+above)"
    L1::MatrixType = zeros(NF, nlayers, nlayers)

    "factor in front of the `div_sum_above` term"
    L2::VectorType = zeros(NF, nlayers)

    "`_sum_above` operator itself"
    L3::MatrixType = zeros(NF, nlayers, nlayers)

    "factor in front of div term in Dlnpₛ/Dt"
    L4::VectorType = zeros(NF, nlayers)

    "for every l the matrix to be inverted"
    S::MatrixType = zeros(NF, nlayers, nlayers)

    "combined inverted operator: S = 1 - ξ²(RL + UW)"
    S⁻¹::TensorType = zeros(NF, trunc + 2, nlayers, nlayers)
end

"""$(TYPEDSIGNATURES)
Generator using the resolution from SpectralGrid."""
function ImplicitPrimitiveEquation(spectral_grid::SpectralGrid; kwargs...)
    (; NF, VectorType, MatrixType, TensorType, trunc, nlayers) = spectral_grid
    return ImplicitPrimitiveEquation{NF, VectorType, MatrixType, TensorType, typeof(trunc)}(;
        trunc, nlayers, kwargs...
    )
end

# function barrier to unpack the constants struct for primitive eq models
function initialize!(
        I::ImplicitPrimitiveEquation,
        dt::Real,
        diagn::DiagnosticVariables,
        model::PrimitiveEquation,
    )
    model.dynamics || return nothing    # escape immediately if no dynamics
    (; geometry, geopotential, atmosphere, adiabatic_conversion) = model
    return initialize!(I, dt, diagn, geometry, geopotential, atmosphere, adiabatic_conversion)
end

"""$(TYPEDSIGNATURES)
Initialize the implicit terms for the PrimitiveEquation models."""
function initialize!(
        implicit::ImplicitPrimitiveEquation,
        dt::Real,                                           # the scaled time step radius*dt
        diagn::DiagnosticVariables,
        geometry::AbstractGeometry,
        geopotential::AbstractGeopotential,
        atmosphere::AbstractAtmosphere,
        adiabatic_conversion::AbstractAdiabaticConversion,
    )
    # option to skip reinitialization at restart
    (implicit.initialized && !implicit.reinitialize) && return nothing

    (; trunc, nlayers, α) = implicit
    (; σ_levels_full, σ_levels_thick) = geometry
    (; R_dry, κ) = atmosphere
    (; Δp_geopot_half, Δp_geopot_full) = geopotential
    (; σ_lnp_A, σ_lnp_B) = adiabatic_conversion

    arch = architecture(implicit.S⁻¹)

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

    (; temp_profile, S⁻¹, L, R, U, W, L0, L1, L2, L3, L4) = implicit

    # use current vertical temperature profile
    temp_profile .= diagn.temp_average

    # return immediately if temp_profile contains NaRs, model blew up in that case
    all(isfinite.(temp_profile)) || return nothing

    ξ = α * dt                        # dt = 2Δt for leapfrog, but = Δt, Δ/2 in first_timesteps!
    implicit.ξ[] = ξ                # also store in Implicit struct

    # index vectors for broadcasting: rows = 1:nlayers (column), cols = (1:nlayers)' (row)
    rows = (1:nlayers)
    cols = (1:nlayers)'

    # DIVERGENCE OPERATORS (called g in Hoskins and Simmons 1975, eq 11 and Appendix 1)
    # R[row, k] = -Δp_geopot_full[k] for row <= k, additionally -Δp_geopot_half[k] for row < k
    R .= .-Δp_geopot_full[cols] .* (rows .<= cols) .- Δp_geopot_half[cols] .* (rows .< cols)

    # U = -R_dry * temp_profile (the R_d*Tₖ∇² term excl eigenvalues from ∇² for divergence)
    U .= .-R_dry .* temp_profile

    # TEMPERATURE OPERATOR (called τ in Hoskins and Simmons 1975, eq 9 and Appendix 1)
    L0 .= inv.(2 .* σ_levels_thick)
    L2 .= κ .* temp_profile .* σ_lnp_A
    L4 .= κ .* temp_profile .* σ_lnp_B

    # L1[k, r]: vertical advection operator
    # k indices for below/above neighbours, clamped to 1:nlayers
    k_above = max.(1, rows .- 1)
    k_below = min.(rows .+ 1, nlayers)
    ΔT_below = temp_profile[k_below] .- temp_profile[rows]
    ΔT_above = temp_profile[rows] .- temp_profile[k_above]
    σₖ = σ_levels_full[rows]
    σₖ_above = σ_levels_full[k_above]
    Δσᵣ = σ_levels_thick[cols]

    L1 .= ΔT_below .* Δσᵣ .* σₖ .- Δσᵣ .* (rows .>= cols) .+
           ΔT_above .* Δσᵣ .* σₖ_above .- Δσᵣ .* ((rows .- 1) .>= cols)

    # L3[r, k]: sum_above operator — lower triangle gets σ_levels_thick[k]
    L3 .= σ_levels_thick[cols] .* (rows .> cols)

    # Combine all operators into L = Diagonal(L0)*L1 + Diagonal(L2)*L3 + Diagonal(L4)
    L .= L0 .* L1 .+ L2 .* L3 .+ L4 .* (rows .== cols)

    # PRESSURE OPERATOR (called πᵣ in Hoskins and Simmons, 1975 Appendix 1)
    W .= .-σ_levels_thick

    # solving the equations above for δD yields
    # δD = SG, with G = G_D + ξRG_T + ξUG_lnps and the operator S
    # S = 1 - ξ²(RL + UW) that has to be inverted to obtain δD from the Gs
    # Compute S⁻¹ for every l via Gauss-Jordan elimination
    S_scratch = similar(S⁻¹)
    launch!(arch, LinearWorkOrder, (trunc + 1,), _implicit_invert_S_kernel!,
        S⁻¹, S_scratch, R, L, U, W, ξ, nlayers)

    return nothing
end

# Compute S⁻¹ for each spectral degree l via Gauss-Jordan elimination
# S = I - ξ²*eigenvalue*(R*L + U*Wᵀ), then invert S per l
# S_scratch[l,:,:] stores the S matrix during elimination,
# while S⁻¹[l,:,:] tracks the identity → inverse transformation
@kernel inbounds = true function _implicit_invert_S_kernel!(
        S⁻¹,                           # Output: (trunc+2) × nlayers × nlayers tensor
        S_scratch,                     # Scratch: same shape as S⁻¹, stores S per l
        @Const(R),                      # Input: divergence operator matrix
        @Const(L),                      # Input: temperature operator matrix
        @Const(U),                      # Input: divergence vector
        @Const(W),                      # Input: pressure vector
        @Const(ξ),                      # Input: semi-implicit time step coefficient
        @Const(nlayers),                # Input: number of layers
    )
    l = @index(Global, Linear)      # spectral degree (1-based)

    NF = eltype(S⁻¹)
    eigenvalue = -l * (l - 1)       # 1-based, -l*(l+1) → -l*(l-1)
    ξ²λ = ξ^2 * eigenvalue

    # Compute S = I - ξ²λ*(R*L + U*Wᵀ) into S_scratch[l,:,:]
    # and initialize S⁻¹[l,:,:] as identity
    for k in 1:nlayers
        for r in 1:nlayers
            RL_kr = zero(NF)
            for j in 1:nlayers
                RL_kr += R[k, j] * L[j, r]
            end
            S_scratch[l, k, r] = (k == r ? one(NF) : zero(NF)) - ξ²λ * (RL_kr + U[k] * W[r])
            S⁻¹[l, k, r] = k == r ? one(NF) : zero(NF)
        end
    end

    # Gauss-Jordan elimination: reduce S_scratch[l,:,:] (S) to I,
    # applying the same row operations to S⁻¹[l,:,:] (I → S⁻¹)
    for pivot in 1:nlayers
        inv_pivot = inv(S_scratch[l, pivot, pivot])

        # Scale pivot row
        for r in 1:nlayers
            S_scratch[l, pivot, r] *= inv_pivot
            S⁻¹[l, pivot, r] *= inv_pivot
        end

        # Eliminate all other rows
        for k in 1:nlayers
            if k != pivot
                factor = S_scratch[l, k, pivot]
                for r in 1:nlayers
                    S_scratch[l, k, r] -= factor * S_scratch[l, pivot, r]
                    S⁻¹[l, k, r] -= factor * S⁻¹[l, pivot, r]
                end
            end
        end
    end
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
    G = diagn.dynamics.a              # reuse work arrays, used for combined tendency G
    geopotential = diagn.dynamics.b   # used for geopotential
    geopotential .= 0

    # Get precomputed l_indices from the spectrum
    l_indices = temp_tend.spectrum.l_indices

    arch = architecture(temp_tend)

    # Single kernel: All implicit correction steps for each spectral mode
    lm_size = size(pres_tend, 1)
    launch!(
        arch, LinearWorkOrder, (lm_size,),
        implicit_primitive_single_kernel!,
        temp_tend, pres_tend, div_tend, G, geopotential,
        div_old, div_new, S⁻¹, R, U, L, W, l_indices, ξ, nlayers
    )

    zero_last_degree!(div_tend)
    zero_last_degree!(pres_tend)
    zero_last_degree!(temp_tend)

    pres_tend.data[1:1] .= 0    # mass conservation

    return nothing
end

# Single kernel that does all steps for one spectral mode
@kernel inbounds = true function implicit_primitive_single_kernel!(
        temp_tend, pres_tend, div_tend, G, geopotential,
        div_old, div_new, @Const(S⁻¹), @Const(R), @Const(U), @Const(L), @Const(W), @Const(l_indices),
        @Const(ξ), @Const(nlayers)
    )
    lm = @index(Global, Linear)

    # Get degree l for this spectral mode
    l = l_indices[lm]

    # Step 1: Move implicit terms of temperature equation from time step i to i-1
    # RHS_expl(Vⁱ) + RHS_impl(Vⁱ⁻¹) = RHS(Vⁱ) + RHS_impl(Vⁱ⁻¹ - Vⁱ)
    for k in 1:nlayers
        temp_tend_val = zero(eltype(temp_tend))
        for r in 1:nlayers
            temp_tend_val += L[k, r] * (div_old[lm, r] - div_new[lm, r])
        end
        temp_tend[lm, k] += temp_tend_val
    end

    for k in 1:nlayers
        # skip 1:k-1 as integration is surface to k
        geopotential_val = zero(eltype(geopotential))
        for r in k:nlayers
            geopotential_val += R[k, r] * temp_tend[lm, r]
        end
        geopotential[lm, k] = geopotential_val
    end

    eigenvalue = -l * (l - 1)  # 1-based, -l*(l+1) → -l*(l-1)

    # Step 2: Calculate the ξ*R*G_T term, vertical integration of geopotential
    # (excl ξ, this is done in step 3)

    # Step 3: Calculate the G = G_D + ξRG_T + ξUG_lnps terms
    # ∇² not part of U so *eigenvalues here
    for k in 1:nlayers
        G[lm, k] = div_tend[lm, k] + ξ * eigenvalue * (U[k] * pres_tend[lm] + geopotential[lm, k])
    end

    # Step 4: Now solve δD = S⁻¹G to correct divergence tendency
    for k in 1:nlayers
        div_val = zero(eltype(div_tend))
        for r in 1:nlayers
            div_val += S⁻¹[l, k, r] * G[lm, r]
        end
        div_tend[lm, k] = div_val
    end

    # Step 5: Semi implicit corrections for temperature and pressure

    # Step 5a: Temperature correction δT = G_T + ξLδD
    for k in 1:nlayers
        temp_correction = zero(eltype(temp_tend))
        for r in 1:nlayers
            temp_correction += ξ * L[k, r] * div_tend[lm, r]
        end
        temp_tend[lm, k] += temp_correction
    end

    # Step 5b: Pressure correction δlnpₛ = G_lnpₛ + ξWδD
    pres_correction = zero(eltype(pres_tend))
    for k in 1:nlayers
        pres_correction += ξ * W[k] * div_tend[lm, k]
    end
    pres_tend[lm] += pres_correction
end
