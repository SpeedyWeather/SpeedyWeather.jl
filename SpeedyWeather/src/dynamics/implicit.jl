abstract type AbstractImplicit <: AbstractModelComponent end

# model.implicit=nothing (for BarotropicModel)
initialize!(::Nothing, dt::Real, ::Variables, ::AbstractModel) = nothing
implicit_correction!(::Variables, ::Nothing, ::AbstractModel) = nothing

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
    implicit.time_step = implicit.α * dt  # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog) from input dt
    return implicit
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
        vars::Variables,
        implicit::ImplicitShallowWater,
        model::ShallowWater
    )

    div_tend = vars.tendencies.div                      # tendency of divergence and interface displacement η
    η_tend = vars.tendencies.η                          # tendency of divergence and interface displacement η
    div_old, div_new = get_steps(vars.prognostic.div)   # divergence at t, t+dt
    η_old, η_new = get_steps(vars.prognostic.η)         # η at t, t+dt

    H = model.atmosphere.layer_thickness        # layer thickness [m], undisturbed, no mountains
    g = model.planet.gravity                    # gravitational acceleration [m/s²]
    ξ = implicit.time_step                      # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog)

    # Get precomputed l_indices from the spectrum
    l_indices = div_tend.spectrum.l_indices

    # GPU kernel launch
    arch = architecture(div_tend)
    launch!(
        arch, SpectralWorkOrder, size(div_tend), implicit_shallow_water_kernel!,
        div_tend, η_tend, div_old, div_new, η_old, η_new, l_indices, H, g, ξ
    )

    zero_last_degree!(div_tend)
    zero_last_degree!(η_tend)
    return nothing
end

@kernel inbounds = true function implicit_shallow_water_kernel!(
        div_tend, η_tend, div_old, div_new, η_old, η_new, l_indices,
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
    G_div = div_tend[lm, k] - g * ∇² * (η_old[lm, k] - η_new[lm, k])
    G_η = η_tend[lm] - H * (div_old[lm, k] - div_new[lm, k])

    # Using the Gs correct the tendencies for semi-implicit time stepping
    S⁻¹ = inv(1 - ξ^2 * H * g * ∇²)  # operator to invert
    div_tend[lm, k] = S⁻¹ * (G_div - ξ * g * ∇² * G_η)
    η_tend[lm] = G_η - ξ * H * div_tend[lm, k]
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
    return ImplicitPrimitiveEquation{NF, VectorType, MatrixType, TensorType}(;
        trunc, nlayers, kwargs...
    )
end

function variables(implicit::ImplicitPrimitiveEquation)
    return (
        GridVariable(:temp_average, VectorDim(implicit.nlayers), desc = "Average vertical temperature profile", units = "K"),
    )
end

# function barrier to unpack the constants struct for primitive eq models
function initialize!(
        I::ImplicitPrimitiveEquation,
        dt::Real,
        vars::Variables,
        model::PrimitiveEquation,
    )
    model.dynamics || return nothing    # escape immediately if no dynamics
    (; geometry, geopotential, atmosphere, adiabatic_conversion) = model
    initialize!(I, dt, vars.grid.temp_average, geometry, geopotential, atmosphere, adiabatic_conversion)
    return nothing
end

"""$(TYPEDSIGNATURES)
Initialize the implicit terms for the PrimitiveEquation models."""
function initialize!(
        implicit::ImplicitPrimitiveEquation,
        dt::Real,                                           # the scaled time step radius*dt
        temp_average::AbstractVector,                       # average vertical temperature profile to construct the operators
        geometry::AbstractGeometry,
        geopotential::AbstractGeopotential,
        atmosphere::AbstractAtmosphere,
        adiabatic_conversion::AbstractAdiabaticConversion,
    )

    NF = eltype(temp_average)

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
    temp_profile, S, S⁻¹, L, R, U, W, L0, L1, L2, L3, L4 =
        on_architecture(
        CPU(), (
            implicit.temp_profile, implicit.S, implicit.S⁻¹, implicit.L,
            implicit.R, implicit.U, implicit.W, implicit.L0, implicit.L1,
            implicit.L2, implicit.L3, implicit.L4,
        )
    )

    # Also transfer geometry and other arrays to CPU
    σ_levels_full, σ_levels_thick, Δp_geopot_half, Δp_geopot_full, σ_lnp_A, σ_lnp_B, temp_average =
        on_architecture(
        CPU(), (
            σ_levels_full, σ_levels_thick, Δp_geopot_half, Δp_geopot_full,
            σ_lnp_A, σ_lnp_B, temp_average,
        )
    )

    # use current vertical temperature profile
    temp_profile .= temp_average

    # return immediately if temp_profile contains NaRs, model blew up in that case
    all(isfinite.(temp_profile)) || return nothing

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

    ξ = α * dt                          # dt = 2Δt for leapfrog, but = Δt, Δ/2 in first_timesteps!
    implicit.ξ[] = ξ                    # also store in Implicit struct

    # DIVERGENCE OPERATORS (called g in Hoskins and Simmons 1975, eq 11 and Appendix 1)
    @inbounds for k in 1:nlayers                # vertical geopotential integration as matrix operator
        R[1:k, k] .= -Δp_geopot_full[k]         # otherwise equivalent to geopotential! with zero orography
        R[1:(k - 1), k] .+= -Δp_geopot_half[k]  # incl the minus but excluding the eigenvalues as with U
    end
    U .= -R_dry * temp_profile                  # the R_d*Tₖ∇² term excl the eigenvalues from ∇² for divergence

    # TEMPERATURE OPERATOR (called τ in Hoskins and Simmons 1975, eq 9 and Appendix 1)
    L0 .= 1 ./ 2σ_levels_thick
    L2 .= κ * temp_profile .* σ_lnp_A    # factor in front of the div_sum_above term
    L4 .= κ * temp_profile .* σ_lnp_B    # factor in front of div term in Dlnps/Dt

    @inbounds for k in 1:nlayers
        Tₖ = temp_profile[k]                    # average temperature at k
        k_above = max(1, k - 1)                 # layer index above
        k_below = min(k + 1, nlayers)           # layer index below
        ΔT_above = Tₖ - temp_profile[k_above]   # temperature difference to layer above
        ΔT_below = temp_profile[k_below] - Tₖ   # and to layer below
        σₖ = σ_levels_full[k]                   # should be Σ_r=1^k Δσᵣ for model top at >0hPa
        σₖ_above = σ_levels_full[k_above]

        for r in 1:nlayers
            L1[k, r] = ΔT_below * σ_levels_thick[r] * σₖ            # vert advection operator below
            L1[k, r] -= k >= r ? σ_levels_thick[r] : zero(NF)

            L1[k, r] += ΔT_above * σ_levels_thick[r] * σₖ_above     # vert advection operator above
            L1[k, r] -= (k - 1) >= r ? σ_levels_thick[r] : zero(NF)
        end

        # _sum_above operator itself
        L3[1:k, k] .= 0                              # fill upper triangle + diagonal with zeros
        L3[(k + 1):end, k] .= σ_levels_thick[k]      # vert integration top to k-1
    end

    L .= Diagonal(L0) * L1 .+ Diagonal(L2) * L3 .+ Diagonal(L4)  # combine all operators into L

    # PRESSURE OPERATOR (called πᵣ in Hoskins and Simmons, 1975 Appendix 1)
    W .= -σ_levels_thick                # the -D̄ term in the log surface pres equation

    # solving the equations above for δD yields
    # δD = SG, with G = G_D + ξRG_T + ξUG_lnps and the operator S
    # S = 1 - ξ²(RL + UW) that has to be inverted to obtain δD from the Gs
    I = LinearAlgebra.I(nlayers)
    @inbounds for l in 1:(trunc + 1)
        eigenvalue = -l * (l - 1)       # 1-based, -l*(l+1) → -l*(l-1)
        S .= I .- ξ^2 * eigenvalue * (R * L .+ U * W')

        # inv(S) but saving memory:
        luS = LinearAlgebra.lu!(S)      # in-place LU decomposition (overwriting S)
        Sinv = L1                       # reuse L1 matrix to store inv(S)
        Sinv .= I                       # use ldiv! so last arg needs to be unity matrix
        LinearAlgebra.ldiv!(luS, Sinv)  # now do S\I = S⁻¹ via LU decomposition
        S⁻¹[l, :, :] .= Sinv            # store in array
    end

    # Transfer computed results back to the original architecture
    # runic: off
    implicit.temp_profile   .= on_architecture(arch, temp_profile)
    implicit.S              .= on_architecture(arch, S)
    implicit.S⁻¹            .= on_architecture(arch, S⁻¹)
    implicit.L              .= on_architecture(arch, L)
    implicit.R              .= on_architecture(arch, R)
    implicit.U              .= on_architecture(arch, U)
    implicit.W              .= on_architecture(arch, W)
    implicit.L0             .= on_architecture(arch, L0)
    implicit.L1             .= on_architecture(arch, L1)
    implicit.L2             .= on_architecture(arch, L2)
    implicit.L3             .= on_architecture(arch, L3)
    implicit.L4             .= on_architecture(arch, L4)
    # runic: on
    return nothing
end

set_initialized!(implicit::ImplicitPrimitiveEquation) = (implicit.initialized = true)

"""$(TYPEDSIGNATURES)
Apply the implicit corrections to dampen gravity waves in the primitive equation models."""
function implicit_correction!(
        vars::Variables,
        implicit::ImplicitPrimitiveEquation,
        model::PrimitiveEquation,
    )

    # escape immediately if explicit
    implicit.α == 0 && return nothing

    (; S⁻¹, R, U, L, W, nlayers) = implicit
    ξ = implicit.ξ[]

    temp_tend = vars.tendencies.temp
    pres_tend = vars.tendencies.pres
    div_tend = vars.tendencies.div
    div_old, div_new = get_steps(vars.prognostic.div)
    G = vars.scratch.a                  # reuse work arrays, used for combined tendency G
    geopotential = vars.scratch.b       # used for geopotential
    geopotential .= 0

    # Get precomputed l_indices from the spectrum
    l_indices = temp_tend.spectrum.l_indices

    arch = architecture(temp_tend)

    # Single kernel: All implicit correction steps for each spectral mode
    launch!(
        arch, LinearWorkOrder, (size(pres_tend, 1),),
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
