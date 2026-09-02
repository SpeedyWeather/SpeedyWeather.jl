export ImplicitPrimitiveEquation

"""
Struct that holds various precomputed arrays for the semi-implicit correction to
prevent gravity waves from amplifying in the primitive equation model.
$(TYPEDFIELDS)"""
@kwdef struct ImplicitPrimitiveEquation{
        NF,             # number format
        VectorType,
        MatrixType,
        TensorType,
        IntType,
        NLAYERS,        # nlayers as a type parameter, see `nlayers_val` below
    } <: AbstractImplicit

    # DIMENSIONS
    "[DERIVED] Spectral resolution"
    truncation::IntType

    "[DERIVED] Number of vertical layers"
    nlayers::IntType

    "[DERIVED] `nlayers` as a `Val` so that the vertical loops are more type stable"
    nlayers_val::Val{NLAYERS} = Val(nlayers)

    # PARAMETERS
    "[OPTION] Time-step coefficient: 0.5 = Crank-Nicolson, 1=backward Euler"
    centering::NF = 1.0

    "[DERIVED] Time step [s] used to initialize. Used to check whether time step has changed and reinitialization is needed."
    Δt::Base.RefValue{NF} = Ref(zero(NF))

    # PRECOMPUTED ARRAYS, to be initialized with initialize!
    "[DERIVED] vertical temperature profile, obtained from diagn on first time step"
    temp_profile::VectorType = zeros(NF, nlayers)

    "[DERIVED] divergence: operator for the geopotential calculation"
    R::MatrixType = zeros(NF, nlayers, nlayers)

    "[DERIVED] divergence: the -RdTₖ∇² term excl the eigenvalues from ∇² for divergence"
    U::VectorType = zeros(NF, nlayers)

    "[DERIVED] temperature: operator for the TₖD + κTₖD(ln pₛ)/Dt term"
    L::MatrixType = zeros(NF, nlayers, nlayers)

    "[DERIVED] pressure: vertical averaging of the -D̄ term in the log surface pressure equation"
    W::VectorType = zeros(NF, nlayers)

    "[DERIVED] components to construct L, 1/ 2Δσ"
    L0::VectorType = zeros(NF, nlayers)

    "[DERIVED] vert advection term in the temperature equation (below+above)"
    L1::MatrixType = zeros(NF, nlayers, nlayers)

    "[DERIVED] factor in front of the `div_sum_above` term"
    L2::VectorType = zeros(NF, nlayers)

    "[DERIVED] `_sum_above` operator itself"
    L3::MatrixType = zeros(NF, nlayers, nlayers)

    "[DERIVED] factor in front of div term in Dlnpₛ/Dt"
    L4::VectorType = zeros(NF, nlayers)

    "[DERIVED] for every l the matrix to be inverted"
    S::MatrixType = zeros(NF, nlayers, nlayers)

    "[DERIVED] combined inverted operator: S = 1 - ξ²(RL + UW)"
    S⁻¹::TensorType = zeros(NF, truncation + 1, nlayers, nlayers)
end

"""$(TYPEDSIGNATURES)
Generator using the resolution from SpectralGrid."""
function ImplicitPrimitiveEquation(spectral_grid::SpectralGrid; kwargs...)
    (; NF, VectorType, MatrixType, TensorType, truncation, nlayers) = spectral_grid
    return ImplicitPrimitiveEquation{NF, VectorType, MatrixType, TensorType, typeof(truncation), nlayers}(;
        truncation, nlayers, kwargs...
    )
end

function variables(implicit::ImplicitPrimitiveEquation)
    return (
        DynamicsVariable(:average_temperature_profile, VectorDim(implicit.nlayers), desc = "Average vertical temperature profile", units = "K"),
    )
end

# function barrier to decide whether to initialize or not based on time step
function reinitialize!(
        implicit::ImplicitPrimitiveEquation,
        model::PrimitiveEquation,
        vars::Variables,
    )
    (; time_stepping, geometry, geopotential, atmosphere, adiabatic_conversion) = model
    Δt = time_step(time_stepping, vars.prognostic.clock)
    implicit.Δt[] == Δt && return nothing                   # if time step has not changed no need to reinitialize
    scale = vars.prognostic.scale[]                         # implicit solver needs to be initialized with scaled time step
    Tₖ = vars.dynamics.average_temperature_profile
    initialize!(implicit, Δt / scale, Tₖ, geometry, geopotential, atmosphere, adiabatic_conversion)
    implicit.Δt[] = Δt
    return nothing
end

"""$(TYPEDSIGNATURES)
Initialize the implicit terms for the PrimitiveEquation models."""
function initialize!(
        implicit::ImplicitPrimitiveEquation,
        Δt::Real,                                           # the time step [s], scaled
        temp_average::AbstractVector,                       # average vertical temperature profile to construct the operators
        geometry::AbstractGeometry,
        geopotential::AbstractGeopotential,
        atmosphere::AbstractAtmosphere,
        adiabatic_conversion::AbstractAdiabaticConversion,
    )

    NF = eltype(temp_average)

    (; truncation, nlayers) = implicit
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

    @assert 0.5 <= implicit.centering <= 1 "Centering coefficient must be between 0.5 (centred implicit) and 1 (backward implicit)"
    ξ = implicit.centering * Δt                 # 2Δt for leapfrog, but = Δt, Δ/2 in first_timesteps!

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
    @inbounds for l in 1:truncation     # 1-based degree
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

"""$(TYPEDSIGNATURES)
Apply the implicit corrections to dampen gravity waves in the primitive equation models."""
function implicit_correction!(
        vars::Variables,
        implicit::ImplicitPrimitiveEquation,
        time_stepping::AbstractLeapfrog,
        model::PrimitiveEquation,
    )

    # escape immediately if explicit
    implicit.centering == 0 && return nothing

    (; S⁻¹, R, U, L, W, nlayers_val) = implicit

    # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog)
    # dynamical core uses scaled time step, scale on the fly
    Δt = time_step(time_stepping, vars.prognostic.clock)
    ξ = implicit.centering * Δt / vars.prognostic.scale[]

    # NOTE: the FULL arrays are handed to the kernel together with their step
    # indices, and the step dimension is indexed inside the kernel. Taking
    # `get_step`/`get_tendency_step` views here instead makes Enzyme's type
    # analysis give up on Julia 1.12 (`EnzymeNoTypeError` in `lta_view`), because
    # these variables are themselves already views into a fused parent, so the
    # step accessor would build a view of a view.
    temperature_tendency = vars.tendencies.temperature
    pressure_tendency = vars.tendencies.pressure
    divergence_tendency = vars.tendencies.divergence
    divergence = vars.prognostic.divergence

    temp_step = which_tendency_step(temperature_tendency, time_stepping, implicit)
    pres_step = which_tendency_step(pressure_tendency, time_stepping, implicit)
    div_step = which_tendency_step(divergence_tendency, time_stepping, implicit)

    # Get precomputed l_indices from the spectrum
    l_indices = temperature_tendency.spectrum.l_indices

    arch = architecture(temperature_tendency)

    launch!(
        arch, LinearWorkOrder, (size(pressure_tendency, 1),),
        implicit_primitive_leapfrog_kernel!,
        temperature_tendency, pressure_tendency, divergence_tendency,
        divergence, temp_step, pres_step, div_step,
        S⁻¹, R, U, L, W, l_indices, ξ, nlayers_val
    )

    return nothing
end

# Single kernel that does all steps for one spectral mode
@kernel inbounds = true function implicit_primitive_leapfrog_kernel!(
        temp_tend, pres_tend, div_tend,
        divergence, temp_step, pres_step, div_step,
        S⁻¹, R, U, L, W, l_indices,
        ξ, ::Val{nlayers}
    ) where {nlayers}
    lm = @index(Global, Linear)

    # Get degree l for this spectral mode
    l = l_indices[lm]

    # Step 1: Move implicit terms of temperature equation from time step i to i-1
    # RHS_expl(Vⁱ) + RHS_impl(Vⁱ⁻¹) = RHS(Vⁱ) + RHS_impl(Vⁱ⁻¹ - Vⁱ)
    for k in 1:nlayers
        temp_tend_val = zero(eltype(temp_tend))
        for r in 1:nlayers
            temp_tend_val += L[k, r] * (divergence[lm, r, 1] - divergence[lm, r, 2])
        end
        temp_tend[lm, k, temp_step] += temp_tend_val
    end

    eigenvalue = -l * (l - 1)  # 1-based, -l*(l+1) → -l*(l-1)

    # Step 2: Calculate the ξ*R*G_T term, vertical integration of geopotential
    # (excl ξ, this is done in step 3). Held in an `ntuple` (registers) rather than a
    # scratch array, see the note at the launch site.
    geopotential = ntuple(Val(nlayers)) do k
        # skip 1:k-1 as integration is surface to k
        geopotential_val = zero(eltype(temp_tend))
        # while loop instead of `for r in k:nlayers`: the triangular range with both
        # endpoints dynamic miscompiles on AMDGPU/CDNA (gfx90a/gfx942), see
        # https://github.com/JuliaGPU/AMDGPU.jl/issues/1015
        r = k
        while r <= nlayers
            geopotential_val += R[k, r] * temp_tend[lm, r, temp_step]
            r += 1
        end
        geopotential_val
    end

    # Step 3: Calculate the G = G_D + ξRG_T + ξUG_lnps terms
    # ∇² not part of U so *eigenvalues here
    G = ntuple(Val(nlayers)) do k
        div_tend[lm, k, div_step] + ξ * eigenvalue * (U[k] * pres_tend[lm, pres_step] + geopotential[k])
    end

    # Step 4: Now solve δD = S⁻¹G to correct divergence tendency
    for k in 1:nlayers
        div_val = zero(eltype(div_tend))
        for r in 1:nlayers
            div_val += S⁻¹[l, k, r] * G[r]
        end
        div_tend[lm, k, div_step] = div_val
    end

    # Step 5: Semi implicit corrections for temperature and pressure

    # Step 5a: Temperature correction δT = G_T + ξLδD
    for k in 1:nlayers
        temp_correction = zero(eltype(temp_tend))
        for r in 1:nlayers
            temp_correction += ξ * L[k, r] * div_tend[lm, r, div_step]
        end
        temp_tend[lm, k, temp_step] += temp_correction
    end

    # Step 5b: Pressure correction δlnpₛ = G_lnpₛ + ξWδD
    pres_correction = zero(eltype(pres_tend))
    for k in 1:nlayers
        pres_correction += ξ * W[k] * div_tend[lm, k, div_step]
    end
    pres_tend[lm, pres_step] += pres_correction
end
