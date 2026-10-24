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
        ScheduleType,
        RefNF,          # <: Base.RefValue{NF} or Base.RefValue{<:Number} for Reactant, see `scalar_ref`
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

    "[OPTION] Schedule to recompute the reference temperature profile (and the implicit operators) from the current global-mean temperature"
    schedule::ScheduleType = PeriodicSchedule(every = Day(1))

    "[DERIVED] Time step [s] used to initialize. Used to check whether time step has changed and reinitialization is needed."
    Δt::RefNF = Ref(zero(NF))

    # PRECOMPUTED ARRAYS, to be initialized with initialize!
    "[DERIVED] Reference temperature profile, the global-mean temperature when last (re)initialized"
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
function ImplicitPrimitiveEquation(
        spectral_grid::SpectralGrid;
        schedule::AbstractSchedule = PeriodicSchedule(every = Day(1)),
        kwargs...
    )
    (; NF, VectorType, MatrixType, TensorType, truncation, nlayers, architecture) = spectral_grid
    IntType = typeof(truncation)
    Δt = scalar_ref(architecture, zero(NF))
    return ImplicitPrimitiveEquation{NF, VectorType, MatrixType, TensorType, IntType, nlayers, typeof(schedule), typeof(Δt)}(;
        truncation, nlayers, schedule, Δt, kwargs...
    )
end

function variables(implicit::ImplicitPrimitiveEquation)
    return (
        DynamicsVariable(:average_temperature_profile, VectorDim(implicit.nlayers), desc = "Average vertical temperature profile", units = "K"),
    )
end

"""$(TYPEDSIGNATURES)
Initialize the `schedule` of `implicit` with the time step of `model`.
The implicit operators themselves are only built in the first time steps by `reinitialize!`,
as they need the temperature profile of the initial conditions and the time step of the
leapfrog start."""
function initialize!(implicit::ImplicitPrimitiveEquation, model::PrimitiveEquation)
    initialize!(implicit.schedule, model.time_stepping)
    return nothing
end

# function barrier to decide whether to (re)initialize: when the time step changes (the leapfrog
# start with Δt/2, Δt, then 2Δt, or a new time step) or when the reference temperature profile
# is scheduled to be recomputed from the current global-mean temperature
function reinitialize!(
        implicit::ImplicitPrimitiveEquation,
        model::PrimitiveEquation,
        vars::Variables,
    )
    (; time_stepping, geometry, geopotential, atmosphere, adiabatic_conversion) = model
    (; clock) = vars.prognostic
    Δt = time_step(time_stepping, clock)
    Δt_changed = implicit.Δt[] != Δt
    Δt_changed || isscheduled(implicit.schedule, clock) || return nothing

    scale = vars.prognostic.scale[]                         # implicit solver needs to be initialized with scaled time step
    Tₖ = vars.dynamics.average_temperature_profile          # current global-mean temperature profile
    initialize!(implicit, Δt / scale, Tₖ, geometry, geopotential, atmosphere, adiabatic_conversion)
    implicit.Δt[] = Δt
    return nothing
end

"""$(TYPEDSIGNATURES)
Initialize the implicit terms for the PrimitiveEquation models. The operators are constructed
via broadcasting on the device, without transfers to the CPU, scalar indexing, or data-dependent
control flow so that Reactant can trace this inside the time step. S⁻¹ is then computed per
degree l with `invert_S!`."""
function initialize!(
        implicit::ImplicitPrimitiveEquation,
        Δt::Number,                                         # the time step [s], scaled
        temp_average::AbstractVector,                       # average vertical temperature profile to construct the operators
        geometry::AbstractGeometry,
        geopotential::AbstractGeopotential,
        atmosphere::AbstractAtmosphere,
        adiabatic_conversion::AbstractAdiabaticConversion,
    )
    (; nlayers) = implicit
    (; σ_levels_full, σ_levels_thick) = geometry
    (; R_dry, κ) = atmosphere
    (; Δp_geopot_half, Δp_geopot_full) = geopotential
    (; σ_lnp_A, σ_lnp_B) = adiabatic_conversion
    (; temp_profile, L, R, U, W, L0, L1, L2, L3, L4) = implicit

    # use current vertical temperature profile
    temp_profile .= temp_average

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

    # index vectors for broadcasting: rows = 1:nlayers (column), cols = (1:nlayers)' (row)
    rows = (1:nlayers)
    cols = (1:nlayers)'

    # DIVERGENCE OPERATORS (called g in Hoskins and Simmons 1975, eq 11 and Appendix 1)
    # vertical geopotential integration as matrix operator, otherwise equivalent to geopotential! with zero orography
    # R[row, k] = -Δp_geopot_full[k] for row <= k, additionally -Δp_geopot_half[k] for row < k
    R .= .-Δp_geopot_full[cols] .* (rows .<= cols) .- Δp_geopot_half[cols] .* (rows .< cols)
    U .= .-R_dry .* temp_profile                # the R_d*Tₖ∇² term excl the eigenvalues from ∇² for divergence

    # TEMPERATURE OPERATOR (called τ in Hoskins and Simmons 1975, eq 9 and Appendix 1)
    L0 .= 1 ./ (2 .* σ_levels_thick)
    L2 .= κ .* temp_profile .* σ_lnp_A          # factor in front of the div_sum_above term
    L4 .= κ .* temp_profile .* σ_lnp_B          # factor in front of div term in Dlnps/Dt

    # L1[k, r]: vertical advection operator, layer indices above/below clamped to 1:nlayers
    k_above = max.(1, rows .- 1)
    k_below = min.(rows .+ 1, nlayers)
    ΔT_below = temp_profile[k_below] .- temp_profile[rows]
    ΔT_above = temp_profile[rows] .- temp_profile[k_above]
    σₖ = σ_levels_full[rows]                    # should be Σ_r=1^k Δσᵣ for model top at >0hPa
    σₖ_above = σ_levels_full[k_above]
    Δσᵣ = σ_levels_thick[cols]

    L1 .= ΔT_below .* Δσᵣ .* σₖ .- Δσᵣ .* (rows .>= cols) .+
        ΔT_above .* Δσᵣ .* σₖ_above .- Δσᵣ .* ((rows .- 1) .>= cols)

    # L3[r, k]: _sum_above operator itself, lower triangle with σ_levels_thick[k]
    L3 .= σ_levels_thick[cols] .* (rows .> cols)

    # combine all operators into L = Diagonal(L0)*L1 + Diagonal(L2)*L3 + Diagonal(L4)
    L .= L0 .* L1 .+ L2 .* L3 .+ L4 .* (rows .== cols)

    # PRESSURE OPERATOR (called πᵣ in Hoskins and Simmons, 1975 Appendix 1)
    W .= .-σ_levels_thick                       # the -D̄ term in the log surface pres equation

    # solving the equations above for δD yields
    # δD = SG, with G = G_D + ξRG_T + ξUG_lnps and the operator S
    # S = 1 - ξ²(RL + UW) that has to be inverted to obtain δD from the Gs
    invert_S!(implicit, ξ, architecture(implicit.S⁻¹))
    return nothing
end

# S⁻¹ via LU decomposition on the CPU, but Reactant traces `initialize!` inside the time step
# so no transfer to the CPU is possible, use Gauss-Jordan elimination in a kernel instead
invert_S!(implicit::ImplicitPrimitiveEquation, ξ, ::AbstractArchitecture) = invert_S_lu!(implicit, ξ)
invert_S!(implicit::ImplicitPrimitiveEquation, ξ, ::ReactantDevice) = invert_S_gauss_jordan!(implicit, ξ)

"""$(TYPEDSIGNATURES)
Compute S⁻¹ = (1 - ξ²λₗ(RL + UWᵀ))⁻¹ for every degree l, with λₗ the eigenvalues of ∇², via
LU decomposition on the CPU, then transfer S⁻¹ back to the device. These are small
(nlayers × nlayers) matrices, so CPU computation is more efficient."""
function invert_S_lu!(implicit::ImplicitPrimitiveEquation, ξ)
    (; truncation, nlayers) = implicit
    arch = architecture(implicit.S⁻¹)
    temp_profile, S, S⁻¹, R, L, U, W = on_architecture(
        CPU(), (implicit.temp_profile, implicit.S, implicit.S⁻¹, implicit.R, implicit.L, implicit.U, implicit.W)
    )

    # return immediately if temp_profile contains NaRs, model blew up in that case
    all(isfinite, temp_profile) || return nothing

    I = LinearAlgebra.I(nlayers)
    RL_UW = R * L .+ U * W'             # same for every degree l
    S⁻¹_l = similar(S)                  # inverse for one degree l
    @inbounds for l in 1:truncation     # 1-based degree
        eigenvalue = -l * (l - 1)       # 1-based, -l*(l+1) → -l*(l-1)
        S .= I .- ξ^2 * eigenvalue * RL_UW

        # inv(S) but saving memory:
        luS = LinearAlgebra.lu!(S)      # in-place LU decomposition (overwriting S)
        S⁻¹_l .= I                      # use ldiv! so last arg needs to be unity matrix
        LinearAlgebra.ldiv!(luS, S⁻¹_l) # now do S\I = S⁻¹ via LU decomposition
        S⁻¹[l, :, :] .= S⁻¹_l           # store in array, S⁻¹[truncation + 1, :, :] remains unused
    end

    implicit.S⁻¹ .= on_architecture(arch, S⁻¹)
    return nothing
end

"""$(TYPEDSIGNATURES)
Compute S⁻¹ as in `invert_S_lu!` but on the device via Gauss-Jordan elimination (without pivoting)
in a kernel, one thread per degree l, so that Reactant can trace it."""
function invert_S_gauss_jordan!(implicit::ImplicitPrimitiveEquation, ξ)
    (; truncation, nlayers, S⁻¹, R, L, U, W) = implicit

    # S⁻¹[truncation + 1, :, :] remains unused as in `invert_S_lu!`
    S_scratch = similar(S⁻¹)
    launch!(
        architecture(S⁻¹), LinearWorkOrder, (truncation,), implicit_invert_S_kernel!,
        S⁻¹, S_scratch, R, L, U, W, ξ, nlayers
    )
    return nothing
end

# Compute S⁻¹ for each spectral degree l via Gauss-Jordan elimination, S = I - ξ²*eigenvalue*(R*L + U*Wᵀ)
# S_scratch[l, :, :] stores S during elimination, S⁻¹[l, :, :] tracks the identity → inverse transformation
@kernel inbounds = true function implicit_invert_S_kernel!(
        S⁻¹,                            # Output: (truncation+1) × nlayers × nlayers tensor
        S_scratch,                      # Scratch: same shape as S⁻¹, stores S per l
        @Const(R),                      # Input: divergence operator matrix
        @Const(L),                      # Input: temperature operator matrix
        @Const(U),                      # Input: divergence vector
        @Const(W),                      # Input: pressure vector
        @Const(ξ),                      # Input: semi-implicit time step coefficient
        @Const(nlayers),                # Input: number of layers
    )
    l = @index(Global, Linear)          # spectral degree (1-based)

    NF = eltype(S⁻¹)
    eigenvalue = -l * (l - 1)           # 1-based, -l*(l+1) → -l*(l-1)
    ξ²λ = ξ^2 * eigenvalue

    # S = I - ξ²λ*(R*L + U*Wᵀ) into S_scratch[l, :, :] and initialize S⁻¹[l, :, :] as identity
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

    # Gauss-Jordan elimination: reduce S_scratch[l, :, :] to I, applying the same row operations
    # to S⁻¹[l, :, :] (I → S⁻¹)
    for pivot in 1:nlayers
        pivot⁻¹ = 1 / S_scratch[l, pivot, pivot]    # TODO: `inv` fails in kernels on Reactant+CUDA https://github.com/EnzymeAD/Reactant.jl/issues/2560

        # scale pivot row
        for r in 1:nlayers
            S_scratch[l, pivot, r] *= pivot⁻¹
            S⁻¹[l, pivot, r] *= pivot⁻¹
        end

        # eliminate all other rows
        for k in 1:nlayers
            factor = ifelse(k == pivot, zero(NF), S_scratch[l, k, pivot])
            for r in 1:nlayers
                S_scratch[l, k, r] -= factor * S_scratch[l, pivot, r]
                S⁻¹[l, k, r] -= factor * S⁻¹[l, pivot, r]
            end
        end
    end
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
        for r in k:nlayers
            geopotential_val += @inbounds R[k, r] * temp_tend[lm, r, temp_step]
        end
        geopotential_val
    end

    # Step 3: Calculate the G = G_D + ξRG_T + ξUG_lnps terms
    # ∇² not part of U so *eigenvalues here
    G = ntuple(Val(nlayers)) do k
        @inbounds div_tend[lm, k, div_step] + ξ * eigenvalue * (U[k] * pres_tend[lm, pres_step] + geopotential[k])
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
