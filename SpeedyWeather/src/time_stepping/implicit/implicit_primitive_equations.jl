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
        B,              # Bool
        RefV,           # Base.RefValue{NF}
    } <: AbstractImplicit

    # DIMENSIONS
    "[DERIVED] Spectral resolution"
    trunc::IntType

    "[DERIVED] Number of vertical layers"
    nlayers::IntType

    # PARAMETERS
    "[OPTION] Time-step coefficient: 0.5 = Crank-Nicolson, 1=backward Euler"
    centering::NF = 1.0

    "[DERIVED] Time step [s] used to initialize. Used to check whether time step has changed and reinitialization is needed."
    Δt::RefV = Ref(zero(NF))

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
    S⁻¹::TensorType = zeros(NF, trunc + 2, nlayers, nlayers)
end

"""$(TYPEDSIGNATURES)
Generator using the resolution from SpectralGrid."""
function ImplicitPrimitiveEquation(spectral_grid::SpectralGrid; kwargs...)
    (; NF, VectorType, MatrixType, TensorType, trunc, nlayers) = spectral_grid
    return ImplicitPrimitiveEquation{NF, VectorType, MatrixType, TensorType, typeof(trunc), Bool, Base.RefValue{NF}}(;
        trunc, nlayers, kwargs...
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
    @trace if implicit.Δt[] != Δt                   # if time step has not changed no need to reinitialize
        scale = vars.prognostic.scale[]             # implicit solver needs to be initialized with scaled time step
        Tₖ = vars.dynamics.average_temperature_profile
        initialize!(implicit, Δt / scale, Tₖ, geometry, geopotential, atmosphere, adiabatic_conversion)
        implicit.Δt[] = Δt
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Initialize the implicit terms for the PrimitiveEquation models."""
function initialize!(
        implicit::ImplicitPrimitiveEquation{NF},
        Δt::Real,                                           # the time step [s], scaled
        temp_average::AbstractVector,                       # average vertical temperature profile to construct the operators
        geometry::AbstractGeometry,
        geopotential::AbstractGeopotential,
        atmosphere::AbstractAtmosphere,
        adiabatic_conversion::AbstractAdiabaticConversion,
    ) where {NF}
    (; trunc, nlayers) = implicit
    (; σ_levels_full, σ_levels_thick) = geometry
    (; R_dry, κ) = atmosphere
    (; Δp_geopot_half, Δp_geopot_full) = geopotential
    (; σ_lnp_A, σ_lnp_B) = adiabatic_conversion

    arch = architecture(implicit.S)

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
    temp_profile .= temp_average

    # return immediately if temp_profile contains NaRs, model blew up in that case
    # TODO: reactive when issues with Reactant resolved
    # all(isfinite.(temp_profile)) || return nothing

    @assert 0.5 <= implicit.centering <= 1 "Centering coefficient must be between 0.5 (centred implicit) and 1 (backward implicit)"
    ξ = implicit.centering * Δt                 # 2Δt for leapfrog, but = Δt, Δ/2 in first_timesteps!

    # index vectors for broadcasting: rows = 1:nlayers (column), cols = (1:nlayers)' (row)
    rows = (1:nlayers)
    cols = (1:nlayers)'

    # DIVERGENCE OPERATORS (called g in Hoskins and Simmons 1975, eq 11 and Appendix 1)
    # R[row, k] = -Δp_geopot_full[k] for row <= k, additionally -Δp_geopot_half[k] for row < k
    R .= .-Δp_geopot_full[cols] .* (rows .<= cols) .- Δp_geopot_half[cols] .* (rows .< cols)

    # U = -R_dry * temp_profile (the R_d*Tₖ∇² term excl eigenvalues from ∇² for divergence)
    U .= .-R_dry .* temp_profile

    # TEMPERATURE OPERATOR (called τ in Hoskins and Simmons 1975, eq 9 and Appendix 1)
    L0 .= 1 ./ (2 .* σ_levels_thick)
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
        inv_pivot = 1/S_scratch[l, pivot, pivot]   #TODO: `inv` isn't compatible with Reactant yet, add it back once that's done

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

    (; S⁻¹, R, U, L, W, nlayers) = implicit
    
    # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog)
    # dynamical core uses scaled time step, scale on the fly
    Δt = time_step(time_stepping, vars.prognostic.clock)       
    ξ = implicit.centering * Δt / vars.prognostic.scale[]

    temp_tend = get_tendency_step(vars.tendencies.temperature, time_stepping, implicit)
    pres_tend = get_tendency_step(vars.tendencies.pressure, time_stepping, implicit)
    div_tend = get_tendency_step(vars.tendencies.divergence, time_stepping, implicit)
    div_old, div_new = get_steps(vars.prognostic.divergence)
    G = vars.scratch.a                  # reuse work arrays, used for combined tendency G
    geopotential = vars.scratch.b       # used for geopotential
    geopotential .= 0

    # Get precomputed l_indices from the spectrum
    l_indices = temp_tend.spectrum.l_indices

    arch = architecture(temp_tend)

    # Single kernel: All implicit correction steps for each spectral mode
    launch!(
        arch, LinearWorkOrder, (size(pres_tend, 1),),
        implicit_primitive_leapfrog_kernel!,
        temp_tend, pres_tend, div_tend, G, geopotential,
        div_old, div_new, S⁻¹, R, U, L, W, l_indices, ξ, nlayers
    )

    return nothing
end

# Single kernel that does all steps for one spectral mode
@kernel inbounds = true function implicit_primitive_leapfrog_kernel!(
        temp_tend, pres_tend, div_tend, G, geopotential,
        div_old, div_new, S⁻¹, R, U, L, W, l_indices,
        ξ, nlayers
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