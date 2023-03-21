# BAROTROPIC MODEL
"""
    nothing = initialize_implicit!(::Real,::Barotropic)

Just passes, as implicit terms are not used in the barotropic model."""
initialize_implicit!(::Real,::Barotropic) = nothing

# !!! structs are defined in define_implicit.jl !!!

# SHALLOW WATER MODEL
"""
    I = Implicit(P::Parameters{<:ShallowWater})

Zero generator function for an `ImplicitShallowWater` struct, which holds precomputed arrays
for the implicit correction in the shallow water model. Actual precomputation happens in
initialize_implicit!."""
function Implicit(P::Parameters{<:ShallowWater})    # shallow water model only
    
    @unpack NF,trunc = P

    # 0-initialize, actual initialization depends on time step, done in initialize_implicit!
    ξH₀ = [zero(NF)]            # time step ξ, layer thickness at rest H₀ (in vec for mutability)
    g∇² = zeros(NF,trunc+2)     # gravity times Laplace operator
    ξg∇² = zeros(NF,trunc+2)    # time step ξ times gravity times Laplace operator 
    S⁻¹ = zeros(NF,trunc+2)     # combined operator to be inverted

    return ImplicitShallowWater(ξH₀,g∇²,ξg∇²,S⁻¹)
end

"""
    initialize_implicit!(dt::Real,M::BarotropicModel)

Update the implicit terms in `M` for the shallow water model as they depend on the time step `dt`."""
function initialize_implicit!(  dt::Real,               # time step
                                model::ShallowWater)    # update Implicit struct in M

    @unpack implicit_α = model.parameters               # = [0,0.5,1], time step fraction for implicit
    @unpack eigenvalues = model.spectral_transform      # = -l*(l+1), degree l of harmonics
    @unpack ξH₀,g∇²,ξg∇²,S⁻¹ = model.implicit           # pull precomputed arrays to be updated
    @unpack layer_thickness = model.constants           # shallow water layer thickness [m]
    @unpack gravity = model.constants                   # gravitational acceleration [m/s²]                  

    # implicit time step between i-1 and i+1
    # α = 0   means the gravity wave terms are evaluated at i-1 (forward)
    # α = 0.5 evaluates at i+1 and i-1 (centered implicit)
    # α = 1   evaluates at i+1 (backward implicit)
    # α ∈ [0.5,1] are also possible which controls the strength of the gravity wave dampening.
    # α = 0.5 slows gravity waves and prevents them from amplifying
    # α > 0.5 will dampen the gravity waves within days to a few timesteps (α=1)

    ξ = implicit_α*dt               # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog) from input dt
    ξH₀[1] = ξ*layer_thickness      # update ξ*H₀ with new ξ, in vec for mutability

    # loop over degree l of the harmonics (implicit terms are independent of order m)
    @inbounds for l in eachindex(g∇²,ξg∇²,S⁻¹,eigenvalues)
        g∇²[l] = gravity*eigenvalues[l]     # doesn't actually change with dt
        ξg∇²[l] = ξ*g∇²[l]                  # update ξg∇² with new ξ
        S⁻¹[l] = inv(1 - ξH₀[1]*ξg∇²[l])   # update 1/(1-ξ²gH₀∇²) with new ξ
    end
end

"""
    implicit_correction!(   diagn::DiagnosticVariablesLayer,
                            progn::PrognosticVariablesLeapfrog,
                            surface::SurfaceVariables,
                            pres::SurfaceLeapfrog,
                            M::ShallowWaterModel)

Apply correction to the tendencies in `diag` to prevent the gravity waves from amplifying.
The correction is implicitly evaluated using the parameter `implicit_α` to switch between
forward, centered implicit or backward evaluation of the gravity wave terms."""
function implicit_correction!(  diagn::DiagnosticVariablesLayer{NF},
                                progn::PrognosticVariablesLeapfrog{NF},
                                surface::SurfaceVariables{NF},
                                pres::SurfaceLeapfrog{NF},
                                model::ShallowWater,
                                ) where NF

    @unpack div_tend = diagn.tendencies     # divergence tendency
    div_old = progn.leapfrog[1].div         # divergence at t
    div_new = progn.leapfrog[2].div         # divergence at t+dt
    pres_old, pres_new = pres.leapfrog      # pressure/η at t,t+dt
    @unpack pres_tend = surface             # tendency of pressure/η

    @unpack g∇²,ξg∇²,S⁻¹ = model.implicit
    ξH₀ = model.implicit.ξH₀[1]                 # unpack as it's stored in a vec for mutation
    H₀ = model.constants.layer_thickness

    @boundscheck size(div_old) == size(div_new) || throw(BoundsError)
    @boundscheck size(pres_old) == size(pres_new) || throw(BoundsError)
    @boundscheck size(div_tend) == size(pres_new) || throw(BoundsError)
    @boundscheck size(pres_tend) == size(div_new) || throw(BoundsError)
    lmax,mmax = size(div_tend) .- (2,1)

    @boundscheck length(S⁻¹) == lmax+2 || throw(BoundsError)
    @boundscheck length(ξg∇²) == lmax+2 || throw(BoundsError)
    @boundscheck length(g∇²) == lmax+2 || throw(BoundsError)

    lm = 0
    @inbounds for m in 1:mmax+1
        for l in m:lmax+1
            lm += 1     # single index lm corresponding to harmonic l,m with a LowerTriangularMatrix
            
            # calculate the G = N(Vⁱ) + NI(Vⁱ⁻¹ - Vⁱ) term.
            # Vⁱ is a prognostic variable at time step i
            # N is the right hand side of ∂V\∂t = N(V)
            # NI is the part of N that's calculated semi-implicitily: N = NE + NI
            G_div = div_tend[lm] - g∇²[l]*(pres_old[lm] - pres_new[lm])
            G_η   = pres_tend[lm] - H₀*(div_old[lm] - div_new[lm])

            # using the Gs correct the tendencies for semi-implicit time stepping
            div_tend[lm] = S⁻¹[l]*(G_div - ξg∇²[l]*G_η)
            pres_tend[lm] = G_η - ξH₀*div_tend[lm]
        end
        lm += 1     # loop skips last row
    end
end

# PRIMITIVE EQUATION MODEL
"""
    I = Implicit(P::Parameters{<:PrimitiveEquation})

Zero generator function for an `ImplicitPrimitiveEq` struct, which holds precomputed arrays for
the implicit correction in the primitive equation model. Actual precomputation happens in initialize_implicit!."""
function Implicit(P::Parameters{<:PrimitiveEquation})    # primitive equation only
    
    @unpack NF,trunc,nlev = P

    # initialize with zeros only, actual initialization depends on time step, done in initialize_implicit!
    ξ = zeros(NF,1)             # time step 2α*Δt packed in a vector for mutability
    R = zeros(NF,nlev,nlev)     # divergence: operator for the geopotential calculation
    U = zeros(NF,nlev)          # divergence: the -RdTₖ∇² term excl the eigenvalues from ∇² for divergence
    L = zeros(NF,nlev,nlev)     # temperature: operator for the TₖD + κTₖDlnps/Dt term
    W = zeros(NF,nlev)          # pressure: vertical averaging of the -D̄ term in the log surface pres equation
    S⁻¹ = zeros(NF,trunc+1,nlev,nlev)   # combined inverted operator: S = 1 - ξ²(RL + UW)
    return ImplicitPrimitiveEq(ξ,R,U,L,W,S⁻¹)
end

function initialize_implicit!(  dt::Real,                   # the scaled time step radius*dt
                                model::PrimitiveEquation)

    @unpack S⁻¹,L,R,U,W = model.implicit
    @unpack nlev, σ_levels_full, σ_levels_thick, temp_ref_profile = model.geometry
    @unpack Δp_geopot_half, Δp_geopot_full, σ_lnp_A, σ_lnp_B = model.geometry
    @unpack R_dry, κ = model.constants
    @unpack eigenvalues, lmax = model.spectral_transform
    α = model.parameters.implicit_α

    # set up R, U, L, W operators from
    # δD = G_D + ξ(RδT + Uδlnps)        divergence D correction
    # δT = G_T + ξLδD                   temperature T correction
    # δlnps = G_lnps + ξWδD             log surface pressure lnps correction
    # 
    # G_X is the uncorrected explicit tendency calculated as RHS_expl(Xⁱ) + RHS_impl(Xⁱ⁻¹)
    # with RHS_expl being the nonlinear terms calculated from the centered time step i
    # and RHS_impl are the linear terms that are supposed to be calcualted semi-implicitly
    # however, they have sofar only been evaluated explicitly at time step i-1
    # and are subject to be corrected to δX following the equations above
    # R, U, L, W are linear operators that are therefore defined here and inverted
    # to obtain δD first, and then δT and δlnps through substitution

    ξ = α*dt                            # dt = 2Δt for leapfrog, but = Δt, Δ/2 in first_timesteps!
    model.implicit.ξ[1] = ξ             # also store in Implicit struct
    
    # DIVERGENCE OPERATORS (called g in Hoskins and Simmons 1975, eq 11 and Appendix 1)
    for k in 1:nlev                         # vertical geopotential integration as matrix operator
        R[1:k,k] .= -Δp_geopot_full[k]      # otherwise equivalent to geopotential! with zero orography
        R[1:k-1,k] .+= -Δp_geopot_half[k]   # incl the minus but excluding the eigenvalues as with U 
    end
    U .= -R_dry*temp_ref_profile        # the R_d*Tₖ∇² term excl the eigenvalues from ∇² for divergence
    
    # TEMPERATURE OPERATOR (called τ in Hoskins and Simmons 1975, eq 9 and Appendix 1)
    L0 = 1 ./ 2σ_levels_thick
    L1 = zero(L)                        # vert advection term in the temperature equation (below)
    L2 = zero(L)                        # vert advection term in the temperature equation (above)
    L3 = κ*temp_ref_profile.*σ_lnp_A    # factor in front of the div_sum_above term
    L4 = zero(L)                        # _sum_above operator itself
    L5 = κ*temp_ref_profile.*σ_lnp_B    # factor in front of div term in Dlnps/Dt

    for k in 1:nlev
        Tₖ = temp_ref_profile[k]                    # reference temperature at k
        k_above = max(1,k-1)                        # layer index above
        k_below = min(k+1,nlev)                     # layer index below
        ΔT_above = Tₖ - temp_ref_profile[k_above]   # temperature difference to layer above
        ΔT_below = temp_ref_profile[k_below] - Tₖ   # and to layer below
        σₖ = σ_levels_full[k]                       # should be Σ_r=1^k Δσᵣ for model top at >0hPa
        σₖ_above = σ_levels_full[k_above]

        for r in 1:nlev
            L1[k,r] = ΔT_below*σ_levels_thick[r]*σₖ         # vert advection operator below
            L1[k,r] -= k>=r ? σ_levels_thick[r] : 0

            L2[k,r] = ΔT_above*σ_levels_thick[r]*σₖ_above   # vert advection operator above
            L2[k,r] -= (k-1)>=r ? σ_levels_thick[r] : 0
        end

        L4[1:k,k] .= 0                  # fill upper triangle + diagonal with zeros
        L4[k+1:end,k] .= σ_levels_thick[k]                  # vert integration top to k-1
    end

    L .= Diagonal(L0)*(L1+L2) + Diagonal(L3)*L4 + Diagonal(L5)  # combine all operators into L

    # PRESSURE OPERATOR (called πᵣ in Hoskins and Simmons, 1975 Appendix 1)
    W .= -σ_levels_thick                # the -D̄ term in the log surface pres equation

    # solving the equations above for δD yields
    # δD = SG, with G = G_D + ξRG_T + ξUG_lnps and the operator S
    # S = 1 - ξ²(RL + UW) that has to be inverted to obtain δD from the Gs
    for l in 1:lmax+1
        S = LinearAlgebra.I(nlev) - ξ^2*eigenvalues[l]*(R*L + U*W')
        S⁻¹[l,:,:] .= inv(S)
    end
end

function implicit_correction!(  diagn::DiagnosticVariables{NF},
                                progn::PrognosticVariables,
                                model::PrimitiveEquation,
                                ) where NF

    # escape immediately if explicit
    model.parameters.implicit_α == 0 && return nothing   

    @unpack nlev = model.geometry
    @unpack eigenvalues, lmax, mmax = model.spectral_transform
    @unpack Δp_geopot_half, Δp_geopot_full = model.geometry     # = R*Δlnp on half or full levels
    @unpack S⁻¹,R,U,L,W = model.implicit
    ξ = model.implicit.ξ[1]
    
    # MOVE THE IMPLICIT TERMS OF THE TEMPERATURE EQUATION FROM TIME STEP i TO i-1
    # geopotential and linear pressure gradient (divergence equation) are already evaluated at i-1
    # so is the -D̄ term for surface pressure in tendencies!
    for k in 1:nlev
        @unpack div_tend, temp_tend = diagn.layers[k].tendencies
        for r in 1:nlev
            div_old = progn.layers[r].leapfrog[1].div   # divergence at i-1
            div_new = progn.layers[r].leapfrog[2].div   # divergence at i

            # RHS_expl(Vⁱ) + RHS_impl(Vⁱ⁻¹) = RHS(Vⁱ) + RHS_impl(Vⁱ⁻¹ - Vⁱ)
            # for temperature tendency do the latter as its cheaper.
            @. temp_tend += L[k,r]*(div_old-div_new)
            # @. temp_tend += L[k,r]*div_old    # for the former
        end
    end       
    
    # SEMI IMPLICIT CORRECTIONS FOR DIVERGENCE
    @unpack pres_tend = diagn.surface
    @inbounds for k in nlev:-1:1    # loop from bottom layer to top for geopotential calculation
        
        # calculate the combined tendency G = G_D + ξRG_T + ξUG_lnps to solve for divergence δD
        G = diagn.layers[k].dynamics_variables.a        # reuse work arrays
        geopot = diagn.layers[k].dynamics_variables.b
        @unpack div_tend, temp_tend = diagn.layers[k].tendencies

        # 1. the ξRG_T term, vertical integration of geopotential (excl ξ, this is done in 2.)
        # R is not used here as it's cheaper to reuse the geopotential from k+1 than
        # to multiply with the entire upper triangular matrix R which recalculates
        # the geopotential for k from all lower levels k...nlev
        if k == nlev
            @. geopot = Δp_geopot_full[k]*temp_tend        # surface geopotential without orography 
        else
            temp_tend_k1 = diagn.layers[k+1].tendencies.temp_tend   # temp tendency from layer below
            geopot_k1 = diagn.layers[k+1].dynamics_variables.b      # geopotential from layer below
            @. geopot = geopot_k1 + Δp_geopot_half[k+1]*temp_tend_k1 + Δp_geopot_full[k]*temp_tend
        end

        # # alternative way to calculate the geopotential
        # for r in 1:nlev
        #     @unpack temp_tend = diagn.layers[r].tendencies
        #     @. geopot += R[k,r]*temp_tend
        # end

        # 2. the G = G_D + ξRG_T + ξUG_lnps terms using geopot from above 
        lm = 0
        for m in 1:mmax+1               # loops over all columns/order m
            for l in m:lmax+1           # but skips the lmax+2 degree (1-based)
                lm += 1                 # single index lm corresponding to harmonic l,m
                                        # ∇² not part of U so *eigenvalues here
                G[lm] = div_tend[lm] + ξ*eigenvalues[l]*(U[k]*pres_tend[lm] - geopot[lm])      
            end
            lm += 1         # skip last row, LowerTriangularMatrices are of size lmax+2 x mmax+1
        end

        # div_tend is now in G, fill with zeros here so that it can be used as an accumulator
        # in the δD = S⁻¹G calculation below
        fill!(div_tend,0)
    end

    # NOW SOLVE THE δD = S⁻¹G to correct divergence tendency
    @inbounds for k in 1:nlev
        @unpack div_tend = diagn.layers[k].tendencies

        for r in 1:nlev
            G = diagn.layers[r].dynamics_variables.a        # reuse work arrays

            lm = 0
            for m in 1:mmax+1       # loops over all columns/order m
                for l in m:lmax+1   # but skips the lmax+2 degree (1-based)
                    lm += 1         # single index lm corresponding to harmonic l,m
                    div_tend[lm] += S⁻¹[l,k,r]*G[lm]
                end
                lm += 1             # skip last row, LowerTriMatrices are of size lmax+2 x mmax+1
            end
        end
    end

    # SEMI IMPLICIT CORRECTIONS FOR PRESSURE AND TEMPERATURE, insert δD to get δT, δlnpₛ
    @inbounds for k in 1:nlev
        @unpack div_tend, temp_tend = diagn.layers[k].tendencies
        @. pres_tend += ξ*W[k]*div_tend         # δlnpₛ = G_lnpₛ + ξWδD

        for r in 1:nlev
            @unpack div_tend = diagn.layers[r].tendencies
            @. temp_tend += ξ*L[k,r]*div_tend   # δT = G_T + ξLδD
        end
    end
end