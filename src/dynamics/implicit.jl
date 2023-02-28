# BAROTROPIC MODEL
"""
    nothing = initialize_implicit!(::Real,::Barotropic)

Just passes, as implicit terms are not used in the barotropic model."""
initialize_implicit!(::Real,::Barotropic) = nothing

# !!! structs are defined in define_implicit.jl !!!

# SHALLOW WATER MODEL
"""
    I = Implicit(P::Parameters{<:ShallowWater})

Zero generator function for an `ImplicitShallowWater` struct, which holds precomputed arrays for
the implicit correction in the shallow water model. Actual precomputation happens in initialize_implicit!."""
function Implicit(P::Parameters{<:ShallowWater})    # shallow water model only
    
    @unpack NF,trunc = P

    # initialize with zeros only, actual initialization depends on time step, done in initialize_implicit!
    ξH₀ = [zero(NF)]          # wrap in vec to make mutable 
    g∇² = zeros(NF,trunc+2)
    ξg∇² = zero(g∇²)
    div_impl = zero(g∇²)

    return ImplicitShallowWater(ξH₀,g∇²,ξg∇²,div_impl)
end

"""
    initialize_implicit!(dt::Real,M::BarotropicModel)

Update the implicit terms in `M` for the shallow water model as they depend on the time step `dt`."""
function initialize_implicit!(  dt::Real,               # time step
                                model::ShallowWater)    # update Implicit struct in M

    @unpack implicit_α = model.parameters               # = [0,0.5,1], time step fraction for implicit
    @unpack eigenvalues = model.spectral_transform      # = -l*(l+1), degree l of harmonics
    @unpack ξH₀,g∇²,ξg∇²,div_impl = model.implicit      # pull precomputed arrays to be updated
    @unpack layer_thickness = model.constants           # shallow water layer thickness [m]
    @unpack gravity = model.constants                   # gravitational acceleration [m/s²]                  

    # implicit time step between i-1 and i+1
    # α = 0   means the gravity wave terms are evaluated at i-1 (forward)
    # α = 0.5 evaluates at i+1 and i-1 (centered implicit)
    # α = 1   evaluates at i+1 (backward implicit)
    # α ∈ [0.5,1] are also possible which controls the strength of the gravity wave dampening.
    # α = 0.5 only prevents the waves from amplifying
    # α > 0.5 will dampen the gravity waves within days to a few timesteps (α=1)

    ξ = implicit_α*dt               # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog) from input dt
    ξH₀[1] = ξ*layer_thickness      # update ξ*H₀ with new ξ, in vec for mutability

    # loop over degree l of the harmonics (implicit terms are independent of order m)
    @inbounds for l in eachindex(g∇²,ξg∇²,div_impl,eigenvalues)
        g∇²[l] = gravity*eigenvalues[l]         # doesn't actually change with dt
        ξg∇²[l] = ξ*g∇²[l]                      # update ξg∇² with new ξ
        div_impl[l] = inv(1 - ξH₀[1]*ξg∇²[l])   # update 1/(1-ξ²gH₀∇²) with new ξ
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

    @unpack g∇²,ξg∇²,div_impl = model.implicit
    ξH₀ = model.implicit.ξH₀[1]                 # unpack as it's stored in a vec for mutation
    H₀ = model.constants.layer_thickness

    @boundscheck size(div_old) == size(div_new) || throw(BoundsError)
    @boundscheck size(pres_old) == size(pres_new) || throw(BoundsError)
    @boundscheck size(div_tend) == size(pres_new) || throw(BoundsError)
    @boundscheck size(pres_tend) == size(div_new) || throw(BoundsError)
    lmax,mmax = size(div_tend) .- (2,1)

    @boundscheck length(div_impl) == lmax+2 || throw(BoundsError)
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
            div_tend[lm] = (G_div - ξg∇²[l]*G_η)*div_impl[l]
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
    ξ = zeros(NF,1)             # time step 2α*dt packed in a vector for mutability
    L = zeros(NF,nlev)          # operator for the +TₖD term (reference temperature profile)
    R = zeros(NF,nlev,nlev)     # operator for the geopotential calculation
    U = zeros(NF,nlev)          # the -RdTₖ∇² term excl the eigenvalues from ∇² for divergence
    W = zeros(NF,nlev)          # vertical averaging of the -D̄ term in the log surface pres equation
    S⁻¹ = zeros(NF,trunc+1,nlev,nlev)   # combined inverted operator: S = 1 + ξ²(RL + UW)
    return ImplicitPrimitiveEq(ξ,L,R,U,W,S⁻¹,)
end

function initialize_implicit!(  dt::Real,
                                model::PrimitiveEquation)

    @unpack S⁻¹,L,R,U,W = model.implicit
    @unpack nlev, σ_levels_thick, temp_ref_profile = model.geometry
    @unpack Δp_geopot_half, Δp_geopot_full = model.geometry     # = R*Δlnp on half or full levels
    @unpack R_dry = model.constants
    @unpack eigenvalues, lmax = model.spectral_transform
    α = model.parameters.implicit_α

    # set up W,L,U,R operators from
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
    # here we absorb the time step ξ directly into the operators R <- ξR etc

    ξ = α*dt    # dt = 2Δt for leapfrog, but = Δt, Δ/2 in first_timesteps!
    @. W = -ξ*σ_levels_thick            # the -D̄ term in the log surface pres equation
    @. L = ξ*temp_ref_profile           # the DTₖ term in the temperature equation

    @. U = ξ*R_dry*temp_ref_profile     # the R_d*Tₖ∇² term excl the eigenvalues from ∇² for divergence

    for k in 1:nlev                     # set up vertical geopotential integration as matrix operator
        R[1:k,k] .= Δp_geopot_full[k]   # but otherwise equivalent to geopotential! with zero orography
        R[1:k-1,k] .+= Δp_geopot_half[k]# but excluding the eigenvalues as with U
    end

    @. R *= ξ                           # include timestep ξ
    model.implicit.ξ[1] = ξ             # also store in Implicit struct

    # solving the equations above for δD yields
    # δD = SG, with G = G_D + ξRG_T + ξUG_lnps and the operator S
    # S = 1 + ξ²(RL + UW) that has to be inverted to obtain δD from the Gs
    S = zero(R)
    for l in 1:lmax+1
        # include (neg) eigenvalues for -∇² here
        S .= LinearAlgebra.I(nlev) .- eigenvalues[l]*(R*LinearAlgebra.Diagonal(L) .+ U*W')
        S⁻¹[l,:,:] .= inv(S)
    end
end

function implicit_correction!(  diagn::DiagnosticVariables{NF},
                                model::PrimitiveEquation,
                                ) where NF

    @unpack nlev = model.geometry
    @unpack eigenvalues, lmax, mmax = model.spectral_transform
    @unpack Δp_geopot_half, Δp_geopot_full = model.geometry     # = R*Δlnp on half or full levels
    @unpack S⁻¹,L,R,U,W = model.implicit
    ξ = model.implicit.ξ[1]

    # SEMI IMPLICIT CORRECTIONS FOR DIVERGENCE
    @unpack pres_tend = diagn.surface
    for k in nlev:-1:1      # loop from bottom layer to top for geopotential calculation
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

        # 2. the G = G_D + ξRG_T + ξUG_lnps terms using geopot from above 
        lm = 0
        @inbounds for m in 1:mmax+1     # loops over all columns/order m
            for l in m:lmax+1           # but skips the lmax+2 degree (1-based)
                lm += 1     # single index lm corresponding to harmonic l,m within a LowerTriangularMatrix
                            # -∇² not part of U so -eigenvalues here
                G[lm] = div_tend[lm] - eigenvalues[l]*(U[k]*pres_tend[lm] + ξ*geopot[lm])    
            end
            lm += 1         # skip last row, LowerTriangularMatrices are of size lmax+2 x mmax+1
        end

        # div_tend is now in G, fill with zeros here so that it can be used as an accumulator
        # in the δD = S⁻¹G calculation
        fill!(div_tend,0)
    end

    # NOW SOLVE THE δD = S⁻¹G to correct divergence tendency
    for k in 1:nlev
        @unpack div_tend, temp_tend = diagn.layers[k].tendencies
        for k2 in 1:nlev
            G = diagn.layers[k2].dynamics_variables.a        # reuse work arrays

            lm = 0
            @inbounds for m in 1:mmax+1     # loops over all columns/order m
                for l in m:lmax+1           # but skips the lmax+2 degree (1-based)
                    lm += 1     # single index lm corresponding to harmonic l,m within a LowerTriangularMatrix
                    div_tend[lm] += S⁻¹[l,k,k2]*G[lm]    
                end
                lm += 1         # skip last row, LowerTriangularMatrices are of size lmax+2 x mmax+1
            end
        end

        # SEMI IMPLICIT CORRECTIONS FOR PRESSURE AND TEMPERATURE
        @. pres_tend += div_tend*W[k]   # δlnpₛ = G_lnpₛ + ξWδD, W <- ξW here
        @. temp_tend += div_tend*L[k]   # δT = G_T + ξLδD, L <- ξL here
    end
end