abstract type AbstractImplicit <: AbstractModelComponent end

# BAROTROPIC MODEL (no implicit needed)
export NoImplicit
struct NoImplicit <: AbstractImplicit end
NoImplicit(SG::SpectralGrid) = NoImplicit()
initialize!(::NoImplicit, args...) = nothing
implicit_correction!(::DiagnosticVariables, ::PrognosticVariables, ::NoImplicit) = nothing

# SHALLOW WATER MODEL
export ImplicitShallowWater

"""
Struct that holds various precomputed arrays for the semi-implicit correction to
prevent gravity waves from amplifying in the shallow water model.
$(TYPEDFIELDS)"""
@kwdef struct ImplicitShallowWater{NF<:AbstractFloat} <: AbstractImplicit

    # DIMENSIONS
    trunc::Int

    "[OPTION] coefficient for semi-implicit computations to filter gravity waves, 0.5 <= α <= 1"
    α::NF = 1

    # PRECOMPUTED ARRAYS, to be initiliased with initialize!
    H::Base.RefValue{NF} = Ref(zero(NF))        # layer_thickness
    ξH::Base.RefValue{NF} = Ref(zero(NF))       # = 2αΔt*layer_thickness, store in RefValue for mutability
    g∇²::Vector{NF} = zeros(NF, trunc+2)        # = gravity*eigenvalues
    ξg∇²::Vector{NF} = zeros(NF, trunc+2)       # = 2αΔt*gravity*eigenvalues
    S⁻¹::Vector{NF} = zeros(NF, trunc+2)        # = 1 / (1-ξH*ξg∇²), implicit operator
end

"""
$(TYPEDSIGNATURES)
Generator using the resolution from `spectral_grid`."""
function ImplicitShallowWater(spectral_grid::SpectralGrid; kwargs...)
    (; NF, trunc) = spectral_grid
    return ImplicitShallowWater{NF}(; trunc, kwargs...)
end

# function barrier to unpack the constants struct for shallow water
function initialize!(I::ImplicitShallowWater, dt::Real, ::DiagnosticVariables, model::ShallowWater)
    initialize!(I, dt, model.planet, model.atmosphere)
end

"""
$(TYPEDSIGNATURES)
Update the implicit terms in `implicit` for the shallow water model as they depend on the time step `dt`."""
function initialize!(   
    implicit::ImplicitShallowWater,
    dt::Real,                   # time step used [s]
    planet::AbstractPlanet,
    atmosphere::AbstractAtmosphere,
)

    (; α, H, ξH, g∇², ξg∇², S⁻¹) = implicit          # precomputed arrays to be updated
    (; gravity) = planet                             # gravitational acceleration [m/s²]
    (; layer_thickness) = atmosphere                 # shallow water layer thickness [m]

    # implicit time step between i-1 and i+1
    # α = 0   means the gravity wave terms are evaluated at i-1 (forward)
    # α = 0.5 evaluates at i+1 and i-1 (centered implicit)
    # α = 1   evaluates at i+1 (backward implicit)
    # α ∈ [0.5, 1] are also possible which controls the strength of the gravity wave dampening.
    # α = 0.5 slows gravity waves and prevents them from amplifying
    # α > 0.5 will dampen the gravity waves within days to a few timesteps (α=1)

    ξ = α*dt                   # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog) from input dt
    H[] = layer_thickness      # update H the undisturbed layer thickness without mountains
    ξH[] = ξ*layer_thickness   # update ξ*H with new ξ, in RefValue for mutability

    # loop over degree l of the harmonics (implicit terms are independent of order m)
    @inbounds for l in eachindex(g∇², ξg∇², S⁻¹)
        eigenvalue = -l*(l-1)               # =∇², with without 1/radius², 1-based -l*(l+1) → -l*(l-1)
        g∇²[l] = gravity*eigenvalue         # doesn't actually change with dt
        ξg∇²[l] = ξ*g∇²[l]                  # update ξg∇² with new ξ
        S⁻¹[l] = inv(1 - ξH[]*ξg∇²[l])      # update 1/(1-ξ²gH∇²) with new ξ
    end
end

"""
$(TYPEDSIGNATURES)
Apply correction to the tendencies in `diagn` to prevent the gravity waves from amplifying.
The correction is implicitly evaluated using the parameter `implicit.α` to switch between
forward, centered implicit or backward evaluation of the gravity wave terms."""
function implicit_correction!(  diagn::DiagnosticVariables,
                                progn::PrognosticVariables,
                                implicit::ImplicitShallowWater)

    (; div_tend, pres_tend) = diagn.tendencies # tendency of divergence and pressure/η
    div_old = progn.div[1]      # divergence at t
    div_new = progn.div[2]      # divergence at t+dt
    pres_old = progn.pres[1]    # pressure/η at t
    pres_new = progn.pres[2]    # pressure/η at t+dt

    (; g∇², ξg∇², S⁻¹) = implicit
    H = implicit.H[]              # unpack as it's stored in a RefValue for mutation
    ξH = implicit.ξH[]            # unpack as it's stored in a RefValue for mutation

    lmax, mmax = matrix_size(div_tend)[1:2] .- (2, 1)
    @boundscheck length(S⁻¹) == lmax+2 || throw(BoundsError)
    @boundscheck length(ξg∇²) == lmax+2 || throw(BoundsError)
    @boundscheck length(g∇²) == lmax+2 || throw(BoundsError)

    for k in eachmatrix(div_tend)
        lm = 0
        for m in 1:mmax+1
            for l in m:lmax+1
                lm += 1     # single index lm corresponding to harmonic l, m with a LowerTriangularMatrix
                
                # calculate the G = N(Vⁱ) + NI(Vⁱ⁻¹ - Vⁱ) term.
                # Vⁱ is a prognostic variable at time step i
                # N is the right hand side of ∂V\∂t = N(V)
                # NI is the part of N that's calculated semi-implicitily: N = NE + NI
                G_div = div_tend[lm, k] - g∇²[l]*(pres_old[lm] - pres_new[lm])
                G_η   = pres_tend[lm] - H*(div_old[lm, k] - div_new[lm, k])

                # using the Gs correct the tendencies for semi-implicit time stepping
                div_tend[lm, k] = S⁻¹[l]*(G_div - ξg∇²[l]*G_η)
                pres_tend[lm] = G_η - ξH*div_tend[lm, k]
            end
            lm += 1     # loop skips last row
        end
    end
end

export ImplicitPrimitiveEquation

"""
Struct that holds various precomputed arrays for the semi-implicit correction to
prevent gravity waves from amplifying in the primitive equation model.
$(TYPEDFIELDS)"""
@kwdef struct ImplicitPrimitiveEquation{NF<:AbstractFloat} <: AbstractImplicit
    
    # DIMENSIONS
    "spectral resolution"
    trunc::Int

    "number of vertical layers"
    nlayers::Int

    # PARAMETERS
    "time-step coefficient: 0=explicit, 0.5=centred implicit, 1=backward implicit"
    α::NF = 1

    # PRECOMPUTED ARRAYS, to be initiliased with initialize!
    "vertical temperature profile, obtained from diagn on first time step"
    temp_profile::Vector{NF} = zeros(NF, nlayers)

    "time step 2α*Δt packed in RefValue for mutability"
    ξ::Base.RefValue{NF} = Ref{NF}(0)       
    
    "divergence: operator for the geopotential calculation"
    R::Matrix{NF} = zeros(NF, nlayers, nlayers)     
    
    "divergence: the -RdTₖ∇² term excl the eigenvalues from ∇² for divergence"
    U::Vector{NF} = zeros(NF, nlayers)
    
    "temperature: operator for the TₖD + κTₖDlnps/Dt term"
    L::Matrix{NF} = zeros(NF, nlayers, nlayers)

    "pressure: vertical averaging of the -D̄ term in the log surface pres equation"
    W::Vector{NF} = zeros(NF, nlayers)
    
    "components to construct L, 1/ 2Δσ"
    L0::Vector{NF} = zeros(NF, nlayers)

    "vert advection term in the temperature equation (below+above)"
    L1::Matrix{NF} = zeros(NF, nlayers, nlayers)

    "factor in front of the div_sum_above term"
    L2::Vector{NF} = zeros(NF, nlayers)

    "_sum_above operator itself"
    L3::Matrix{NF} = zeros(NF, nlayers, nlayers)

    "factor in front of div term in Dlnps/Dt"
    L4::Vector{NF} = zeros(NF, nlayers)

    "for every l the matrix to be inverted"
    S::Matrix{NF} = zeros(NF, nlayers, nlayers)

    "combined inverted operator: S = 1 - ξ²(RL + UW)"
    S⁻¹::Array{NF, 3} = zeros(NF, trunc+1, nlayers, nlayers)   
end

"""$(TYPEDSIGNATURES)
Generator using the resolution from SpectralGrid."""
function ImplicitPrimitiveEquation(spectral_grid::SpectralGrid, kwargs...) 
    (; NF, trunc, nlayers) = spectral_grid
    return ImplicitPrimitiveEquation{NF}(; trunc, nlayers, kwargs...)
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
    diagn::DiagnosticVariables,
    geometry::AbstractGeometry,
    geopotential::AbstractGeopotential,
    atmosphere::AbstractAtmosphere,
    adiabatic_conversion::AbstractAdiabaticConversion,
)

    (; trunc, nlayers, α, temp_profile, S, S⁻¹, L, R, U, W, L0, L1, L2, L3, L4) = implicit
    (; σ_levels_full, σ_levels_thick) = geometry
    (; R_dry, κ) = atmosphere
    (; Δp_geopot_half, Δp_geopot_full) = geopotential
    (; σ_lnp_A, σ_lnp_B) = adiabatic_conversion

    # use current vertical temperature profile                                     
    temp_profile .= diagn.temp_average
        
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

    ξ = α*dt                        # dt = 2Δt for leapfrog, but = Δt, Δ/2 in first_timesteps!
    implicit.ξ[] = ξ                # also store in Implicit struct
    
    # DIVERGENCE OPERATORS (called g in Hoskins and Simmons 1975, eq 11 and Appendix 1)
    @inbounds for k in 1:nlayers                # vertical geopotential integration as matrix operator
        R[1:k, k] .= -Δp_geopot_full[k]         # otherwise equivalent to geopotential! with zero orography
        R[1:k-1, k] .+= -Δp_geopot_half[k]      # incl the minus but excluding the eigenvalues as with U 
    end
    U .= -R_dry*temp_profile        # the R_d*Tₖ∇² term excl the eigenvalues from ∇² for divergence
    
    # TEMPERATURE OPERATOR (called τ in Hoskins and Simmons 1975, eq 9 and Appendix 1)
    L0 .= 1 ./ 2σ_levels_thick
    L2 .= κ*temp_profile.*σ_lnp_A    # factor in front of the div_sum_above term                       
    L4 .= κ*temp_profile.*σ_lnp_B    # factor in front of div term in Dlnps/Dt

    @inbounds for k in 1:nlayers
        Tₖ = temp_profile[k]                    # average temperature at k
        k_above = max(1, k-1)                   # layer index above
        k_below = min(k+1, nlayers)             # layer index below
        ΔT_above = Tₖ - temp_profile[k_above]   # temperature difference to layer above
        ΔT_below = temp_profile[k_below] - Tₖ   # and to layer below
        σₖ = σ_levels_full[k]                   # should be Σ_r=1^k Δσᵣ for model top at >0hPa
        σₖ_above = σ_levels_full[k_above]

        for r in 1:nlayers
            L1[k, r] = ΔT_below*σ_levels_thick[r]*σₖ         # vert advection operator below
            L1[k, r] -= k>=r ? σ_levels_thick[r] : 0

            L1[k, r] += ΔT_above*σ_levels_thick[r]*σₖ_above   # vert advection operator above
            L1[k, r] -= (k-1)>=r ? σ_levels_thick[r] : 0
        end

        # _sum_above operator itself
        L3[1:k, k] .= 0                              # fill upper triangle + diagonal with zeros
        L3[k+1:end, k] .= σ_levels_thick[k]          # vert integration top to k-1
    end

    L .= Diagonal(L0)*L1 .+ Diagonal(L2)*L3 .+ Diagonal(L4)  # combine all operators into L

    # PRESSURE OPERATOR (called πᵣ in Hoskins and Simmons, 1975 Appendix 1)
    W .= -σ_levels_thick                # the -D̄ term in the log surface pres equation

    # solving the equations above for δD yields
    # δD = SG, with G = G_D + ξRG_T + ξUG_lnps and the operator S
    # S = 1 - ξ²(RL + UW) that has to be inverted to obtain δD from the Gs
    I = LinearAlgebra.I(nlayers)
    @inbounds for l in 1:trunc+1
        eigenvalue = -l*(l-1)           # 1-based, -l*(l+1) → -l*(l-1)
        S .= I .- ξ^2*eigenvalue*(R*L .+ U*W')

        # inv(S) but saving memory:
        luS = LinearAlgebra.lu!(S)      # in-place LU decomposition (overwriting S)
        Sinv = L1                       # reuse L1 matrix to store inv(S)
        Sinv .= I                       # use ldiv! so last arg needs to be unity matrix
        LinearAlgebra.ldiv!(luS, Sinv)  # now do S\I = S⁻¹ via LU decomposition
        S⁻¹[l, :, :] .= Sinv            # store in array
    end
end

"""$(TYPEDSIGNATURES)
Apply the implicit corrections to dampen gravity waves in the primitive equation models."""
function implicit_correction!(  
    diagn::DiagnosticVariables,
    implicit::ImplicitPrimitiveEquation,
    progn::PrognosticVariables,
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
    div_old, div_new = progn.div    # divergence at i-1 (old), i (new, i.e. current)

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

    for k in 1:nlayers
        for r in k:nlayers      # skip 1:k-1 as integration is surface to k
            for lm in eachharmonic(temp_tend, div_old, div_new)
                # 1. the ξ*R*G_T term, vertical integration of geopotential (excl ξ, this is done in 2.)
                geopot[lm, k] += R[k, r]*temp_tend[lm, r]
            end
        end

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
            lm += 1         # skip last row, LowerTriangularMatrices are of size lmax+2 x mmax+1
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
                temp_tend[lm, k] += ξ*L[k, r]*div_tend[lm, r]
            end
        end

        for lm in eachharmonic(div_tend, temp_tend)
            # δlnpₛ = G_lnpₛ + ξWδD
            pres_tend[lm] += ξ*W[k]*div_tend[lm, k]
        end
    end
end