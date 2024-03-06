abstract type AbstractImplicit <: AbstractModelComponent end

# BAROTROPIC MODEL (no implicit needed)
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
Base.@kwdef struct ImplicitShallowWater{NF<:AbstractFloat} <: AbstractImplicit

    # DIMENSIONS
    trunc::Int

    "coefficient for semi-implicit computations to filter gravity waves"
    α::Float64 = 1

    # PRECOMPUTED ARRAYS, to be initiliased with initialize!
    H::Base.RefValue{NF} = Ref(zero(NF))    # layer_thickness
    ξH::Base.RefValue{NF} = Ref(zero(NF))   # = 2αΔt*layer_thickness, store in RefValue for mutability
    g∇²::Vector{NF} = zeros(NF, trunc+2)     # = gravity*eigenvalues
    ξg∇²::Vector{NF} = zeros(NF, trunc+2)    # = 2αΔt*gravity*eigenvalues
    S⁻¹::Vector{NF} = zeros(NF, trunc+2)     # = 1 / (1-ξH*ξg∇²), implicit operator
end

"""
$(TYPEDSIGNATURES)
Generator using the resolution from `spectral_grid`."""
function ImplicitShallowWater(spectral_grid::SpectralGrid;kwargs...) 
    (;NF, trunc) = spectral_grid
    return ImplicitShallowWater{NF}(;trunc, kwargs...)
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

    (;α, H, ξH, g∇², ξg∇², S⁻¹) = implicit          # precomputed arrays to be updated                
    (;gravity) = planet                             # gravitational acceleration [m/s²]
    (;layer_thickness) = atmosphere                 # shallow water layer thickness [m]

    # implicit time step between i-1 and i+1
    # α = 0   means the gravity wave terms are evaluated at i-1 (forward)
    # α = 0.5 evaluates at i+1 and i-1 (centered implicit)
    # α = 1   evaluates at i+1 (backward implicit)
    # α ∈ [0.5, 1] are also possible which controls the strength of the gravity wave dampening.
    # α = 0.5 slows gravity waves and prevents them from amplifying
    # α > 0.5 will dampen the gravity waves within days to a few timesteps (α=1)

    ξ = α*dt                    # new implicit timestep ξ = α*dt = 2αΔt (for leapfrog) from input dt
    H[] = layer_thickness      # update H the undisturbed layer thickness without mountains
    ξH[] = ξ*layer_thickness   # update ξ*H with new ξ, in RefValue for mutability

    # loop over degree l of the harmonics (implicit terms are independent of order m)
    @inbounds for l in eachindex(g∇², ξg∇², S⁻¹)
        eigenvalue = -l*(l-1)               # =∇², with without 1/radius², 1-based -l*l(l+1) → -l*(l-1)
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
function implicit_correction!(  diagn::DiagnosticVariablesLayer{NF},
                                progn::PrognosticLayerTimesteps{NF},
                                diagn_surface::SurfaceVariables{NF},
                                progn_surface::PrognosticSurfaceTimesteps{NF},
                                implicit::ImplicitShallowWater) where NF

    (;div_tend) = diagn.tendencies          # divergence tendency
    div_old = progn.timesteps[1].div        # divergence at t
    div_new = progn.timesteps[2].div        # divergence at t+dt
    pres_old = progn_surface.timesteps[1].pres  # pressure/η at t
    pres_new = progn_surface.timesteps[2].pres  # pressure/η at t+dt
    (;pres_tend) = diagn_surface            # tendency of pressure/η

    (;g∇², ξg∇², S⁻¹) = implicit
    H = implicit.H[]              # unpack as it's stored in a RefValue for mutation
    ξH = implicit.ξH[]            # unpack as it's stored in a RefValue for mutation

    lmax, mmax = size(div_tend) .- (2, 1)
    @boundscheck length(S⁻¹) == lmax+2 || throw(BoundsError)
    @boundscheck length(ξg∇²) == lmax+2 || throw(BoundsError)
    @boundscheck length(g∇²) == lmax+2 || throw(BoundsError)

    lm = 0
    @inbounds for m in 1:mmax+1
        for l in m:lmax+1
            lm += 1     # single index lm corresponding to harmonic l, m with a LowerTriangularMatrix
            
            # calculate the G = N(Vⁱ) + NI(Vⁱ⁻¹ - Vⁱ) term.
            # Vⁱ is a prognostic variable at time step i
            # N is the right hand side of ∂V\∂t = N(V)
            # NI is the part of N that's calculated semi-implicitily: N = NE + NI
            G_div = div_tend[lm] - g∇²[l]*(pres_old[lm] - pres_new[lm])
            G_η   = pres_tend[lm] - H*(div_old[lm] - div_new[lm])

            # using the Gs correct the tendencies for semi-implicit time stepping
            div_tend[lm] = S⁻¹[l]*(G_div - ξg∇²[l]*G_η)
            pres_tend[lm] = G_η - ξH*div_tend[lm]
        end
        lm += 1     # loop skips last row
    end
end

export ImplicitPrimitiveEquation

"""
Struct that holds various precomputed arrays for the semi-implicit correction to
prevent gravity waves from amplifying in the primitive equation model.
$(TYPEDFIELDS)"""
Base.@kwdef struct ImplicitPrimitiveEquation{NF<:AbstractFloat} <: AbstractImplicit
    
    # DIMENSIONS
    "spectral resolution"
    trunc::Int

    "number of vertical levels"
    nlev::Int

    # PARAMETERS
    "time-step coefficient: 0=explicit, 0.5=centred implicit, 1=backward implicit"
    α::Float64 = 1

    # PRECOMPUTED ARRAYS, to be initiliased with initialize!
    "vertical temperature profile, obtained from diagn"
    temp_profile::Vector{NF} = zeros(NF, nlev)

    "time step 2α*Δt packed in RefValue for mutability"
    ξ::Base.RefValue{NF} = Ref{NF}(0)       
    
    "divergence: operator for the geopotential calculation"
    R::Matrix{NF} = zeros(NF, nlev, nlev)     
    
    "divergence: the -RdTₖ∇² term excl the eigenvalues from ∇² for divergence"
    U::Vector{NF} = zeros(NF, nlev)
    
    "temperature: operator for the TₖD + κTₖDlnps/Dt term"
    L::Matrix{NF} = zeros(NF, nlev, nlev)

    "pressure: vertical averaging of the -D̄ term in the log surface pres equation"
    W::Vector{NF} = zeros(NF, nlev)
    
    "components to construct L, 1/ 2Δσ"
    L0::Vector{NF} = zeros(NF, nlev)

    "vert advection term in the temperature equation (below+above)"
    L1::Matrix{NF} = zeros(NF, nlev, nlev)

    "factor in front of the div_sum_above term"
    L2::Vector{NF} = zeros(NF, nlev)

    "_sum_above operator itself"
    L3::Matrix{NF} = zeros(NF, nlev, nlev)

    "factor in front of div term in Dlnps/Dt"
    L4::Vector{NF} = zeros(NF, nlev)

    "for every l the matrix to be inverted"
    S::Matrix{NF} = zeros(NF, nlev, nlev)

    "combined inverted operator: S = 1 - ξ²(RL + UW)"
    S⁻¹::Array{NF, 3} = zeros(NF, trunc+1, nlev, nlev)   
end

"""$(TYPEDSIGNATURES)
Generator using the resolution from SpectralGrid."""
function ImplicitPrimitiveEquation(spectral_grid::SpectralGrid, kwargs...) 
    (;NF, trunc, nlev) = spectral_grid
    return ImplicitPrimitiveEquation{NF}(;trunc, nlev, kwargs...)
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

    (;trunc, nlev, α, temp_profile, S, S⁻¹, L, R, U, W, L0, L1, L2, L3, L4) = implicit
    (;σ_levels_full, σ_levels_thick) = geometry
    (;R_dry, κ) = atmosphere
    (;Δp_geopot_half, Δp_geopot_full) = geopotential
    (;σ_lnp_A, σ_lnp_B) = adiabatic_conversion

    for k in 1:nlev    
        # use current vertical temperature profile                                     
        temp_profile[k] = diagn.layers[k].temp_average[]   
        
        # return immediately if temp_profile contains NaRs, model blew up in that case
        if !isfinite(temp_profile[k]) return nothing end
    end

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

    ξ = α*dt                        # dt = 2Δt for leapfrog, but = Δt, Δ/2 in first_timesteps!
    implicit.ξ[] = ξ                # also store in Implicit struct
    
    # DIVERGENCE OPERATORS (called g in Hoskins and Simmons 1975, eq 11 and Appendix 1)
    @inbounds for k in 1:nlev               # vertical geopotential integration as matrix operator
        R[1:k, k] .= -Δp_geopot_full[k]      # otherwise equivalent to geopotential! with zero orography
        R[1:k-1, k] .+= -Δp_geopot_half[k]   # incl the minus but excluding the eigenvalues as with U 
    end
    U .= -R_dry*temp_profile        # the R_d*Tₖ∇² term excl the eigenvalues from ∇² for divergence
    
    # TEMPERATURE OPERATOR (called τ in Hoskins and Simmons 1975, eq 9 and Appendix 1)
    L0 .= 1 ./ 2σ_levels_thick
    L2 .= κ*temp_profile.*σ_lnp_A    # factor in front of the div_sum_above term                       
    L4 .= κ*temp_profile.*σ_lnp_B    # factor in front of div term in Dlnps/Dt

    @inbounds for k in 1:nlev
        Tₖ = temp_profile[k]                    # average temperature at k
        k_above = max(1, k-1)                    # layer index above
        k_below = min(k+1, nlev)                 # layer index below
        ΔT_above = Tₖ - temp_profile[k_above]   # temperature difference to layer above
        ΔT_below = temp_profile[k_below] - Tₖ   # and to layer below
        σₖ = σ_levels_full[k]                   # should be Σ_r=1^k Δσᵣ for model top at >0hPa
        σₖ_above = σ_levels_full[k_above]

        for r in 1:nlev
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
    # δD = SG, with G = G_D + ξRG_T + ξUG_lnps and the operator S
    # S = 1 - ξ²(RL + UW) that has to be inverted to obtain δD from the Gs
    I = LinearAlgebra.I(nlev)
    @inbounds for l in 1:trunc+1
        eigenvalue = -l*(l-1)           # 1-based, -l*(l+1) → -l*(l-1)
        S .= I .- ξ^2*eigenvalue*(R*L .+ U*W')

        # inv(S) but saving memory:
        luS = LinearAlgebra.lu!(S)      # in-place LU decomposition (overwriting S)
        Sinv = L1                       # reuse L1 matrix to store inv(S)
        Sinv .= I                       # use ldiv! so last arg needs to be unity matrix
        LinearAlgebra.ldiv!(luS, Sinv)   # now do S\I = S⁻¹ via LU decomposition
        S⁻¹[l, :, :] .= Sinv              # store in array
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

    # (;Δp_geopot_half, Δp_geopot_full) = model.geometry     # = R*Δlnp on half or full levels
    (;nlev, trunc, S⁻¹, R, U, L, W) = implicit
    ξ = implicit.ξ[]
    
    # MOVE THE IMPLICIT TERMS OF THE TEMPERATURE EQUATION FROM TIME STEP i TO i-1
    # geopotential and linear pressure gradient (divergence equation) are already evaluated at i-1
    # so is the -D̄ term for surface pressure in tendencies!
    @floop for k in 1:nlev
        (;temp_tend) = diagn.layers[k].tendencies       # unpack temp_tend
        for r in 1:nlev
            div_old = progn.layers[r].timesteps[1].div   # divergence at i-1
            div_new = progn.layers[r].timesteps[2].div   # divergence at i

            # RHS_expl(Vⁱ) + RHS_impl(Vⁱ⁻¹) = RHS(Vⁱ) + RHS_impl(Vⁱ⁻¹ - Vⁱ)
            # for temperature tendency do the latter as its cheaper.
            @. temp_tend += L[k, r]*(div_old-div_new)
            # @. temp_tend += L[k, r]*div_old    # for the former
        end
    end       
    
    # SEMI IMPLICIT CORRECTIONS FOR DIVERGENCE
    (;pres_tend) = diagn.surface
    @floop for k in 1:nlev    # loop from bottom layer to top for geopotential calculation
        
        # calculate the combined tendency G = G_D + ξRG_T + ξUG_lnps to solve for divergence δD
        G = diagn.layers[k].dynamics_variables.a        # reuse work arrays, used for combined tendency G
        geopot = diagn.layers[k].dynamics_variables.b   # used for geopotential
        (;div_tend) = diagn.layers[k].tendencies        # unpack div_tend

        # 1. the ξ*R*G_T term, vertical integration of geopotential (excl ξ, this is done in 2.)
        @inbounds for r in k:nlev                       # skip 1:k-1 as integration is surface to k
            (;temp_tend) = diagn.layers[r].tendencies   # unpack temp_tend
            @. geopot += R[k, r]*temp_tend
        end

        # alternative way to calculate the geopotential (not thread-safe because geopot from below is reused)
        # R is not used here as it's cheaper to reuse the geopotential from k+1 than
        # to multiply with the entire upper triangular matrix R which recalculates
        # the geopotential for k from all lower levels k...nlev
        # TODO swap sign here and not in + geopot[lm] further down
        # if k == nlev
        #     @. geopot = Δp_geopot_full[k]*temp_tend        # surface geopotential without orography 
        # else
        #     temp_tend_k1 = diagn.layers[k+1].tendencies.temp_tend   # temp tendency from layer below
        #     geopot_k1 = diagn.layers[k+1].dynamics_variables.b      # geopotential from layer below
        #     @. geopot = geopot_k1 + Δp_geopot_half[k+1]*temp_tend_k1 + Δp_geopot_full[k]*temp_tend
        # end

        # 2. the G = G_D + ξRG_T + ξUG_lnps terms using geopot from above 
        lm = 0
        @inbounds for m in 1:trunc+1    # loops over all columns/order m
            for l in m:trunc+1          # but skips the lmax+2 degree (1-based)
                lm += 1                 # single index lm corresponding to harmonic l, m
                                        # ∇² not part of U so *eigenvalues here
                eigenvalue = -l*(l-1)   # 1-based, -l*(l+1) → -l*(l-1)
                G[lm] = div_tend[lm] + ξ*eigenvalue*(U[k]*pres_tend[lm] + geopot[lm])      
            end
            lm += 1         # skip last row, LowerTriangularMatrices are of size lmax+2 x mmax+1
        end

        # div_tend is now in G, fill with zeros here so that it can be used as an accumulator
        # in the δD = S⁻¹G calculation below
        fill!(div_tend, 0)
    end

    # NOW SOLVE THE δD = S⁻¹G to correct divergence tendency
    @floop for k in 1:nlev
        (;div_tend) = diagn.layers[k].tendencies        # unpack div_tend

        @inbounds for r in 1:nlev
            G = diagn.layers[r].dynamics_variables.a    # reuse work arrays

            lm = 0
            for m in 1:trunc+1       # loops over all columns/order m
                for l in m:trunc+1   # but skips the lmax+2 degree (1-based)
                    lm += 1         # single index lm corresponding to harmonic l, m
                    div_tend[lm] += S⁻¹[l, k, r]*G[lm]
                end
                lm += 1             # skip last row, LowerTriMatrices are of size lmax+2 x mmax+1
            end
        end
    end

    # SEMI IMPLICIT CORRECTIONS FOR PRESSURE AND TEMPERATURE, insert δD to get δT, δlnpₛ
    @floop for k in 1:nlev
        (;temp_tend) = diagn.layers[k].tendencies       # unpack temp_tend
        @inbounds for r in 1:nlev
            (;div_tend) = diagn.layers[r].tendencies    # unpack div_tend
            @. temp_tend += ξ*L[k, r]*div_tend           # δT = G_T + ξLδD
        end
    end

    for k in 1:nlev     # not thread safe
        (;div_tend) = diagn.layers[k].tendencies        # unpack div_tend
        @. pres_tend += ξ*W[k]*div_tend                 # δlnpₛ = G_lnpₛ + ξWδD
    end
end