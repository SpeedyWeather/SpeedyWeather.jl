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

    @unpack NF, implicit_α = model.parameters           # = [0,0.5,1], time step fraction for implicit
    @unpack eigenvalues = model.spectral_transform      # = -l*(l+1), degree l of harmonics
    @unpack ξH₀,ξg∇²,div_impl = model.implicit          # pull precomputed arrays to be updated
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

    @inbounds for i in eachindex(ξg∇²,div_impl,eigenvalues)
        ξg∇²[i] = gravity*eigenvalues[i]        # update ξ∇² with new ξ
        div_impl[i] = inv(1 - ξH₀[1]*ξg∇²[i])   # update 1/(1-ξ²gH₀∇²) with new ξ
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
initialize_implicit!(::Real,::PrimitiveEquation) = nothing

"""
    I = Implicit(P::Parameters{<:PrimitiveEquation})

Zero generator function for an `ImplicitPrimitiveEq` struct, which holds precomputed arrays for
the implicit correction in the primitive equation model. Actual precomputation happens in initialize_implicit!."""
function Implicit(P::Parameters{<:PrimitiveEquation})    # primitive equation only
    
    @unpack NF,trunc,nlev = P

    # initialize with zeros only, actual initialization depends on time step, done in initialize_implicit!
    S⁻¹ = zeros(0,0)
    L = zeros(0,0)
    R = zeros(0,0)
    U = zeros(0)
    W = zeros(0)

    return ImplicitPrimitiveEq(S⁻¹,L,R,U,W)
end