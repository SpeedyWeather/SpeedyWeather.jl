"""
    nothing = initialize_implicit!(::Real,::BarotropicModel)

Just passes, as implicit terms are not used in the barotropic model."""
initialize_implicit!(::Real,::BarotropicModel) = nothing

"""
    initialize_implicit!(dt::Real,M::BarotropicModel)

Update the implicit terms in `M` for the shallow water model as they depend on the time step `dt`."""
function initialize_implicit!(  dt::Real,                   # time step to update implicit terms with
                                M::ShallowWaterModel{NF}    # update Implicit struct in M
                                ) where NF

    @unpack implicit_α = M.parameters                       # = [0,0.5,1], time step fraction for implicit
    @unpack eigenvalues = M.spectral_transform              # =-l*(l+1), degree l of harmonics
    @unpack ξH₀,ξg∇²,div_impl = M.implicit                  # pull precomputed arrays to be updated
    @unpack layer_thickness = M.constants                   # shallow water layer thickness [m]
    @unpack radius_earth, gravity = M.constants             # gravitational acceleration [m/s²]                  

    α = convert(NF,implicit_α)                              # to avoid promotions
    ξ = α*convert(NF,dt)                                    # new implicit timestep ξ = α*dt = 2αΔt from input dt
    ξH₀[1] = ξ*layer_thickness                              # update ξ*H₀ with new ξ

    @inbounds for i in eachindex(ξg∇²,div_impl,eigenvalues)
        ξg∇²[i] = ξ*gravity*eigenvalues[i]                  # update precomputed ξ∇² with new ξ
        div_impl[i] = inv(1 - ξH₀[1]*ξg∇²[i])               # update precomputed 1/(1-ξ²gH₀∇²) with new ξ
    end
end

initialize_implicit!(::Real,::PrimitiveEquationModel) = nothing

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
                                M::ShallowWaterModel{NF},
                                ) where NF

    @unpack div_tend = diagn.tendencies     # divergence tendency
    div_old = progn.leapfrog[1].div         # divergence at t
    div_new = progn.leapfrog[2].div         # divergence at t+dt
    pres_old, pres_new = pres.leapfrog      # pressure/η at t,t+dt
    @unpack pres_tend = surface             # tendency of pressure/η

    @unpack g∇²,ξg∇²,div_impl = M.implicit
    ξH₀ = M.implicit.ξH₀[1]                 # unpack as it's stored in a vec for mutation
    H₀ = M.constants.layer_thickness

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
            G_div = div_tend[lm] - g∇²[l]*(pres_old[lm] - pres_new[lm])
            G_η   = pres_tend[lm] - H₀*(div_old[lm] - div_new[lm])
            div_tend[lm] = (G_div - ξg∇²[l]*G_η)*div_impl[l]
            pres_tend[lm] = G_η - ξH₀*div_tend[lm]
        end
        lm += 1     # loop skips last row
    end
end