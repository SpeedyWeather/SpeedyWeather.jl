initialize_implicit!(::Real,::BarotropicModel) = nothing

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
    ξH₀ = M.implicit.ξH₀[1]
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
    for m in 1:mmax+1
        for l in m:lmax+1
            lm += 1
            G_div = div_tend[lm] - g∇²[l]*(pres_old[lm] - pres_new[lm])
            G_η   = pres_tend[lm] - H₀*(div_old[lm] - div_new[lm])
            δdiv = (G_div - ξg∇²[l]*G_η)*div_impl[l]
            div_tend[lm] = δdiv
            pres_tend[lm] = G_η - ξH₀*δdiv
        end
        lm += 1     # loop skips last row
    end
end