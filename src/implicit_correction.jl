initialize_implicit!(Δt::Real,M::BarotropicModel) = nothing

function initialize_implicit!(  dt::Real,                   # time step to update implicit terms with
                                M::ShallowWaterModel{NF}    # update Implicit struct in M
                                ) where NF

    @unpack implicit_α = M.parameters                       # = [0,0.5,1], time step fraction for implicit
    @unpack eigenvalues = M.geospectral.spectral_transform  # =-l*(l+1), degree l of harmonics
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

function implicit_correction!(  diagn::DiagnosticVariables{NF},
                                progn::PrognosticVariables{NF},
                                M::ShallowWaterModel{NF}
                                ) where NF
    
    @unpack div,pres = progn
    @unpack div_tend, pres_tend = diagn.tendencies
    @unpack g∇²,ξg∇²,div_impl = M.implicit
    ξH₀ = M.implicit.ξH₀[1]

    @unpack lmax,mmax = M.geospectral.spectral_transform
    H₀ = M.constants.layer_thickness

    @boundscheck (lmax+1,mmax+1,2,1) == size(div) || throw(BoundsError)
    @boundscheck (lmax+1,mmax+1,1) == size(div_tend) || throw(BoundsError)
    @boundscheck (lmax+1,mmax+1,2) == size(pres) || throw(BoundsError)
    @boundscheck (lmax+1,mmax+1) == size(pres_tend) || throw(BoundsError)

    k = 1       # only one vertical level for shallow water model
    for m in 1:mmax+1
        for l in m:lmax+1
            G_div = div_tend[l,m,k] - g∇²[l]*(pres[l,m,1] - pres[l,m,2])
            G_η   = pres_tend[l,m,k] - H₀*(div[l,m,1,k] - div[l,m,2,k])
            δdiv = (G_div - ξg∇²[l]*G_η)*div_impl[l]
            div_tend[l,m,k] = δdiv
            pres_tend[l,m,k] = G_η - ξH₀*δdiv
        end
    end
end