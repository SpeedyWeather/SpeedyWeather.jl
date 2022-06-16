struct Implicit{NF<:AbstractFloat}
    ξH₀::Vector{NF}
    ξ∇²::Vector{NF}
    div_impl::Vector{NF}
end

function Implicit(  P::Parameters,
                    C::Constants{NF},
                    G::GeoSpectral{NF}
                    ) where NF
    
    @unpack Δt,gravity = C                          # time step Δt (scaled), gravity g [m/s²]
    @unpack implicit_α = P                          # implicit time step fraction between i-1 and i+1
    @unpack eigen_values = G.spectral_transform     # = -(l(l+1)), degree l of harmonics (0-based)
    @unpack layer_thickness = P                     # = H₀, layer thickness at rest without mountains

    # implicit time step between i-1 and i+1
    # α = 0   means the implicit terms (gravity waves) are evaluated at i-1 (forward)
    # α = 0.5 evaluates at i+1 and i-1 (centered implicit)
    # α = 1   evaluates at i+1 (backward implicit)

    ξ = 2implicit_α*Δt      # = [0,2Δt], time step within the [forward,backward] range for implicit terms
    ξ∇² = ξ*eigen_values    # = -ξl(l+1)
    ξH₀ = ξ*layer_thickness # = ξ*H₀

    # (inverse of) implicit denominator for divergence tendency correction (inverse)
    div_impl = 1 ./ (1 - ξH₀*ξ∇²)   # = 1/(1+l(l+1)*ξ²H₀)

    return Implicit{NF}([ξH₀],ξ∇²,div_impl)
end

initialize_implicit!(Δt::Real,M::BarotropicModel{NF}) where NF = nothing

function initialize_implicit!(  Δt::Real,                   # time step to update implicit terms with
                                M::ShallowWaterModel{NF}    # update Implicit struct in M
                                ) where NF

    @unpack implicit_α, layer_thickness = M.parameters
    @unpack eigen_values = M.geospectral.spectral_transform
    @unpack ξH₀,ξ∇²,div_impl = M.implicit

    ξ = convert(NF,2implicit_α*Δt)  # create a new implicit timestep ξ = 2αΔt based on input Δt
    ξH₀[1] = convert(NF,ξ*layer_thickness)

    @inbounds for i in eachindex(ξ∇²,div_impl,eigen_values)
        ξ∇²[i] = ξ*eigenvalues[i]
        div_impl[i] = inv(1 - ξH₀[1]*ξ∇²[i])
    end
end

function implicit_correction!(  diagn::DiagnosticVariables{NF},
                                progn::PrognosticVariables{NF},
                                M::ShallowWaterModel{NF}
                                ) where NF
    
    @unpack div,pres = progn
    @unpack div_tend, pres_tend = diagn.tendencies
    @unpack ξH₀,ξ∇²,div_impl = M.implicit
    @unpack lmax,mmax = M.geospectral.spectral_transform
    @unpack implicit_α, layer_thickness = M.parameters

    @boundscheck (lmax+1,mmax+1) == size(div)[1:2] || throw(BoundsError)
    @boundscheck (lmax+1,mmax+1) == size(div_tend)[1:2] || throw(BoundsError)
    @boundscheck size(div) == size(pres) || throw(BoundsError)
    @boundscheck size(div_tend) == size(pres_tend) || throw(BoundsError)

    k = 1       # only one vertical level for shallow water model
    @inbounds for m in 1:mmax+1
        for l in m:lmax+1
            G_div = div_tend[l,m,k] - ξ∇²[l]*(pres[l,m,1] - pres[l,m,2])
            G_η  = pres_tend[l,m,k] - ξH₀[1]*(div[l,m,1,k] - div[l,m,2,k])
            δdiv = (G_div - ξ∇²[l]*Gη)*div_impl[l]
            div_tend[l,m,k] = δdiv
            pres_tend[l,m,k] = G_η - ξH₀[1]*δdiv
        end
    end
end

# struct Implicit{M1<:AbstractMatrix, M2<:AbstractMatrix, A<:AbstractArray, V<:AbstractVector}
#     # Constants for backwards implicit biharmonic diffusion
#     dmp1::M1
#     dmp1d::M1
#     dmp1s::M1

#     xa::M2
#     xb::M2
#     xc::M2
#     xd::M2
#     xe::M2
#     xf::A
#     xj::A

#     tref::V
#     tref1::V
#     tref2::V
#     tref3::V

#     elz::M1
#     σ_thick_x::V
# end

# Initialize constants for the implicit gravity wave computation.
# It is assumed that that all implicit steps are of length 2*Δt and use
# the forward/backward parameter α. initialize_implicit has to be re-called
# whenever either of these 2.0 parameters is changed. initialize_implicit should
# be called even if the explicit option is chosen for the gravity wave
# terms (the reference state temperature tref is subtracted from some
# terms anyway to reduce roundoff error; also the constants needed for
# the biharmonic diffusion, which is assumed always to be backwards
# implicit, are defined in initialize_implicit)


# function Implicit(T, geometry::Geometry, constants::Constants, params::Parameters,
#                   horizontal_diffusion::HorizontalDiffusion, Δt)
#     @unpack nlev, mx, nx, σ_levels_half, σ_levels_full, σ_thick = geometry
#     @unpack Rₑ, g, akap, R, γ = constants
#     @unpack α = params
#     @unpack dmp, dmpd, dmps = horizontal_diffusion

#     dmp1 = zeros(T, mx, nx)
#     dmp1d = zeros(T, mx, nx)
#     dmp1s = zeros(T, mx, nx)

#     xa = zeros(T, nlev, nlev)
#     xb = zeros(T, nlev, nlev)
#     xc = zeros(T, nlev, nlev)
#     xd = zeros(T, nlev, nlev)
#     xe = zeros(T, nlev, nlev)
#     xf = zeros(T, nlev, nlev, mx+nx+1)
#     xj = zeros(T, nlev, nlev, mx+nx+1)

#     tref = zeros(T, nlev)
#     tref1 = zeros(T, nlev)
#     tref2 = zeros(T, nlev)
#     tref3 = zeros(T, nlev)

#     elz = zeros(T, mx, nx)
#     σ_thick_x = zeros(T, nlev)

#     # 1. Constants for backwards implicit biharmonic diffusion
#     dmp1  = 1.0./(1.0 .+ dmp*Δt)
#     dmp1d = 1.0./(1.0 .+ dmpd*Δt)
#     dmp1s = 1.0./(1.0 .+ dmps*Δt)

#     # 1. Constants for implicit gravity wave computation
#     # reference atmosphere, function of sigma only
#     γ_g = γ/(1000.0*g)

#     tref = 288.0max.(0.2, σ_levels_full).^(R*γ_g)
#     tref1 = R*tref
#     tref2 = akap*tref
#     tref3 = σ_levels_full.*tref

#     # Other constants
#     xi = Δt*α
#     xxi = xi/(Rₑ^2.0)
#     σ_thick_x = xi*σ_thick

#     for n in 1:nx
#         for m in 1:mx
#             elz[m,n] = T(m + n - 2)*T(m + n - 1)*xxi
#         end
#     end

#     #T(K) = TEX(K)+YA(K,K')*D(K') + XA(K,K')*SIG(K')
#     ya = zeros(T, nlev, nlev)
#     for k in 1:nlev
#         for k1 in 1:nlev
#             ya[k,k1] = -akap*tref[k]*σ_thick[k1]
#         end
#     end

#     for k in 2:nlev
#         xa[k,k-1] = 0.5*(akap*tref[k]/σ_levels_full[k] - (tref[k] - tref[k-1])/σ_thick[k])
#     end

#     for k in 1:nlev-1
#         xa[k,k] = 0.5*(akap*tref[k]/σ_levels_full[k] - (tref[k+1] - tref[k])/σ_thick[k])
#     end

#     #sig(k)=xb(k,k')*d(k')
#     dsum = zeros(T, nlev)
#     dsum[1] = σ_thick[1]
#     for k in 2:nlev
#         dsum[k] = dsum[k-1] + σ_thick[k]
#     end

#     for k in 1:nlev-1
#         for k1 in 1:nlev
#             xb[k,k1] = σ_thick[k1]*dsum[k]
#             if k1 <= k
#                 xb[k,k1] = xb[k,k1] - σ_thick[k1]
#             end
#         end
#     end

#     #t(k)=tex(k)+xc(k,k')*d(k')
#     for k in 1:nlev
#         for k1 in 1:nlev
#             xc[k,k1] = ya[k,k1]
#             for k2 in 1:nlev-1
#                 xc[k,k1] = xc[k,k1] + xa[k,k2]*xb[k2,k1]
#             end
#         end
#     end

#     #P(K)=XD(K,K')*T(K')
#     for k in 1:nlev
#         for k1 in k+1:nlev
#             xd[k,k1] = R*log(σ_levels_half[k1+1]/σ_levels_half[k1])
#         end
#     end
#     for k in 1:nlev
#         xd[k,k] = R*log(σ_levels_half[k+1]/σ_levels_full[k])
#     end

#     #P(K)=YE(K)+XE(K,K')*D(K')
#     for k in 1:nlev
#         for k1 in 1:nlev
#             for k2 in 1:nlev
#                 xe[k,k1] = xe[k,k1] + xd[k,k2]*xc[k2,k1]
#             end
#         end
#     end

#     for l in 1:mx+nx+1
#         xxx = (T(l)*T(l+1))/(Rₑ^2.0)
#         for k in 1:nlev
#             for k1 in 1:nlev
#                 xf[k,k1,l] = xi^2.0*xxx*(R*tref[k]*σ_thick[k1] - xe[k,k1])
#             end
#         end
#         for k in 1:nlev
#             xf[k,k,l] = xf[k,k,l] + 1.0
#         end
#     end

#     for l in 1:mx+nx+1
#         xj[:,:,l] = inv(xf[:,:,l])
#     end

#     for k in 1:nlev
#         for k1 in 1:nlev
#             xc[k,k1] = xc[k,k1]*xi
#         end
#     end

#     Implicit(dmp1, dmp1d, dmp1s, xa, xb, xc, xd, xe, xf, xj, tref, tref1, tref2, tref3,
#              elz, σ_thick_x)
# end

# Correct tendencies for implicit gravity wave model
#  Input/output : D_tend  = divergence tendency
#                 Tₐ_tend = temperature tendency
#                 pₛ_tend = tendency of log(surface pressure)

# """
# Compute the tendencies for implicit gravity wave model
# """
# function implicit_terms!(div_tend, t_tend, ps_tend)
#     ye = zeros(Complex{T}, mx, nx, nlev)
#     yf = zeros(Complex{T}, mx, nx, nlev)

#     for k1 in 1:nlev
#         for k in 1:nlev
#             ye[:,:,k] = ye[:,:,k] + xd[k,k1]*Tₐ_tend[:,:,k1]
#         end
#     end

#     for k in 1:nlev
#         ye[:,:,k] = ye[:,:,k] + tref1[k]*pₛ_tend
#     end

#     for k in 1:nlev
#         for m in 1:mx
#             for n in 1:nx
#                 yf[m,n,k] = D_tend[m,n,k] + elz[m,n]*ye[m,n,k]
#             end
#         end
#     end

#     D_tend *= zero

#     for n in 1:nx
#         for m in 1:mx
#             if (m + n - 2) != 0
#                 for k1 in 1:nlev
#                     D_tend[m,n,:] = D_tend[m,n,:] + xj[:,k1,m+n-2]*yf[m,n,k1]
#                 end
#             end
#         end
#     end

#     for k in 1:nlev
#         pₛ_tend = pₛ_tend - D_tend[:,:,k]*σ_thick_x[k]
#     end

#     for k in 1:nlev
#         for k1 in 1:nlev
#             Tₐ_tend[:,:,k] = Tₐ_tend[:,:,k] + xc[k,k1]*D_tend[:,:,k1]
#         end
#     end
# end