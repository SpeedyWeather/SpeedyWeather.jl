struct Implicit{M1<:AbstractMatrix, M2<:AbstractMatrix, A<:AbstractArray, V<:AbstractVector}
    # Constants for backwards implicit biharmonic diffusion
    dmp1::M1
    dmp1d::M1
    dmp1s::M1

    xa::M2
    xb::M2
    xc::M2
    xd::M2
    xe::M2
    xf::A
    xj::A

    tref::V
    tref1::V
    tref2::V
    tref3::V

    elz::M1
    σ_thick_x::V
end

# Initialize constants for the implicit gravity wave computation.
# It is assumed that that all implicit steps are of length 2*Δt and use
# the forward/backward parameter α. initialize_implicit has to be re-called
# whenever either of these 2.0 parameters is changed. initialize_implicit should
# be called even if the explicit option is chosen for the gravity wave
# terms (the reference state temperature tref is subtracted from some
# terms anyway to reduce roundoff error; also the constants needed for
# the biharmonic diffusion, which is assumed always to be backwards
# implicit, are defined in initialize_implicit)
function Implicit(T, geometry::Geometry, constants::Constants, params::Params,
                  horizontal_diffusion::HorizontalDiffusion, Δt)
    @unpack nlev, mx, nx, σ_half, σ_full, σ_thick = geometry
    @unpack Rₑ, g, akap, R, γ = constants
    @unpack α = params
    @unpack dmp, dmpd, dmps = horizontal_diffusion

    dmp1 = zeros(T, mx, nx)
    dmp1d = zeros(T, mx, nx)
    dmp1s = zeros(T, mx, nx)

    xa = zeros(T, nlev, nlev)
    xb = zeros(T, nlev, nlev)
    xc = zeros(T, nlev, nlev)
    xd = zeros(T, nlev, nlev)
    xe = zeros(T, nlev, nlev)
    xf = zeros(T, nlev, nlev, mx+nx+1)
    xj = zeros(T, nlev, nlev, mx+nx+1)

    tref = zeros(T, nlev)
    tref1 = zeros(T, nlev)
    tref2 = zeros(T, nlev)
    tref3 = zeros(T, nlev)

    elz = zeros(T, mx, nx)
    σ_thick_x = zeros(T, nlev)

    # 1. Constants for backwards implicit biharmonic diffusion
    dmp1  = 1.0./(1.0 .+ dmp*Δt)
    dmp1d = 1.0./(1.0 .+ dmpd*Δt)
    dmp1s = 1.0./(1.0 .+ dmps*Δt)

    # 1. Constants for implicit gravity wave computation
    # reference atmosphere, function of sigma only
    γ_g = γ/(1000.0*g)

    tref = 288.0max.(0.2, σ_full).^(R*γ_g)
    tref1 = R*tref
    tref2 = akap*tref
    tref3 = σ_full.*tref

    # Other constants
    xi = Δt*α
    xxi = xi/(Rₑ^2.0)
    σ_thick_x = xi*σ_thick

    for n in 1:nx
        for m in 1:mx
            elz[m,n] = T(m + n - 2)*T(m + n - 1)*xxi
        end
    end

    #T(K) = TEX(K)+YA(K,K')*D(K') + XA(K,K')*SIG(K')
    ya = zeros(T, nlev, nlev)
    for k in 1:nlev
        for k1 in 1:nlev
            ya[k,k1] = -akap*tref[k]*σ_thick[k1]
        end
    end

    for k in 2:nlev
        xa[k,k-1] = 0.5*(akap*tref[k]/σ_full[k] - (tref[k] - tref[k-1])/σ_thick[k])
    end

    for k in 1:nlev-1
        xa[k,k] = 0.5*(akap*tref[k]/σ_full[k] - (tref[k+1] - tref[k])/σ_thick[k])
    end

    #sig(k)=xb(k,k')*d(k')
    dsum = zeros(T, nlev)
    dsum[1] = σ_thick[1]
    for k in 2:nlev
        dsum[k] = dsum[k-1] + σ_thick[k]
    end

    for k in 1:nlev-1
        for k1 in 1:nlev
            xb[k,k1] = σ_thick[k1]*dsum[k]
            if k1 <= k
                xb[k,k1] = xb[k,k1] - σ_thick[k1]
            end
        end
    end

    #t(k)=tex(k)+xc(k,k')*d(k')
    for k in 1:nlev
        for k1 in 1:nlev
            xc[k,k1] = ya[k,k1]
            for k2 in 1:nlev-1
                xc[k,k1] = xc[k,k1] + xa[k,k2]*xb[k2,k1]
            end
        end
    end

    #P(K)=XD(K,K')*T(K')
    for k in 1:nlev
        for k1 in k+1:nlev
            xd[k,k1] = R*log(σ_half[k1+1]/σ_half[k1])
        end
    end
    for k in 1:nlev
        xd[k,k] = R*log(σ_half[k+1]/σ_full[k])
    end

    #P(K)=YE(K)+XE(K,K')*D(K')
    for k in 1:nlev
        for k1 in 1:nlev
            for k2 in 1:nlev
                xe[k,k1] = xe[k,k1] + xd[k,k2]*xc[k2,k1]
            end
        end
    end

    for l in 1:mx+nx+1
        xxx = (T(l)*T(l+1))/(Rₑ^2.0)
        for k in 1:nlev
            for k1 in 1:nlev
                xf[k,k1,l] = xi^2.0*xxx*(R*tref[k]*σ_thick[k1] - xe[k,k1])
            end
        end
        for k in 1:nlev
            xf[k,k,l] = xf[k,k,l] + 1.0
        end
    end

    for l in 1:mx+nx+1
        xj[:,:,l] = inv(xf[:,:,l])
    end

    for k in 1:nlev
        for k1 in 1:nlev
            xc[k,k1] = xc[k,k1]*xi
        end
    end

    Implicit(dmp1, dmp1d, dmp1s, xa, xb, xc, xd, xe, xf, xj, tref, tref1, tref2, tref3,
             elz, σ_thick_x)
end

# Correct tendencies for implicit gravity wave model
#  Input/output : D_tend  = divergence tendency
#                 Tₐ_tend = temperature tendency
#                 pₛ_tend = tendency of log(surface pressure)
function implicit_terms!(D_tend, Tₐ_tend, pₛ_tend)
    ye = zeros(Complex{T}, mx, nx, nlev)
    yf = zeros(Complex{T}, mx, nx, nlev)

    for k1 in 1:nlev
        for k in 1:nlev
            ye[:,:,k] = ye[:,:,k] + xd[k,k1]*Tₐ_tend[:,:,k1]
        end
    end

    for k in 1:nlev
        ye[:,:,k] = ye[:,:,k] + tref1[k]*pₛ_tend
    end

    for k in 1:nlev
        for m in 1:mx
            for n in 1:nx
                yf[m,n,k] = D_tend[m,n,k] + elz[m,n]*ye[m,n,k]
            end
        end
    end

    D_tend *= zero

    for n in 1:nx
        for m in 1:mx
            if (m + n - 2) != 0
                for k1 in 1:nlev
                    D_tend[m,n,:] = D_tend[m,n,:] + xj[:,k1,m+n-2]*yf[m,n,k1]
                end
            end
        end
    end

    for k in 1:nlev
        pₛ_tend = pₛ_tend - D_tend[:,:,k]*σ_thick_x[k]
    end

    for k in 1:nlev
        for k1 in 1:nlev
            Tₐ_tend[:,:,k] = Tₐ_tend[:,:,k] + xc[k,k1]*D_tend[:,:,k1]
        end
    end
end
