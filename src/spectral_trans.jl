struct SpectralTrans{T<:AbstractFloat}
    leg_weight::Array{T,1}      # Legendre weights
    nsh2::Vector{Int}
    leg_poly::Array{T,2}        # Legendre polynomials
    el2::Array{T,2}             # el2 = l*(l+1)/(r^2)
    el2⁻¹::Array{T,2}           # el2⁻¹ = 1/el2
    el4::Array{T,2}             # el2^2.0, for biharmonic diffusion

    # used to filter out "non-triangular" part of rhomboidal truncation
    trfilt::Array{T,2}

    # Quantities required by functions grad, uvspec, and vds
    gradx::Array{T,1}
    uvdx::Array{T,2}
    uvdym::Array{T,2}
    uvdyp::Array{T,2}
    gradym::Array{T,2}
    gradyp::Array{T,2}
    vddym::Array{T,2}
    vddyp::Array{T,2}
end

function SpectralTrans(T, geometry::Geometry, Rₑ)
    @unpack nlat, trunc, mx, nx = geometry

    # First compute Gaussian latitudes and weights at the nlat/2 points from pole to equator
    # wt contains the Gaussian weights for quadratures
    leg_weight = get_leg_weights(T, geometry)

    nsh2 = zeros(Int, nx)
    for n in 1:nx
        nsh2[n] = 0
        for m in 1:mx
            N = m + n - 2
            if N <= trunc + 1
                nsh2[n] = nsh2[n] + 1
            end
        end
    end

    ε   = zeros(T, mx+1,nx+1)
    ε⁻¹ = zeros(T, mx+1,nx+1)
    for m in 1:mx+1
        for n in 1:nx+1
            if n == nx + 1
                ε[m,n] = 0.0
            elseif n == 1 && m == 1
                ε[m,n] = 0.0
            else
                ε[m,n] = √((T(n + m - 2)^2.0 - T(m - 1)^2.0)/
                    (4.0*T(n + m - 2)^2.0 - 1.0))
            end
            if ε[m,n] > 0.0
                ε⁻¹[m,n] = 1.0/ε[m,n]
            end
        end
    end

    # Generate associated Legendre polynomials
    # get_legendre_poly computes the polynomials at a particular latitiude
    leg_poly = zeros(Complex{T}, mx, nx, div(nlat,2))
    for j in 1:div(nlat,2)
        poly = get_legendre_poly(T, geometry, j, ε, ε⁻¹)
        leg_poly[:,:,j] = poly
    end

    el2 = zeros(T, mx, nx)
    el2⁻¹ = zeros(T, mx, nx)
    el4 = zeros(T, mx, nx)
    trfilt = zeros(Complex{T}, mx, nx)

    for n in 1:nx
        for m in 1:mx
            N = m - 2 + n
            el2[m,n] = T(N*(N + 1))/Rₑ^2.0
            el4[m,n] = el2[m,n]^2.0
            if N <= trunc
                trfilt[m,n] = Complex{T}(1.0)
            else
                trfilt[m,n] = Complex{T}(0.0)
            end
        end
    end

    el2⁻¹[1,1] = 0.0
    el2⁻¹[2:mx,:] = 1.0./el2[2:mx,:]
    el2⁻¹[1,2:nx] = 1.0./el2[1,2:nx]

    gradx = zeros(T, mx)
    uvdx = zeros(T, mx, nx)
    uvdym = zeros(T, mx, nx)
    uvdyp = zeros(T, mx, nx)
    gradym = zeros(T, mx, nx)
    gradyp = zeros(T, mx, nx)
    vddym = zeros(T, mx, nx)
    vddyp = zeros(T, mx, nx)

    for m in 1:mx
        for n in 1:nx
            m1 = m - 1
            m2 = m1 + 1
            el1 = T(m - 2 + n)
            if n == 1
                gradx[m]   = T(m1)/Rₑ
                uvdx[m,1]  = -Rₑ/T(m1 + 1)
                uvdym[m,1] = 0.0
                vddym[m,1] = 0.0
            else
                uvdx[m,n]   = -Rₑ*T(m1)/(el1*(el1 + 1.0))
                gradym[m,n] = (el1 - 1.0)*ε[m2,n]/Rₑ
                uvdym[m,n]  = -Rₑ*ε[m2,n]/el1
                vddym[m,n]  = (el1 + 1.0)*ε[m2,n]/Rₑ
            end
            gradyp[m,n] = (el1 + 2.0)*ε[m2,n+1]/Rₑ
            uvdyp[m,n]  = -Rₑ*ε[m2,n+1]/(el1 + 1.0)
            vddyp[m,n]  = el1*ε[m2,n+1]/Rₑ
        end
    end

    SpectralTrans(leg_weight, nsh2, leg_poly, el2, el2⁻¹, el4, trfilt, gradx, uvdx, uvdym, uvdyp,
                  gradym, gradyp, vddym, vddyp)
end

# Laplacian and inverse Laplacian
∇²(spectral_trans::SpectralTrans, field) = -field.*spectral_trans.el2
∇⁻²(spectral_trans::SpectralTrans, field) = -field.*spectral_trans.el2⁻¹

# Spectral transforms
function spec_to_grid(geometry::Geometry, spectral_trans::SpectralTrans, input; scale=false)
    fourier_inv(geometry, legendre_inv(geometry, spectral_trans, input), scale=scale)
end

function grid_to_spec(geometry::Geometry, spectral_trans::SpectralTrans, input)
    legendre_dir(geometry, spectral_trans, fourier_dir(geometry, input))
end

function grad!(geometry::Geometry, spectral_trans::SpectralTrans, ψ, psdx, psdy)
    @unpack trunc, mx, nx = geometry
    @unpack gradx, gradyp, gradym = spectral_trans

    for n in 1:nx
        psdx[:,n] = gradx.*ψ[:,n]*im
    end

    for m in 1:mx
        psdy[m,1]  =  gradyp[m,1]*ψ[m,2]
        psdy[m,nx] = -gradym[m,nx]*ψ[m,trunc+1]
    end

    for n in 2:trunc+1
        for m in 1:mx
            psdy[m,n] = -gradym[m,n]*ψ[m,n-1] + gradyp[m,n]*ψ[m,n+1]
        end
    end
end

function vds!(geometry::Geometry, spectral_trans::SpectralTrans, ucosm, vcosm, vorm, divm)
    @unpack trunc, mx, nx = geometry
    @unpack gradx, vddym, vddyp = spectral_trans

    zp = zeros(Complex, mx,nx)
    zc = zeros(Complex, mx,nx)

    for n in 1:nx
        zp[:,n] = gradx.*ucosm[:,n]*im
        zc[:,n] = gradx.*vcosm[:,n]*im
    end

    for m in 1:mx
        vorm[m,1]  = zc[m,1] - vddyp[m,1]*ucosm[m,2]
        vorm[m,nx] = vddym[m,nx]*ucosm[m,trunc+1]
        divm[m,1]  = zp[m,1] + vddyp[m,1]*vcosm[m,2]
        divm[m,nx] = -vddym[m,nx]*vcosm[m,trunc+1]
    end

    for n in 2:trunc+1
        for m in 1:mx
            vorm[m,n] =  vddym[m,n]*ucosm[m,n-1] - vddyp[m,n]*ucosm[m,n+1] + zc[m,n]
            divm[m,n] = -vddym[m,n]*vcosm[m,n-1] + vddyp[m,n]*vcosm[m,n+1] + zp[m,n]
        end
    end
end

function uvspec!(geometry::Geometry, spectral_trans::SpectralTrans, vorm, divm, ucosm, vcosm )
    @unpack trunc, mx, nx = geometry
    @unpack uvdx, uvdyp, uvdym = spectral_trans

    zp = uvdx.*vorm*Complex{typeof(vorm[1,1]).types[1]}(0.0 + 1.0im)
    zc = uvdx.*divm*Complex{typeof(vorm[1,1]).types[1]}(0.0 + 1.0im)

    for m in 1:mx
        ucosm[m,1]  =  zc[m,1] - uvdyp[m,1]*vorm[m,2]
        ucosm[m,nx] =  uvdym[m,nx]*vorm[m,trunc+1]
        vcosm[m,1]  =  zp[m,1] + uvdyp[m,1]*divm[m,2]
        vcosm[m,nx] = -uvdym[m,nx]*divm[m,trunc+1]
    end

    for n in 2:trunc+1
        for m in 1:mx
          vcosm[m,n] = -uvdym[m,n]*divm[m,n-1] + uvdyp[m,n]*divm[m,n+1] + zp[m,n]
          ucosm[m,n] =  uvdym[m,n]*vorm[m,n-1] - uvdyp[m,n]*vorm[m,n+1] + zc[m,n]
        end
    end
end

function vdspec!(geometry::Geometry, ug, vg, vorm, divm, kcos)
    @unpack nlat, nlon, cosgr, cosgr2 = geometry

    ug1 = zeros(T, nlon, nlat)
    vg1 = zeros(T, nlon, nlat)

    if kcos
        for j in 1:nlat
            for i in 1:nlon
                ug1[i,j] = ug[i,j]*cosgr[j]
                vg1[i,j] = vg[i,j]*cosgr[j]
            end
        end
    else
        for j in 1:nlat
            for i in 1:nlon
                ug1[i,j] = ug[i,j]*cosgr2[j]
                vg1[i,j] = vg[i,j]*cosgr2[j]
            end
        end
    end

    specu = grid_to_spec(ug1)
    specv = grid_to_spec(vg1)
    vds!(specu, specv, vorm, divm)
end

function truncate(trfilt, input)
    input.*trfilt
end
