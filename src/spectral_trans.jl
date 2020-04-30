struct SpectralTrans{T<:AbstractFloat}
    leg_weight::Array{T,1}      # Legendre weights
    nsh2::Vector{Int}           # What's this?
    leg_poly::Array{T,2}        # Legendre polynomials
    el2::Array{T,2}             # el2 = l*(l+1)/(r^2)
    el2⁻¹::Array{T,2}           # el2⁻¹ = 1/el2
    el4::Array{T,2}             # el2^2.0, for biharmonic diffusion

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

function SpectralTrans{T}(constants::Constants{T},geometry::Geometry{T}) where T
    @unpack nlat, trunc, mx, nx = geometry
    @unpack R_earth = constants
    nlat_half = nlat ÷ 2

    # First compute Gaussian latitudes and weights at the nlat/2 points from pole to equator
    # wt contains the Gaussian weights for quadratures
    leg_weight = get_leg_weights(T, geometry)

    #TODO what's this?
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

    #TODO what's this, separate into function
    ε   = zeros(mx+1,nx+1)
    ε⁻¹ = zeros(mx+1,nx+1)
    for m in 1:mx+1
        for n in 1:nx+1
            if n == nx + 1
                ε[m,n] = 0.0
            elseif n == 1 && m == 1
                ε[m,n] = 0.0
            else
                ε[m,n] = sqrt(((n+m-2)^2 - (m-1)^2)/(4*(n+m-2)^2 - 1))
            end
            if ε[m,n] > 0.0
                ε⁻¹[m,n] = 1.0/ε[m,n]
            end
        end
    end

    # Generate associated Legendre polynomials
    # get_legendre_poly computes the polynomials at a particular latitiude
    leg_poly = zeros(Complex{T}, mx, nx, nlat_half)
    for j in 1:nlat_half
        leg_poly[:,:,j] = get_legendre_poly(j, ε, ε⁻¹, geometry)
    end

    el2   = zeros(mx, nx)
    el2⁻¹ = zeros(mx, nx)
    el4   = zeros(mx, nx)

    for n in 1:nx
        for m in 1:mx
            N = m+n-2
            el2[m,n] = N*(N+1)/R_earth^2
            el4[m,n] = el2[m,n]^2
        end
    end

    el2⁻¹[1] = 0.0                      # don't divide by 0 for el2[1,1]
    el2⁻¹[2:end] = 1.0 ./ el2[2:end]

    gradx   = zeros(mx)          #TODO what's this?
    uvdx    = zeros(mx, nx)
    uvdym   = zeros(mx, nx)
    uvdyp   = zeros(mx, nx)
    gradym  = zeros(mx, nx)
    gradyp  = zeros(mx, nx)
    vddym   = zeros(mx, nx)
    vddyp   = zeros(mx, nx)

    for m in 1:mx
        for n in 1:nx
            m1  = m - 1
            m2  = m1 + 1
            el1 = m + n - 2
            if n == 1
                gradx[m]   = m1/R_earth
                uvdx[m,1]  = -R_earth/(m1 + 1)
                uvdym[m,1] = 0.0
                vddym[m,1] = 0.0
            else
                uvdx[m,n]   = -R_earth*m1/(el1*(el1 + 1.0))
                gradym[m,n] = (el1 - 1.0)*ε[m2,n]/R_earth
                uvdym[m,n]  = -R_earth*ε[m2,n]/el1
                vddym[m,n]  = (el1 + 1.0)*ε[m2,n]/R_earth
            end
            gradyp[m,n] = (el1 + 2.0)*ε[m2,n+1]/R_earth
            uvdyp[m,n]  = -R_earth*ε[m2,n+1]/(el1 + 1.0)
            vddyp[m,n]  = el1*ε[m2,n+1]/R_earth
        end
    end

    SpectralTrans{T}(leg_weight, nsh2, leg_poly, el2, el2⁻¹, el4,
                    gradx, uvdx, uvdym, uvdyp, gradym, gradyp, vddym, vddyp)
end

# Laplacian and inverse Laplacian
#TODO this is presumably faster with loops
∇²(field, spectral_trans::SpectralTrans) = -field.*spectral_trans.el2
∇⁻²(field, spectral_trans::SpectralTrans) = -field.*spectral_trans.el2⁻¹

# Spectral transforms
function spec_to_grid(  input::AbstractMatrix,
                        spectral_trans::SpectralTrans,
                        geometry::Geometry)
    return fourier_inv(legendre_inv(input,spectral_trans,geometry),geometry)
end

function grid_to_spec(  input::AbstractMatrix,
                        spectral_trans::SpectralTrans,
                        geometry::Geometry)
    return legendre(fourier(input,geometry),spectral_trans,geometry)
end

function grad!( ψ::Array{T,2},
                psdx::Array{Complex{T},2},
                psdy::Array{T,2},
                spectral_trans::SpectralTrans{T},
                geometry::Geometry{T}) where {T<:AbstractFloat}

    #TODO boundscheck

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

function vds!(  ucosm::Array{T,2},
                vcosm::Array{T,2},
                vorm::Arry{T,2},
                divm::Array{T,2},
                spectral_trans::SpectralTrans{T},
                geometry::Geometry{T}) where {T<:AbstractFloat}

    #TODO boundscheck

    @unpack trunc, mx, nx = geometry
    @unpack gradx, vddym, vddyp = spectral_trans

    #TODO preallocate in a diagnostics struct
    zp = zeros(Complex{T}, mx,nx)
    zc = zeros(Complex{T}, mx,nx)

    for n in 1:nx
        zp[:,n] = gradx.*ucosm[:,n]*im
        zc[:,n] = gradx.*vcosm[:,n]*im
    end

    for m in 1:mx
        #TODO this has an implicit conversion to complex{T}, issue?
        vorm[m,1]  = zc[m,1] - vddyp[m,1]*ucosm[m,2]
        vorm[m,nx] = vddym[m,nx]*ucosm[m,trunc+1]
        divm[m,1]  = zp[m,1] + vddyp[m,1]*vcosm[m,2]
        divm[m,nx] = -vddym[m,nx]*vcosm[m,trunc+1]
    end

    for n in 2:trunc+1
        for m in 1:mx
            #TODO same here
            vorm[m,n] =  vddym[m,n]*ucosm[m,n-1] - vddyp[m,n]*ucosm[m,n+1] + zc[m,n]
            divm[m,n] = -vddym[m,n]*vcosm[m,n-1] + vddyp[m,n]*vcosm[m,n+1] + zp[m,n]
        end
    end
end

function uvspec!(   vorm::Array{T,2},
                    divm::Array{T,2},
                    ucosm::Array{T,2},
                    vcosm::Array{T,2},
                    spectral_trans::SpectralTrans{T},
                    geometry::Geometry{T}) where {T<:AbstractFloat}

    #TODO boundscheck

    @unpack trunc, mx, nx = geometry
    @unpack uvdx, uvdyp, uvdym = spectral_trans

    #TODO preallocate elsewhere
    zp = uvdx.*vorm*im
    zc = uvdx.*divm*im

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

function vdspec!(   ug::Array{T,2},
                    vg::Array{T,2},
                    vorm::Array{T,2},
                    divm::Array{T,2},
                    kcos::Bool,
                    geometry::Geometry{T}) where {T<:AbstractFloat}

    #TODO boundscheck

    @unpack nlat, nlon, cosgr, cosgr2 = geometry

    #TODO preallocate elsewhere
    ug1 = zeros(T, nlon, nlat)
    vg1 = zeros(T, nlon, nlat)

    # either cosgr or cosgr2
    cosgr = kcos ? cosgr : cosgr2

    for j in 1:nlat
        for i in 1:nlon
            ug1[i,j] = ug[i,j]*cosgr[j]
            vg1[i,j] = vg[i,j]*cosgr[j]
        end
    end

    #TODO add spectral_trans and geometry as arguments
    specu = grid_to_spec(ug1)
    specv = grid_to_spec(vg1)
    vds!(specu, specv, vorm, divm)
end

""" Set the spectral coefficients of the lower right triangle
to zero. """
function truncate!(A::AbstractMatrix)
    m,n = size(A)
    o = zero(eltype(A))

    @inbounds for j in 1:n
        for i in 1:m
            if i+j-1 > m    #TODO dependency on trunc? was i+j-2 > trunc
                A[i,j] = o
            end
        end
    end
end
