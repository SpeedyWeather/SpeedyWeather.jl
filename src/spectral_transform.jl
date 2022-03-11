# import FastGaussQuadrature
# import AssociatedLegendrePolynomials
# import FFTW
# import LinearAlgebra

struct SpectralTransform{NF<:AbstractFloat}
    # SPECTRAL RESOLUTION
    lmax::Int       # Maximum degree l=[0,lmax] of spherical harmonics
    mmax::Int       # Maximum order m=[0,l] of spherical harmonics

    # CORRESPONDING GRID SIZE
    nlon::Int               # number of longitudes
    nlat::Int               # number of latitudes
    nlat_half::Int          # nlat on one hemisphere
    nfreq::Int              # number of fourier frequencies (real FFT)
    
    # CORRESPONDING GRID VECTORS
    colat::Vector{NF}       # Gaussian colatitudes (0,π) North to South Pole 
    cos_colat::Vector{NF}   # Cosine of colatitudes
    sin_colat::Vector{NF}   # Sine of colatitudes
    lon_offset::NF          # Offset of first longitude from prime meridian

    # FFT plans
    rfft_plan::FFTW.rFFTWPlan{NF}           # grid to spectral transform
    brfft_plan::FFTW.rFFTWPlan{Complex{NF}} # spectral to grid transform (inverse)

    # LEGENDRE POLYNOMIALS
    recompute_legendre::Bool        # Pre or recompute Legendre polynomials
    Λ::Matrix{NF}                   # Legendre polynomials for one latitude (requires recomputing)
    Λs::Array{NF,3}                 # Legendre polynomials for all latitudes (all precomputed)
    legendre_weights::Vector{NF}    # Legendre weights (extra normalisation of 1/π included)
end

function SpectralTransform( ::Type{NF},                 # Number format NF
                            nlon::Int,                  # Number of longitudes
                            nlat::Int,                  # Number of latitudes
                            trunc::Int,                 # Spectral truncation
                            recompute_legendre::Bool) where NF

    # SPECTRAL RESOLUTION
    lmax = trunc                # Maximum degree l=[0,lmax] of spherical harmonics
    mmax = trunc                # Maximum order m=[0,l] of spherical harmonics
    
    # SIZE OF GRID
    nlat_half = (nlat+1) ÷ 2    # number of latitudes in one hemisphere
    nfreq = nlon÷2 + 1          # number of fourier frequencies (real FFTs)

    # PLAN THE FFTs
    rfft_plan = FFTW.plan_rfft(zeros(NF,nlon))
    brfft_plan = FFTW.plan_brfft(zeros(Complex{NF},nfreq),nlon)

    # GAUSSIAN LATITUDES (0,π) North to South
    colat = π .- acos.(FastGaussQuadrature.gausslegendre(nlat)[1])

    # EQUI-DISTANCE LATITUDES (0,π) North to South
    sin_colat = sin.(colat)
    cos_colat = cos.(colat)
    lon_offset = π/nlon

    # PREALLOCATE LEGENDRE POLYNOMIALS
    Λ = zeros(lmax+1,mmax+1)                # Legendre polynomials for one latitude

    # allocate memory for polynomials at all latitudes or allocate dummy array if precomputed
    # for recomputed only Λ is used, not Λs
    b = ~recompute_legendre                 # true for precomputed
    Λs = zeros(b*lmax + 1, b*mmax + 1, b*(nlat_half-1) + 1)

    if recompute_legendre == false          # then precompute all polynomials
        for ilat in 1:nlat_half             # only one hemisphere due to symmetry
            AssociatedLegendrePolynomials.λlm!(@view(Λs[:,:,ilat]), lmax, mmax, cos_colat[ilat])
        end
    end

    legendre_weights = FastGaussQuadrature.gausslegendre(nlat)[2][1:nlat_half]
    legendre_weights *= π/nlat              # extra normalisation for forward transform included

    # conversion to NF happens here
    SpectralTransform{NF}(  lmax,mmax,
                            nlon,nlat,nlat_half,nfreq,
                            colat,cos_colat,sin_colat,lon_offset,
                            rfft_plan,brfft_plan,
                            recompute_legendre,Λ,Λs,
                            legendre_weights)
end

SpectralTransform(args...) = SpectralTransform(Float64,args...)
SpectralTransform(P::Parameters) = SpectralTransform(P.NF,P.nlon,P.nlat,P.trunc,P.recompute_legendre)

struct GeoSpectral{NF<:AbstractFloat}
    geometry::Geometry{NF}
    spectral::SpectralTransform{NF}
end

function GeoSpectral(P::Parameters)
    G = Geometry(P)
    S = SpectralTransform(P)
    return GeoSpectral{P.NF}(G,S)
end

function get_legendre_polynomials!( Λ::Matrix{NF},
                                    Λs::Array{NF,3},
                                    ilat::Int,
                                    cos_colat::NF,
                                    recompute_legendre::Bool) where NF
    if recompute_legendre
        # Recalculate the (normalized) λ_l^m(cos(colat)) factors of the ass. Legendre polynomials
        lmax,mmax = size(Λ) .- 1
        AssociatedLegendrePolynomials.λlm!(Λ, lmax, mmax, cos_colat)
    else    # copy over precomputed values
        copyto!(Λ,@view(Λs[:,:,ilat]))
    end
end

function gridded!(  map::AbstractMatrix{NF},                    # gridded output
                    alms::AbstractMatrix{Complex{NF}},          # spectral coefficients input
                    S::SpectralTransform{NF}                    # precomputed parameters struct
                    ) where {NF<:AbstractFloat}                 # number format NF

    @unpack lmax, mmax, nlon, nlat, nlat_half = S
    @unpack cos_colat, nfreq, lon_offset = S
    @unpack recompute_legendre, Λ, Λs = S
    @unpack brfft_plan = S

    # 1-based indexing: coefficients l=0:lmax, m=0:l in alms
    @boundscheck (lmax, mmax) == size(alms) .- 1 || throw(BoundsError)
    @boundscheck (nlon, nlat) == size(map) || throw(BoundsError)

    # preallocate work arrays
    gn = zeros(Complex{NF}, nfreq)          # phase factors for northern latitudes
    gs = zeros(Complex{NF}, nfreq)          # phase factors for southern latitudes

    @inbounds for ilat in 1:nlat_half       # loop over northern latitudes only due to symmetry
        ilat_s = nlat - ilat + 1            # southern latitude index

        # Recalculate or use precomputed Legendre polynomials
        get_legendre_polynomials!(Λ,Λs,ilat,cos_colat[ilat],recompute_legendre)

        # inverse Legendre transform by looping over wavenumbers l,m
        for m in 1:mmax+1                   # Σ_{m=0}^{mmax}, but 1-based index
            accn = zero(Complex{NF})        # accumulators for northern/southern hemisphere
            accs = zero(Complex{NF})        

            for l in m:lmax+1                       # Σ_{l=m}^{lmax}, but 1-based index
                term = alms[l,m] * Λ[l,m]
                accn += term
                accs += isodd(l+m) ? -term : term   # flip sign for southern odd wavenumbers
            end
            w = cis((m-1)*lon_offset)               # longitude offset rotation
            gn[m] += accn*w                         # no aliasing here as we are always close
            gs[m] += accs*w                         # to the triangular truncation
        end

        # Fourier transform in zonal direction
        LinearAlgebra.mul!(@view(map[:,ilat]),  brfft_plan,gn)  # Northern latitude
        LinearAlgebra.mul!(@view(map[:,ilat_s]),brfft_plan,gs)  # Southern latitude

        fill!(gn, zero(Complex{NF}))    # set phase factors back to zero
        fill!(gs, zero(Complex{NF}))
    end

    return map
end

function gridded(alms::AbstractMatrix{Complex{NF}}  # spectral coefficients
                ) where NF                          # number format NF

    # check for square matrix of coefficients
    lmax, mmax = size(alms) .- 1
    @boundscheck lmax == mmax || throw(BoundsError)

    # get grid size from spectral resolution
    nlon = roundup_fft(3*lmax+1)    # number of longitudes from triangular truncation
    nlat = nlon÷2                   # number of latitudes
    recompute_legendre = true       # saves memory

    S = SpectralTransform(NF,nlon,nlat,lmax,recompute_legendre)
    return gridded(alms,S)          # now execute the in-place version
end

function gridded(   alms::AbstractMatrix{Complex{NF}},  # spectral coefficients
                    S::SpectralTransform{NF}            # struct for spectral transform parameters
                    ) where NF                          # number format NF

    # check for square matrix of coefficients
    lmax, mmax = size(alms) .- 1
    @boundscheck lmax == mmax || throw(BoundsError)
    @boundscheck lmax == S.lmax || throw(BoundsError)

    output = Matrix{NF}(undef,S.nlon,S.nlat)    # preallocate output
    gridded!(output,alms,S)                     # now execute the in-place version
    return output
end

function spectral!( alms::AbstractMatrix{Complex{NF}},
                    map::AbstractMatrix{NF},
                    S::SpectralTransform{NF}
                    ) where {NF<:AbstractFloat}

    @unpack lmax, mmax, nlon, nlat, nfreq, nlat_half = S
    @unpack cos_colat, sin_colat, lon_offset = S
    @unpack recompute_legendre, Λ, Λs, legendre_weights = S
    @unpack rfft_plan = S
    
    @boundscheck (nlon, nlat) == size(map) || throw(BoundsError)
    @boundscheck (lmax, mmax) == size(alms) .- 1 || throw(BoundsError)

    # preallocate work warrays
    fn = zeros(Complex{NF},nfreq)   # Fourier-transformed northern latitude
    fs = zeros(Complex{NF},nfreq)   # Fourier-transformed southern latitude

    # partial sums are accumulated in alms, force zeros initially.
    fill!(alms,zero(Complex{NF}))   

    @inbounds for ilat in 1:nlat_half   # loop over northern latitudes only due to symmetry
        ilat_s = nlat - ilat + 1        # corresponding southern latitude index

        # Fourier transform in zonal direction
        LinearAlgebra.mul!(fn, rfft_plan, @view(map[:,ilat]))       # Northern latitude
        LinearAlgebra.mul!(fs, rfft_plan, @view(map[:,ilat_s]))     # Southern latitude

        # Legendre transform in meridional direction
        # Recalculate or use precomputed Legendre polynomials
        get_legendre_polynomials!(Λ,Λs,ilat,cos_colat[ilat],recompute_legendre)
        # sin_colat_ΔΩ = sin_colat[ilat] * ΔΩ
        legendre_weight = legendre_weights[ilat]

        for m in 1:mmax+1                                           # Σ_{m=0}^{mmax}, but 1-based index
            # w = sin_colat_ΔΩ * cis((m-1) * -lon_offset)
            w = legendre_weight * cis((m-1) * -lon_offset)
            an = fn[m] * w
            as = fs[m] * w
            for l in m:lmax+1
                c = isodd(l+m) ? an - as : an + as
                alms[l,m] += c * Λ[l,m]
            end
        end
    end

    return alms
end

function spectral(  map::AbstractMatrix{NF} # gridded field nlon x nlat
                    ) where NF              # number format NF

    # check grid is compatible with triangular spectral truncation
    nlon, nlat = size(map)
    @boundscheck nlon == 2nlat || throw(BoundsError)

    # make assumptions about the spectral resolution from triangular truncation
    trunc = ceil(Int,nlon/3-1)
    recompute_legendre = true       # saves memory

    # allocate spectral transform struct
    S = SpectralTransform(NF,nlon,nlat,trunc,recompute_legendre)
    return spectral(map,S)
end

function spectral(  map::AbstractMatrix{NF},    # gridded field nlon x nlat
                    S::SpectralTransform{NF}    # spectral transform struct
                    ) where NF                  # number format NF

    # check grid is compatible with triangular spectral truncation
    nlon, nlat = size(map)
    @boundscheck nlon == 2nlat || throw(BoundsError)
    @boundscheck nlon == S.nlon || throw(BoundsError)

    output = Matrix{Complex{NF}}(undef,S.lmax+1,S.mmax+1)
    return spectral!(output,map,S)  # in-place version
end

# function show_leakage(l,m,trunc=31)
#     alms = zeros(ComplexF64,trunc+1,trunc+1)
#     alms[l+1,m+1] = 1

#     map = gridded(alms)
#     alms2 = spectral(map)
#     imshow(log10.(abs.(alms-alms2)))
#     colorbar(label="Absolute error, log10(abs(alms-alms2))")
#     xlabel("order m")
#     ylabel("degree l")
#     title("Spectral leakage from a(l=$l,m=$m) = 1", loc="left")
#     return alms2[l+1,m+1]
# end

# function SpectralTransform(P::Parameters,G::Geometry)

#     @unpack nlon, nlat, nlon_half, nlat_half = G
#     @unpack coslat_NH = G
#     @unpack R_earth, trunc = P

#     # SIZE OF SPECTRAL GRID
#     mx = trunc+1
#     nx = trunc+1

#     # PLAN THE FFTs
#     # rfft_plan = plan_rfft(rand(NF,nlon))
#     # irfft_plan = plan_irfft(rand(Complex{NF},nlon_half+1),nlon)

#     # LEGENDRE WEIGHTS from pole to equator (=first half or array)
#     leg_weight = FastGaussQuadrature.gausslegendre(nlat)[2][1:nlat_half]

#     # Spectral packing of speedy is m',n', where m' = m, n' = m+n with m,n being the
#     # conventional wavenumbers. Due to Julia's 1-based indexing subtract two, as
#     # the Legendre polynomials start with
#     nsh2 = zeros(Int, nx)
#     for n in 1:nx
#         for m in 1:mx
#             if m + n - 2 <= trunc + 1
#                 nsh2[n] = nsh2[n] + 1
#             end
#         end
#     end

#     # Epsilon-factors for the recurrence relation of the associated Legendre polynomials.
#     ε,ε⁻¹ = ε_recurrence(mx,nx)

#     # Generate associated Legendre polynomials
#     # get_legendre_poly computes the polynomials at a particular latitiude
#     leg_poly = zeros(mx, nx, nlat_half)
#     for j in 1:nlat_half
#         leg_poly[:,:,j] = legendre_polynomials(j,ε,ε⁻¹,mx,nx,G)
#     end

#     # leg_poly = zeros(mx, nx, nlat_half)
#     # for j in 1:nlat_half
#     #     leg_poly[:,:,j] = AssociatedLegendrePolynomials.λlm(0:mx-1,0:nx-1,coslat_NH[j])
#     # end

#     # LAPLACIANS for harmonic & biharmonic diffusion
#     ∇²,∇⁻²,∇⁴ = Laplacians(mx,nx,R_earth)

#     gradx   = zeros(mx)          #TODO what's this?
#     uvdx    = zeros(mx, nx)
#     uvdym   = zeros(mx, nx)
#     uvdyp   = zeros(mx, nx)
#     gradym  = zeros(mx, nx)
#     gradyp  = zeros(mx, nx)
#     vddym   = zeros(mx, nx)
#     vddyp   = zeros(mx, nx)

#     for m in 1:mx
#         for n in 1:nx
#             m1  = m - 1
#             m2  = m1 + 1
#             el1 = m + n - 2
#             if n == 1
#                 gradx[m]   = m1/R_earth
#                 uvdx[m,1]  = -R_earth/(m1 + 1)
#                 uvdym[m,1] = 0.0
#                 vddym[m,1] = 0.0
#             else
#                 uvdx[m,n]   = -R_earth*m1/(el1*(el1 + 1.0))
#                 gradym[m,n] = (el1 - 1.0)*ε[m2,n]/R_earth
#                 uvdym[m,n]  = -R_earth*ε[m2,n]/el1
#                 vddym[m,n]  = (el1 + 1.0)*ε[m2,n]/R_earth
#             end
#             gradyp[m,n] = (el1 + 2.0)*ε[m2,n+1]/R_earth
#             uvdyp[m,n]  = -R_earth*ε[m2,n+1]/(el1 + 1.0)
#             vddyp[m,n]  = el1*ε[m2,n+1]/R_earth
#         end
#     end

#     SpectralTransform{P.NF}(    trunc,mx,nx,
#                                 # rfft_plan,irfft_plan,
#                                 leg_weight,nsh2,leg_poly,∇²,∇⁻²,∇⁴,
#                                 gradx,uvdx,uvdym,uvdyp,gradym,gradyp,vddym,vddyp)
# end

# """
# Laplacian operator in spectral space via element-wise matrix-matrix multiplication.
# """
# function ∇²(A::Array{Complex{NF},2},
#             G::GeoSpectral{NF}) where {NF<:AbstractFloat}
#     return -G.spectral.∇².*A
# end

# """
# In-place version of ∇².
# """
# function ∇²!(   Out::Array{Complex{NF},2},
#                 In::Array{Complex{NF},2},
#                 G::GeoSpectral{NF}) where {NF<:AbstractFloat}

#     @unpack ∇² = G.spectral

#     mx,nx = size(In)
#     @boundscheck (mx,nx) == size(Out) || throw(BoundsError())
#     @boundscheck (mx,nx) == size(∇²) || throw(BoundsError())

#     for n in 1:nx
#         for m in 1:mx
#             @inbounds Out[m,n] = -∇²[m,n]*In[m,n]
#         end
#     end
# end

# """
# Inverse Laplacian in spectral space via element-wise matrix-matrix multiplication.
# """
# function ∇⁻²(   A::Array{Complex{NF},2},
#                 G::GeoSpectral{NF}) where {NF<:AbstractFloat}
#     return -G.spectral.∇⁻².*A
# end

# """
# Transform a spectral array into grid-point space.
# """
# function gridded(  input::Array{Complex{NF},2},
#                    G::GeoSpectral{NF}) where {NF<:AbstractFloat}
#     return fourier_inverse(legendre_inverse(input,G),G)
# end

# """
# Transform a gridded array into spectral space.
# """
# function spectral(  input::Array{NF,2},
#                     G::GeoSpectral{NF}) where {NF<:AbstractFloat}
#     return legendre(fourier(input,G),G)
# end

# function grad!( ψ::Array{NF,2},
#                 psdx::Array{Complex{NF},2},
#                 psdy::Array{NF,2},
#                 G::GeoSpectral{NF}) where {NF<:AbstractFloat}

#     #TODO boundscheck

#     @unpack trunc, mx, nx = G.spectral
#     @unpack gradx, gradyp, gradym = G.spectral

#     for n in 1:nx
#         psdx[:,n] = gradx.*ψ[:,n]*im
#     end

#     for m in 1:mx
#         psdy[m,1]  =  gradyp[m,1]*ψ[m,2]
#         psdy[m,nx] = -gradym[m,nx]*ψ[m,trunc+1]
#     end

#     for n in 2:trunc+1
#         for m in 1:mx
#             psdy[m,n] = -gradym[m,n]*ψ[m,n-1] + gradyp[m,n]*ψ[m,n+1]
#         end
#     end
# end

# function vds!(  ucosm::Array{NF,2},
#                 vcosm::Array{NF,2},
#                 vorm::Array{NF,2},
#                 divm::Array{NF,2},
#                 G::GeoSpectral{NF}) where {NF<:AbstractFloat}

#     #TODO boundscheck

#     @unpack trunc, mx, nx = G.spectral
#     @unpack gradx, vddym, vddyp = G.spectral

#     #TODO preallocate in a diagnosticvars struct
#     zp = zeros(Complex{NF}, mx,nx)
#     zc = zeros(Complex{NF}, mx,nx)

#     for n in 1:nx
#         zp[:,n] = gradx.*ucosm[:,n]*im
#         zc[:,n] = gradx.*vcosm[:,n]*im
#     end

#     for m in 1:mx
#         #TODO this has an implicit conversion to complex{NF}, issue?
#         vorm[m,1]  = zc[m,1] - vddyp[m,1]*ucosm[m,2]
#         vorm[m,nx] = vddym[m,nx]*ucosm[m,trunc+1]
#         divm[m,1]  = zp[m,1] + vddyp[m,1]*vcosm[m,2]
#         divm[m,nx] = -vddym[m,nx]*vcosm[m,trunc+1]
#     end

#     for n in 2:trunc+1
#         for m in 1:mx
#             #TODO same here
#             vorm[m,n] =  vddym[m,n]*ucosm[m,n-1] - vddyp[m,n]*ucosm[m,n+1] + zc[m,n]
#             divm[m,n] = -vddym[m,n]*vcosm[m,n-1] + vddyp[m,n]*vcosm[m,n+1] + zp[m,n]
#         end
#     end
# end

# function uvspec!(   vorm::Array{NF,2},
#                     divm::Array{NF,2},
#                     ucosm::Array{NF,2},
#                     vcosm::Array{NF,2},
#                     G::GeoSpectral{NF}) where {NF<:AbstractFloat}

#     #TODO boundscheck

#     @unpack trunc, mx, nx = G.spectral
#     @unpack uvdx, uvdyp, uvdym = G.spectral

#     #TODO preallocate elsewhere
#     zp = uvdx.*vorm*im
#     zc = uvdx.*divm*im

#     for m in 1:mx
#         ucosm[m,1]  =  zc[m,1] - uvdyp[m,1]*vorm[m,2]
#         ucosm[m,nx] =  uvdym[m,nx]*vorm[m,trunc+1]
#         vcosm[m,1]  =  zp[m,1] + uvdyp[m,1]*divm[m,2]
#         vcosm[m,nx] = -uvdym[m,nx]*divm[m,trunc+1]
#     end

#     for n in 2:trunc+1
#         for m in 1:mx
#           vcosm[m,n] = -uvdym[m,n]*divm[m,n-1] + uvdyp[m,n]*divm[m,n+1] + zp[m,n]
#           ucosm[m,n] =  uvdym[m,n]*vorm[m,n-1] - uvdyp[m,n]*vorm[m,n+1] + zc[m,n]
#         end
#     end
# end

# function vdspec!(   ug::Array{NF,2},
#                     vg::Array{NF,2},
#                     vorm::Array{NF,2},
#                     divm::Array{NF,2},
#                     kcos::Bool,
#                     G::GeoSpectral{NF}) where {NF<:AbstractFloat}

#     #TODO boundscheck

#     @unpack nlat, nlon, cosgr, cosgr2 = G.geometry

#     #TODO preallocate elsewhere
#     ug1 = zeros(NF, nlon, nlat)
#     vg1 = zeros(NF, nlon, nlat)

#     # either cosgr or cosgr2
#     cosgr = kcos ? cosgr : cosgr2

#     for j in 1:nlat
#         for i in 1:nlon
#             ug1[i,j] = ug[i,j]*cosgr[j]
#             vg1[i,j] = vg[i,j]*cosgr[j]
#         end
#     end

#     #TODO add spectral_trans and geometry as arguments
#     specu = spectral(ug1,G)
#     specv = spectral(vg1,G)
#     vds!(specu, specv, vorm, divm)
# end

"""Truncate spectral field by seting the spectral coefficients of the upper right triangle to zero. """
function spectral_truncation!(  alms::AbstractMatrix{Complex{NF}},  # spectral field to be truncated
                                trunc::Int                          # truncate to total wave number `trunc`
                                ) where {NF<:AbstractFloat}         # number format NF
    
    lmax,mmax = size(alms) .- 1    # degree l, order m or the legendre polynomials

    @inbounds for m in 1:mmax+1
        for l in 1:lmax+1
            if m > l || l > trunc+1
                alms[l,m] = zero(Complex{NF})
            end
        end
    end
    return alms
end

spectral_truncation!(alms::AbstractMatrix) = spectral_truncation!(alms,size(alms)[1])

function spectral_truncation(   alms::AbstractMatrix{Complex{NF}},  # spectral field to be truncated
                                trunc::Int                          # truncate to degree and order trunc
                                ) where {NF<:AbstractFloat}         # number format NF
    
    lmax,mmax = size(alms) .- 1    # degree l, order m or the legendre polynomials
    @boundscheck lmax == mmax || throw(BoundsError)
    trunc > lmax && return spectral_interpolation(alms,trunc)

    alms_trunc = Matrix{Complex{NF}}(undef,trunc+1,trunc+1)
    copyto!(alms_trunc,@view(alms[1:trunc+1,1:trunc+1]))
    spectral_truncation!(alms_trunc,trunc)
    return alms_trunc
end

function spectral_interpolation(    alms::AbstractMatrix{Complex{NF}},  # spectral field to be truncated
                                    trunc::Int                          # truncate to degree and order trunc
                                    ) where {NF<:AbstractFloat}         # number format NF
    
    lmax,mmax = size(alms) .- 1    # degree l, order m or the legendre polynomials
    @boundscheck lmax == mmax || throw(BoundsError)
    trunc <= lmax && return spectral_truncation(alms,trunc)

    alms_interp = zeros(Complex{NF},trunc+1,trunc+1)
    copyto!(@view(alms_interp[1:lmax+1,1:mmax+1]),alms)
    return alms_interp
end

# """Spectral truncation with unpacking Geospectral struct."""
# function spectral_truncation!(  A::AbstractArray{Complex{NF},2},
#                                 G::GeoSpectral{NF}
#                                 ) where NF
#     spectral_truncation!(A,G.spectral.trunc)    # unpack GeoSpectral struct
# end

# """Spectral truncation of a grid-point field with memory allocation."""
# function spectral_truncation(   input::AbstractArray{NF,2},
#                                 G::GeoSpectral{NF}
#                                 ) where NF
#     input_spectral = spectral(input,G)          # allocates memory
#     spectral_truncation!(input_spectral,G)      # in-place truncation

#     # allocates memory to return spectrally truncated gridded field
#     return gridded(input_spectral, G)       
# end

# """In-place version of spectral trunction of a grid-point field."""
# function spectral_truncation!(  input::AbstractArray{NF,2},
#                                 input_spectral::AbstractArray{Complex{NF},2},
#                                 G::GeoSpectral{NF}
#                                 ) where NF
#     spectral!(input_spectral,input,G)       # in-place spectral transform from input to input_spectral
#     spectral_truncation!(input_spectral,G)  # in-place truncation
#     gridded!(input,input_spectral,G)        # in-place backtransform
# end

# function spectral_truncation!(  input::Array{NF,2},
#                                 G::GeoSpectral{NF}
#                                 ) where NF
#     spectral_truncation!(input,spectral(input,G),G)
# end


# """
# Computes Laplacian matrix-operators (for element-wise multiplication).
# Laplacian, Laplacian squared and inverse Laplacian in spectral space
# correspond to a multiplication with the total wavenumber:

#     ∇²_n^m = N*(N+1)/R_earth^2

# with N = m+n-2 the total wavenumber -2 due to 1-based indexing. The biharmonic
# operator is ∇⁴ = (∇²)², the inverse Laplacian is ∇⁻² = 1 ./ ∇².
# """
# function Laplacians(mx::Integer,nx::Integer,R_earth::Real)
#     ∇²   = zeros(mx, nx)
#     ∇⁻²  = zeros(mx, nx)
#     ∇⁴   = zeros(mx, nx)

#     for n in 1:nx
#         for m in 1:mx
#             # total wavenumber is m+n, -2 due to Julia's 1-based indexing
#             N = m+n-2
#             ∇²[m,n] = N*(N+1)/R_earth^2
#             ∇⁴[m,n] = ∇²[m,n]^2
#         end
#     end

#     # inverse Laplacian, the first coefficient being zero corresponds
#     # to some (?, TODO) boundary conditions
#     ∇⁻²[1,1] = 0.0                  # don't divide by 0
#     ∇⁻²[2:end] = 1 ./ ∇⁻²[2:end]    # all other elements in matrix

#     return ∇²,∇⁻²,∇⁴
# end
