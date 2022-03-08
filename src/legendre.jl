"""
Epsilon-factors for the recurrence relation of the normalized associated
Legendre polynomials.

    ε_n^m = sqrt(((n+m-2)^2 - (m-1)^2)/(4*(n+m-2)^2 - 1))
    ε⁻¹_n^m = 1/ε_n^m   if ε_n^m != 0 else 0.

with m,n being the wavenumbers of the associated Legendre polynomial P_n^m.
Due to the spectral packing in speedy and Julia's 1-based indexing we have
substituted

    m -> m-1
    n -> n+m-2.

compared to the more conventional

    ε_n^m = sqrt( (n^2-m^2) / (4n^2 - 1) )

    From Krishnamurti, Bedi, Hardiker, 2014. Introduction to
    global spectral modelling, Chapter 6.5 Recurrence Relations, Eq. (6.37)
"""
function ε_recurrence(mx::Integer,nx::Integer)
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

    return ε,ε⁻¹
end


"""
Calculate the Legendre polynomials. TODO reference
"""
function legendre_polynomials(  j::Int,
                                ε::AbstractMatrix,
                                ε⁻¹::AbstractMatrix,
                                mx::Int,
                                nx::Int,
                                G::Geometry)

    @unpack coslat_NH, sinlat_NH = G

    alp = zeros(mx+1,nx)

    x = sinlat_NH[j]
    y = coslat_NH[j]

    # Start recursion with n = 1 (m = l) diagonal
    alp[1,1] = sqrt(0.5)
    for m in 2:mx+1
        alp[m,1] = sqrt(0.5*(2m - 1)/(m-1))*y*alp[m-1,1]
    end

    # Continue with other elements
    for m in 1:mx+1
        alp[m,2] = (x*alp[m,1])*ε⁻¹[m,2]
    end

    for n in 3:nx
        for m in 1:mx+1
            alp[m,n] = (x*alp[m,n-1] - ε[m,n-1]*alp[m,n-2])*ε⁻¹[m,n]
        end
    end

    # Zero polynomials with absolute values smaller than 10^(-30)
    # small = 1e-30
    # for n in 1:nx
    #     for m in 1:mx+1
    #         if abs(alp[m,n]) <= small
    #             alp[m,n] = 0.0
    #         end
    #     end
    # end

    # pick off the required polynomials
    return alp[1:mx,1:nx]
end

"""
Computes the inverse Legendre transform.
"""
function legendre_inverse(  input::Array{Complex{NF},2},
                            G::GeoSpectral{NF}) where {NF<:AbstractFloat}

    @unpack leg_weight, nsh2, leg_poly = G.spectral
    @unpack trunc, mx, nx = G.spectral
    @unpack nlat, nlat_half = G.geometry

    # Initialize output array
    output = zeros(Complex{NF}, mx, nlat)

    # Loop over Northern Hemisphere, computing odd and even decomposition of incoming field
    for j in 1:nlat_half
        j1 = nlat + 1 - j

        # Initialise arrays TODO preallocate them in a separate struct
        even = zeros(Complex{NF}, mx)
        odd  = zeros(Complex{NF}, mx)

        # Compute even decomposition
        for n in 1:2:nx
            for m in 1:nsh2[n]
                even[m] = even[m] + input[m,n]*leg_poly[m,n,j]
            end
        end

        # Compute odd decomposition
        for n in 2:2:nx
            for m in 1:nsh2[n]
                odd[m] = odd[m] + input[m,n]*leg_poly[m,n,j]
            end
        end

        # Compute Southern Hemisphere
        output[:,j1] = even + odd

        # Compute Northern Hemisphere
        output[:,j]  = even - odd
    end
    return output
end

"""
Computes the Legendre transform
"""
function legendre(  input::Array{Complex{NF},2},
                    G::GeoSpectral{NF}) where {NF<:AbstractFloat}

    @unpack leg_weight, nsh2, leg_poly = G.spectral
    @unpack trunc, mx, nx = G.spectral
    @unpack nlat, nlat_half = G.geometry

    # Initialise output array
    output = zeros(Complex{NF}, mx, nx)

    even = zeros(Complex{NF}, mx, nlat_half)
    odd  = zeros(Complex{NF}, mx, nlat_half)

    # Loop over Northern Hemisphere, computing odd and even decomposition of
    # incoming field. The Legendre weights (leg_weight) are applied here
    for j in 1:nlat_half
        # Corresponding Southern Hemisphere latitude
        j1 = nlat + 1 - j

        even[:,j] = (input[:,j1] + input[:,j])*leg_weight[j]
        odd[:,j]  = (input[:,j1] - input[:,j])*leg_weight[j]
    end

    # The parity of an associated Legendre polynomial is the same
    # as the parity of n' - m'. n', m' are the actual total wavenumber and zonal
    # wavenumber, n and m are the indices used for SPEEDY's spectral packing.
    # m' = m - 1 and n' = m + n - 2, therefore n' - m' = n - 1

    # Loop over coefficients corresponding to even associated Legendre polynomials
    for n in 1:2:trunc+1
        for m in 1:nsh2[n]
            output[m,n] = dot(leg_poly[m,n,:], even[m,:])
        end
    end

    # Loop over coefficients corresponding to odd associated Legendre polynomials
    for n in 2:2:trunc+1
        for m in 1:nsh2[n]
            output[m,n] = dot(leg_poly[m,n,:], odd[m,:])
        end
    end

    return output
end

import FFTW
import AssociatedLegendrePolynomials
import LinearAlgebra

function synthesize(alms::AbstractMatrix{Complex{NF}},          # spectral coefficients
                    nlon::Integer,                              # number of longitudes
                    nlat::Integer,                              # number of latitudes
                    brfft_plan::FFTW.rFFTWPlan{Complex{NF}}     # FFT plan
                    ) where {NF<:AbstractFloat}                 # number format NF

    
    lmax, mmax = size(alms) .- 1        # 1-based indexing: coefficient in alms from 0:lmax, 0:mmax

    colat = NF(π)/nlat/2:NF(π)/nlat:π   # colatitudes from north to south pole
    lon_offset = NF(π/nlon)             # offset of the first longitude from prime meridian
    nfreq = nlon÷2 + 1                  # length real-only half of FFT axis; see `rfft()`
    nlat_nh = (nlat+1) ÷ 2              # number of latitudes in northern hemisphere

    # preallocate
    Λ = Matrix{NF}(undef, size(alms)...)    # legendre polynomials (1 latitude only)
    gn = zeros(Complex{NF}, nfreq)          # phase factors for northern latitudes
    gs = zeros(Complex{NF}, nfreq)          # phase factors for southern latitudes
    map = zeros(NF, nlon, nlat)             # output array

    @inbounds for ilat in 1:nlat_nh         # loop over northern latitudes only due to symmetry
        ilat_s = nlat - ilat + 1            # southern latitude index

        # Recalculate the (normalized) λ_l^m(cos(colat)) factors of the ass. Legendre polynomials
        AssociatedLegendrePolynomials.λlm!(Λ, lmax, mmax, cos(colat[ilat]))

        # inverse Legendre transform by looping over wavenumbers l,m
        for m in 1:mmax+1                   # Σ_{m=0}^{mmax}, but 1-based index
            accn = zero(Complex{NF})        # accumulators for northern/southern hemisphere
            accs = zero(Complex{NF})        

            for l in m:lmax+1                       # Σ_{l=m}^{lmax}, but 1-based index
                term = alms[l,m] * Λ[l,m]
                accn += term
                accs += isodd(l+m) ? -term : term   # flip sign for southern odd wavenumbers
            end
            accn, accs = (accn, accs) .* cis((m-1)*lon_offset)    # longitude offset rotation
            gn[m] += accn
            gs[m] += accs
        end

        map[:,ilat]   = brfft_plan*gn   # Fourier transform in zonal direction
        map[:,ilat_s] = brfft_plan*gs

        fill!(gn, zero(Complex{NF}))    # set phase factors back to zero
        fill!(gs, zero(Complex{NF}))
    end
    return map
end

function synthesize(alms::Matrix{Complex{NF}},      # spectral coefficients
                    nlon::Int,                      # number of longitudes
                    nlat::Int,                      # number of latitudes
                    ) where {NF<:AbstractFloat}     # number format NF

    nfreq = nlon÷2 + 1                              # length real-only half of FFT axis; see `rfft()`
    gn = zeros(Complex{NF}, nfreq)                  # phase factors for northern ring
    brfft_plan = FFTW.plan_brfft(gn,nlon)           # plan the FFT
    return synthesize(alms,nlon,nlat,brfft_plan)    # call with planned FFT 
end

function analyze(   map::AbstractMatrix{NF},
                    lmax::Integer,
                    mmax::Integer,
                    rfft_plan::FFTW.rFFTWPlan{NF}) where {NF<:AbstractFloat}
    
    nlon, nlat = size(map)      # number of longitudes, latitudes
    nlon_half = nlon÷2 + 1      # real-symmetric FFT's Nyquist length (index)
    nlat_nh = (nlat+1) ÷ 2      # number of rings in northern hemisphere

    colat = NF(π)/nlat/2:NF(π)/nlat:π   # colatitudes from north to south pole
    lon0 = NF(π)/nlon                   # offset of the first longitude from prime meridian
    ΔΩ = 2NF(π)^2 / (nlat*nlon)

    alms = zeros(Complex{NF},lmax+1,mmax+1)
    Λ = zeros(NF,lmax+1,mmax+1)
    Λw = AssociatedLegendrePolynomials.Work(
            AssociatedLegendrePolynomials.λlm!, Λ, zero(NF))
    f₁ = zeros(Complex{NF},nlon_half)   # northern ring
    f₂ = zeros(Complex{NF},nlon_half)   # southern ring

    @inbounds for ilat in 1:nlat_nh     # loop over northern latitudes only due to symmetry
        ilat_s = nlat - ilat + 1        # corresponding southern latitude index

        # Fourier transform in zonal direction
        LinearAlgebra.mul!(f₁, rfft_plan, @view(map[:,ilat]))       # Northern latitude
        LinearAlgebra.mul!(f₂, rfft_plan, @view(map[:,ilat_s]))     # Southern latitude

        # Legendre transform in meridional direction
        slat, clat = sincos(colat[ilat])
        AssociatedLegendrePolynomials.unsafe_legendre!(Λw, Λ, lmax, mmax, clat)

        # Σ_{m=0}^{mmax}, but 1-based index
        for m in 1:mmax+1                                                     
            a₁, a₂ = (f₁[m], f₂[m]) .* (slat * ΔΩ * cis((m-1) * -lon0))
            for l in m:lmax+1
                c = isodd(l+m) ? a₁ - a₂ : a₁ + a₂
                alms[l,m] += c * Λ[l,m]
            end
        end
    end
    return alms
end

function analyze(   map::AbstractMatrix{NF},
                    lmax::Integer
                    ) where {NF<:AbstractFloat}

    nlon,_ = size(map)                          # number of longitudes
    rfft_plan = FFTW.plan_rfft(zeros(NF,nlon))
    return analyze(map,lmax,lmax,rfft_plan)
end