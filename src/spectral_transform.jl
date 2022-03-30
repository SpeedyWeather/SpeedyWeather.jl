"""SpectralTransform struct that contains all parameters and preallocated arrays
for the spectral transform."""
struct SpectralTransform{NF<:AbstractFloat}
    # SPECTRAL RESOLUTION
    lmax::Int       # Maximum degree l=[0,lmax] of spherical harmonics
    mmax::Int       # Maximum order m=[0,l] of spherical harmonics
    nfreq::Int      # Number of fourier frequencies (real FFT)

    # CORRESPONDING GRID SIZE
    nlon::Int               # Number of longitudes
    nlat::Int               # Number of latitudes
    nlat_half::Int          # nlat on one hemisphere
    
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
    legendre_weights::Vector{NF}    # Legendre weights (extra normalisation of π/nlat included)
end

"""
    S = SpectralTransform(NF,nlon,nlat,trunc,recompute_legendre)

Generator function for a SpectralTransform struct. From number of longitudes `nlon`,
number of latitudes `nlat` and spectral truncation `trunc` this function sets up the
spectral resolution, plans the Fourier transforms, retrieves the Gaussian colatitudes,
and preallocates the Legendre polynomials (if recompute_legendre == false) and legendre
weights.
"""
function SpectralTransform( ::Type{NF},                 # Number format NF
                            nlon::Int,                  # Number of longitudes
                            nlat::Int,                  # Number of latitudes
                            trunc::Int,                 # Spectral truncation
                            recompute_legendre::Bool) where NF

    # SPECTRAL RESOLUTION
    lmax = trunc                # Maximum degree l=[0,lmax] of spherical harmonics
    mmax = trunc                # Maximum order m=[0,l] of spherical harmonics
    nfreq = nlon÷2 + 1          # Number of fourier frequencies (real FFTs)

    # SIZE OF GRID
    nlat_half = (nlat+1) ÷ 2    # number of latitudes in one hemisphere
    
    # PLAN THE FFTs
    rfft_plan = FFTW.plan_rfft(zeros(NF,nlon))
    brfft_plan = FFTW.plan_brfft(zeros(Complex{NF},nfreq),nlon)

    # GAUSSIAN COLATITUDES (0,π) North to South
    nodes = FastGaussQuadrature.gausslegendre(nlat)[1]  # zeros of the Legendre polynomial
    colat = π .- acos.(nodes)                           # corresponding colatitudes
    sin_colat = sin.(colat)                             
    cos_colat = cos.(colat)
    lon_offset = π/nlon                                 # offset of first longitude from prime meridian

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
    SpectralTransform{NF}(  lmax,mmax,nfreq,
                            nlon,nlat,nlat_half,
                            colat,cos_colat,sin_colat,lon_offset,
                            rfft_plan,brfft_plan,
                            recompute_legendre,Λ,Λs,
                            legendre_weights)
end

"""Generator function for a SpectralTransform struct in case the number format is not provided.
Use Float64 as default."""
SpectralTransform(args...) = SpectralTransform(Float64,args...)

"Generator function for a SpectralTransform struct pulling in parameters from a Parameters struct."
SpectralTransform(P::Parameters) = SpectralTransform(P.NF,P.nlon,P.nlat,P.trunc,P.recompute_legendre)

"""
    G = Geospectral{NF}(geometry,spectraltrans)

Struct that holds both a Geometry struct `geometry` and a SpectralTransform struct `spectraltrans`."""
struct GeoSpectral{NF<:AbstractFloat}
    geometry::Geometry{NF}
    spectral::SpectralTransform{NF}
end

"""
    G = GeoSpectral(P)

Generator function for a GeoSpectral struct based on a Parameters struct `P`."""
function GeoSpectral(P::Parameters)
    G = Geometry(P)
    S = SpectralTransform(P)
    return GeoSpectral{P.NF}(G,S)
end

"""
    get_legendre_polynomials!(Λ,Λs,ilat,cos_colat,recompute_legendre)

Base on `recompute_legendre` (true/false) this function either updates the Legendre polynomials
`Λ` for a given latitude `ilat, cos_colat` by recomputation (`recompute_legendre == true`), or
`Λ` is changed by copying over the polynomials from precomputed `Λs`. Recomputation takes usually
longer, but precomputation requires a large amount of memory for high resolution."""
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

"""
    gridded!(map,alms,S)

Backward or inverse spectral transform (spectral to grid space) from coefficients `alms` and SpectralTransform
struct `S` into the preallocated output `map`. Uses a planned inverse Fast Fourier Transform for efficiency in the
zonal direction and a Legendre Transform in the meridional direction exploiting symmetries for effciency.
Either recomputes the Legendre polynomials to save memory, or uses precomputed polynomials from `S` depending on
`S.recompute_legendre`."""
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
                term = alms[l,m] * Λ[l,m]           # Legendre polynomials in Λ
                accn += term
                accs += isodd(l+m) ? -term : term   # flip sign for southern odd wavenumbers
            end
            w = cis((m-1)*lon_offset)               # longitude offset rotation
            gn[m] += accn*w                         # no aliasing here as we are always close
            gs[m] += accs*w                         # to the triangular truncation
        end

        # Inverse Fourier transform in zonal direction
        LinearAlgebra.mul!(@view(map[:,ilat]),  brfft_plan,gn)  # Northern latitude
        LinearAlgebra.mul!(@view(map[:,ilat_s]),brfft_plan,gs)  # Southern latitude

        fill!(gn, zero(Complex{NF}))    # set phase factors back to zero
        fill!(gs, zero(Complex{NF}))
    end

    return map
end

"""
    map = gridded(alms)

Backward or inverse spectral transform (spectral to grid space) from coefficients `alms`. Based on the size
of `alms` the corresponding grid space resolution is retrieved based on triangular truncation and a 
SpectralTransform struct `S` is allocated to execute `gridded(alms,S)`."""
function gridded(alms::AbstractMatrix{Complex{NF}}  # spectral coefficients
                ) where NF                          # number format NF

    lmax, mmax = size(alms) .- 1                    # l = 0,lmax hence -1
    @boundscheck lmax == mmax || throw(BoundsError) # check for square matrix of coefficients

    # get grid size from spectral resolution via triangular_truncation
    nlon, nlat = triangular_truncation(lmax)        # number of longitudes, number of latitudes
    recompute_legendre = true                       # saves memory as Legendre precomputation is
                                                    # unnecessary as S is not stored

    S = SpectralTransform(NF,nlon,nlat,lmax,recompute_legendre)
    return gridded(alms,S)          # now execute the in-place version
end

"""
    map = gridded(alms,S)

Backward or inverse spectral transform (spectral to grid space) from coefficients `alms` and the 
SpectralTransform struct `S`. Allocates the output `map` with Gaussian latitudes and executes
`gridded!(map,alms,S)`."""
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

"""
    spectral!(alms,map,S)

Forward spectral transform (grid to spectral space) from the gridded field `map` on a regular Gaussian
grid (with Gaussian latitudes). Uses a planned real-valued Fast Fourier Transform in the zonal direction,
and a Legendre Transform in the meridional direction exploiting symmetries.Either recomputes the Legendre
polynomials `Λ` for each latitude on one hemisphere or uses precomputed polynomials from `S.Λs`, depending
on `S.recompute_legendre`. Further uses Legendre weights on Gaussian latitudes for a leakage-free
transform."""
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
        legendre_weight = legendre_weights[ilat]                    # weights normalised with π/nlat

        for m in 1:mmax+1                                           # Σ_{m=0}^{mmax}, but 1-based index
            w = legendre_weight * cis((m-1) * -lon_offset)          # apply Legendre weights on Gaussian
            an = fn[m] * w                                          # latitudes for leakage-free transform
            as = fs[m] * w
            for l in m:lmax+1
                c = isodd(l+m) ? an - as : an + as                  # odd/even wavenumbers
                alms[l,m] += c * Λ[l,m]                             # Legendre polynomials in Λ
            end
        end
    end

    return alms
end

"""
    alms = spectral(map)

Forward spectral transform (grid to spectral space) from the gridded field `map` on a regular Gaussian
grid (with Gaussian latitudes) into the spectral coefficients of the Legendre polynomials `alms`. Based
on the size of `map` this function retrieves the corresponding spectral resolution via triangular
truncation and sets up a SpectralTransform struct `S` to execute `spectral(map,S)`."""
function spectral(  map::AbstractMatrix{NF} # gridded field nlon x nlat
                    ) where NF              # number format NF

    # check grid is compatible with triangular spectral truncation
    nlon, nlat = size(map)
    @boundscheck nlon == 2nlat || throw(BoundsError)

    # Spectral resolution via triangular truncation
    trunc = triangular_truncation(nlon,nlat)    # largest trunc that satisfies the constraints
    recompute_legendre = true                   # saves memory

    # allocate spectral transform struct
    S = SpectralTransform(NF,nlon,nlat,trunc,recompute_legendre)
    return spectral(map,S)
end

"""
    alms = spectral(map,S)

Forward spectral transform (grid to spectral space) from the gridded field `map` on a regular Gaussian
grid (with Gaussian latitudes) and the SpectralTransform struct `S` into the spectral coefficients of
the Legendre polynomials `alms`. This function allocates `alms` and executes `spectral!(alms,map,S)`."""
function spectral(  map::AbstractMatrix{NF},    # gridded field nlon x nlat
                    S::SpectralTransform{NF}    # spectral transform struct
                    ) where NF                  # number format NF

    # check grid is compatible with triangular spectral truncation
    nlon, nlat = size(map)
    @boundscheck nlon == 2nlat || throw(BoundsError)
    @boundscheck nlon == S.nlon || throw(BoundsError)

    alms = Matrix{Complex{NF}}(undef,S.lmax+1,S.mmax+1)
    return spectral!(alms,map,S)                # in-place version
end