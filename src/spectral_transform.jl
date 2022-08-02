"""SpectralTransform struct that contains all parameters and preallocated arrays
for the spectral transform."""
struct SpectralTransform{NF<:AbstractFloat}
    # SPECTRAL RESOLUTION
    lmax::Int       # Maximum degree l=[0,lmax] of spherical harmonics
    mmax::Int       # Maximum order m=[0,l] of spherical harmonics
    nfreq::Int      # Number of fourier frequencies (real FFT)
    radius::Real    # radius of the sphere/Earth

    # CORRESPONDING GRID SIZE
    nlon::Int               # Number of longitudes
    nlat::Int               # Number of latitudes
    nlat_half::Int          # nlat on one hemisphere
    
    # CORRESPONDING GRID VECTORS
    colat::Vector{NF}       # Gaussian colatitudes (0,π) North to South Pole 
    cos_colat::Vector{NF}   # Cosine of colatitudes
    sin_colat::Vector{NF}   # Sine of colatitudes
    lon_offset::NF          # Offset of first longitude from prime meridian

    # NORMALIZATION
    norm_sphere::NF         # normalization of the l=0,m=0 mode
    norm_forward::NF        # normalization of the Legendre weights for forward transform

    # FFT plans
    rfft_plan::FFTW.rFFTWPlan{NF}           # grid to spectral transform
    brfft_plan::FFTW.rFFTWPlan{Complex{NF}} # spectral to grid transform (inverse)

    # FFT work arrays (currently not used)
    gn::Vector{Complex{NF}}         # phase factors north
    gs::Vector{Complex{NF}}         # phase factors south
    fn::Vector{Complex{NF}}         # Fourier-transformed latitude north
    fs::Vector{Complex{NF}}         # Fourier-transformed latitude south

    # LEGENDRE POLYNOMIALS
    recompute_legendre::Bool                # Pre or recompute Legendre polynomials
    Λ::LowerTriangularMatrix{NF}            # Legendre polynomials for one latitude (requires recomputing)
    Λs::Vector{LowerTriangularMatrix{NF}}   # Legendre polynomials for all latitudes (all precomputed)
    legendre_weights::Vector{NF}            # Legendre weights (extra normalisation of π/nlat included)

    # RECURSION FACTORS
    ϵlms::LowerTriangularMatrix{NF}         # precomputed for meridional gradients gradients grad_y1, grad_y2

    # GRADIENT MATRICES (on unit sphere, no 1/radius-scaling included)
    grad_x ::Vector{Complex{NF}}            # = i*m but precomputed
    grad_y1::LowerTriangularMatrix{NF}      # precomputed meridional gradient factors, term 1
    grad_y2::LowerTriangularMatrix{NF}      # term 2

    # GRADIENT MATRICES FOR U,V -> Vorticity,Divergence
    grad_y_vordiv1::LowerTriangularMatrix{NF}
    grad_y_vordiv2::LowerTriangularMatrix{NF}

    # GRADIENT MATRICES FOR Vorticity,Divergence -> U,V
    vordiv_to_uv_x::LowerTriangularMatrix{NF}
    vordiv_to_uv1::LowerTriangularMatrix{NF}
    vordiv_to_uv2::LowerTriangularMatrix{NF}

    # EIGENVALUES (on unit sphere, no 1/radius²-scaling included)
    eigenvalues  ::Vector{NF}               # = -l*(l+1), degree l of spherical harmonic
    eigenvalues⁻¹::Vector{NF}               # = -1/(l*(l+1))
end

"""
    S = SpectralTransform(NF,nlon,nlat,trunc,recompute_legendre)

Generator function for a SpectralTransform struct. From number of longitudes `nlon`,
number of latitudes `nlat` and spectral truncation `trunc` this function sets up the
spectral resolution, plans the Fourier transforms, retrieves the Gaussian colatitudes,
and preallocates the Legendre polynomials (if recompute_legendre == false) and legendre
weights.
"""
function SpectralTransform( ::Type{NF},     # Number format NF
                            nlon::Int,      # Number of longitudes
                            nlat::Int,      # Number of latitudes
                            trunc::Int,     # Spectral truncation
                            radius::Real,   # radius of sphere/Earth
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

    # FFT work arrays (currently not used)
    gn = zeros(Complex{NF},nfreq)       # phase factors north
    gs = zeros(Complex{NF},nfreq)       # phase factors south
    fn = zeros(Complex{NF},nfreq)       # Fourier-transformed latitude north
    fs = zeros(Complex{NF},nfreq)       # Fourier-transformed latitude south

    # GAUSSIAN COLATITUDES (0,π) North to South
    nodes = FastGaussQuadrature.gausslegendre(nlat)[1]  # zeros of the Legendre polynomial
    colat = π .- acos.(nodes)                           # corresponding colatitudes
    sin_colat = sin.(colat)                             
    cos_colat = cos.(colat)
    lon_offset = π/nlon                                 # offset of first longitude from prime meridian

    # NORMALIZATION
    norm_sphere = 2sqrt(π)      # norm_sphere at l=0,m=0 translates to 1s everywhere in grid space
    norm_forward = π/nlat       # normalization for forward transform to be baked into the Legendre weights

    # PREALLOCATE LEGENDRE POLYNOMIALS, lmax+2 for one more degree l for meridional gradient recursion
    Λ = zeros(LowerTriangularMatrix,lmax+2,mmax+1)  # Legendre polynomials for one latitude

    # allocate memory in Λs for polynomials at all latitudes or allocate dummy array if precomputed
    # Λs is of size (lmax+2) x (mmax+1) x nlat_half unless recomputed, one more degree l as before
    # for recomputed only Λ is used, not Λs, create dummy array of size 1x1x1 instead
    b = ~recompute_legendre                 # true for precomputed
    # Λs = zeros(b*lmax + 1 + b, b*mmax + 1, b*(nlat_half-1) + 1) 
    Λs = [zeros(LowerTriangularMatrix,b*lmax+1+b,b*mmax+1) for _ in 1:b*(nlat_half-1)+1]

    if recompute_legendre == false          # then precompute all polynomials
        for ilat in 1:nlat_half             # only one hemisphere due to symmetry
            AssociatedLegendrePolynomials.λlm!(Λs[ilat], lmax+1, mmax, cos_colat[ilat])
        end
    end

    # Legendre weights
    legendre_weights = FastGaussQuadrature.gausslegendre(nlat)[2][1:nlat_half]
    legendre_weights *= norm_forward        # extra normalisation for forward transform included

    # RECURSION FACTORS
    ϵlms = get_recursion_factors(lmax+1,mmax)

    # GRADIENTS (on unit sphere, hence 1/radius-scaling is omitted)
    grad_x = [im*m for m in 0:mmax+1]       # zonal gradient (precomputed currently not used)

    # meridional gradient for scalars (coslat scaling included)
    grad_y1 = zeros(LowerTriangularMatrix,lmax+2,mmax+1)          # term 1, mul with harmonic l-1,m
    grad_y2 = zeros(LowerTriangularMatrix,lmax+2,mmax+1)          # term 2, mul with harmonic l+1,m

    for m in 0:mmax                         # 0-based degree l, order m
        for l in m:lmax+1           
            grad_y1[l+1,m+1] = -(l-1)*ϵlms[l+1,m+1]
            grad_y2[l+1,m+1] = (l+2)*ϵlms[l+2,m+1]
        end
    end

    # meridional gradient used to get from u,v/coslat to vorticity and divergence
    grad_y_vordiv1 = zeros(LowerTriangularMatrix,lmax+2,mmax+1)   # term 1, mul with harmonic l-1,m
    grad_y_vordiv2 = zeros(LowerTriangularMatrix,lmax+2,mmax+1)   # term 2, mul with harmonic l+1,m

    for m in 0:mmax                         # 0-based degree l, order m
        for l in m:lmax+1           
            grad_y_vordiv1[l+1,m+1] = (l+1)*ϵlms[l+1,m+1]
            grad_y_vordiv2[l+1,m+1] = l*ϵlms[l+2,m+1]
        end
    end

    # zonal integration (sort of) to get from vorticity and divergence to u,v*coslat
    vordiv_to_uv_x = LowerTriangularMatrix([-m/(l*(l+1)) for l in 0:lmax+1, m in 0:mmax])
    vordiv_to_uv_x[1,1] = 0

    # meridional integration (sort of) to get from vorticity and divergence to u,v*coslat
    vordiv_to_uv1 = zeros(LowerTriangularMatrix,lmax+2,mmax+1)    # term 1, to be mul with harmonic l-1,m
    vordiv_to_uv2 = zeros(LowerTriangularMatrix,lmax+2,mmax+1)    # term 2, to be mul with harmonic l+1,m

    for m in 0:mmax                         # 0-based degree l, order m
        for l in m:lmax+1           
            vordiv_to_uv1[l+1,m+1] = ϵlms[l+1,m+1]/l
            vordiv_to_uv2[l+1,m+1] = ϵlms[l+2,m+1]/(l+1)
        end
    end

    vordiv_to_uv1[1,1] = 0                  # remove NaN from 0/0

    # EIGENVALUES (on unit sphere, hence 1/radius²-scaling is omitted)
    eigenvalues = [-l*(l+1) for l in 0:lmax+1]
    eigenvalues⁻¹ = inv.(eigenvalues)
    eigenvalues⁻¹[1] = 0                    # set the integration constant to 0
        
    # conversion to NF happens here
    SpectralTransform{NF}(  lmax,mmax,nfreq,radius,
                            nlon,nlat,nlat_half,
                            colat,cos_colat,sin_colat,lon_offset,
                            norm_sphere,norm_forward,
                            rfft_plan,brfft_plan,
                            gn,gs,fn,fs,
                            recompute_legendre,Λ,Λs,legendre_weights,
                            ϵlms,grad_x,grad_y1,grad_y2,
                            grad_y_vordiv1,grad_y_vordiv2,vordiv_to_uv_x,
                            vordiv_to_uv1,vordiv_to_uv2,
                            eigenvalues,eigenvalues⁻¹)
end

"""Generator function for a SpectralTransform struct in case the number format is not provided.
Use Float64 as default."""
SpectralTransform(args...) = SpectralTransform(Float64,args...)

"Generator function for a SpectralTransform struct pulling in parameters from a Parameters struct."
function SpectralTransform(P::Parameters)
    T = triangular_truncation(trunc=P.trunc)
    return SpectralTransform(P.NF,T.nlon,T.nlat,P.trunc,P.radius_earth,P.recompute_legendre)
end

# don't recompute triangular truncation if Geometry is already available
SpectralTransform(P::Parameters,G::Geometry) = SpectralTransform(P.NF,G.nlon,G.nlat,P.trunc,
                                                                P.radius_earth,P.recompute_legendre)

"""
    ϵ = ϵ(NF,l,m) 

Recursion factors `ϵ` as a function of degree `l` and order `m` (0-based) of the spherical harmonics.
ϵ(l,m) = sqrt((l^2-m^2)/(4*l^2-1)) and then converted to number format NF."""
function ϵlm(::Type{NF},l::Int,m::Int) where NF
    return convert(NF,sqrt((l^2-m^2)/(4*l^2-1)))
end

"""
    ϵ = ϵ(l,m) 

Recursion factors `ϵ` as a function of degree `l` and order `m` (0-based) of the spherical harmonics.
ϵ(l,m) = sqrt((l^2-m^2)/(4*l^2-1)) with default number format Float64."""
ϵlm(l::Int,m::Int) = ϵlm(Float64,l,m)

"""
    get_recursion_factors(  ::Type{NF}, # number format NF
                            lmax::Int,  # max degree l of spherical harmonics (0-based here)
                            mmax::Int   # max order m of spherical harmonics
                            ) where {NF<:AbstractFloat}
        
Returns a matrix of recursion factors `ϵ` up to degree `lmax` and order `mmax` of number format `NF`."""
function get_recursion_factors( ::Type{NF}, # number format NF
                                lmax::Int,  # max degree l of spherical harmonics (0-based here)
                                mmax::Int   # max order m of spherical harmonics
                                ) where {NF<:AbstractFloat}

    # preallocate array with one more l for meridional gradients
    ϵlms = zeros(LowerTriangularMatrix{NF},lmax+2,mmax+1)      
    for m in 1:mmax+1                   # loop over 1-based l,m
        for l in m:lmax+2
            ϵlms[l,m] = ϵlm(NF,l-1,m-1) # convert to 0-based l,m for function call
        end
    end
    return ϵlms
end

# if number format not provided use Float64
get_recursion_factors(lmax::Int,mmax::Int) = get_recursion_factors(Float64,lmax,mmax)

"""
    Λ_ilat = get_legendre_polynomials!(Λ,Λs,ilat,cos_colat,recompute_legendre)

Base on `recompute_legendre` (true/false) this function either updates the Legendre polynomials
`Λ` for a given latitude `ilat, cos_colat` by recomputation (`recompute_legendre == true`), or
`Λ` is changed by creating a view on the corresponding latitude in precomputed `Λs`.
Recomputation takes usually longer, but precomputation requires a large amount of memory for high resolution."""
function get_legendre_polynomials!( Λ::LowerTriangularMatrix{NF},           # Out: Polynomials at lat
                                    Λs::Vector{LowerTriangularMatrix{NF}},  # Precomputed polynomials
                                    ilat::Int,                              # latitude index
                                    cos_colat::NF,                          # cosine of colatitude
                                    recompute_legendre::Bool                # recompute ignores Λs, but uses cos_colat
                                    ) where NF                              # Number format NF

    if recompute_legendre
        # Recalculate the (normalized) λ_l^m(cos(colat)) factors of the ass. Legendre polynomials
        lmax,mmax = size(Λ) .- 1
        AssociatedLegendrePolynomials.λlm!(Λ, lmax, mmax, cos_colat)
        return Λ
    else    # view on precomputed values
        return Λs[ilat]
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
                    alms::LowerTriangularMatrix{Complex{NF}},   # spectral coefficients input
                    S::SpectralTransform{NF}                    # precomputed parameters struct
                    ) where {NF<:AbstractFloat}                 # number format NF

    lmax, mmax = size(alms) .- 1            # maximum degree l, order m of spherical harmonics
    nlon, nlat = size(map)                  # number of longitudes, latitudes in grid space
    nlat_half = (nlat+1) ÷ 2                # half the number of longitudes
    nfreq = nlon÷2 + 1                      # Number of fourier frequencies (real FFTs)

    @unpack cos_colat, lon_offset = S
    @unpack recompute_legendre, Λ, Λs = S
    @unpack brfft_plan = S

    @boundscheck (lmax, mmax) <= size(Λ) .- 1 || throw(BoundsError) # -1 for 0-based lmax,mmax
    @boundscheck mmax+1 <= nfreq || throw(BoundsError)
    @boundscheck nlat == length(cos_colat) || throw(BoundsError)

    # preallocate work arrays
    gn = zeros(Complex{NF}, nfreq)          # phase factors for northern latitudes
    gs = zeros(Complex{NF}, nfreq)          # phase factors for southern latitudes

    @inbounds for ilat in 1:nlat_half       # loop over northern latitudes only due to symmetry
        ilat_s = nlat - ilat + 1            # southern latitude index

        # Recalculate or use precomputed Legendre polynomials Λ
        Λ_ilat = get_legendre_polynomials!(Λ,Λs,ilat,cos_colat[ilat],recompute_legendre)

        # inverse Legendre transform by looping over wavenumbers l,m
        lm = 0                              # single index for non-zero l,m indices in LowerTriangularMatrix
        for m in 1:mmax+1                   # Σ_{m=0}^{mmax}, but 1-based index
            accn = zero(Complex{NF})        # accumulators for northern
            accs = zero(Complex{NF})        # and southern hemisphere

            for l in m:lmax+1               # Σ_{l=m}^{lmax}, but 1-based index
                lm += 1                     # next non-zero coeffs
                term = alms[lm] * Λ_ilat[lm]# Legendre polynomials in Λ at latitude
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

        fill!(gn, zero(Complex{NF}))        # set phase factors back to zero
        fill!(gs, zero(Complex{NF}))
    end

    return map
end

"""
    map = gridded(alms)

Backward or inverse spectral transform (spectral to grid space) from coefficients `alms`. Based on the size
of `alms` the corresponding grid space resolution is retrieved based on triangular truncation and a 
SpectralTransform struct `S` is allocated to execute `gridded(alms,S)`."""
function gridded(   alms::AbstractMatrix{Complex{NF}};  # spectral coefficients
                    recompute_legendre::Bool=true       # saves memory
                    ) where NF                          # number format NF

    _, mmax = size(alms) .- 1                           # -1 for 0-based degree l, order m

    # get grid size from spectral resolution via triangular_truncation
    # use mmax instead of lmax in case lmax = mmax + 1 (required in the meridional gradient recursion)
    T = triangular_truncation(trunc=mmax)           
    @unpack nlon, nlat = T                              # number of longitudes, number of latitudes
    
    radius = 1                                          # only needed for SpectralTransform() but not used
    S = SpectralTransform(NF,nlon,nlat,mmax,radius,recompute_legendre)
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

    output = Matrix{NF}(undef,S.nlon,S.nlat)    # preallocate output
    input = LowerTriangularMatrix(alms)         # drop the upper triangle entries
    gridded!(output,input,S)                    # now execute the in-place version
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
function spectral!( alms::LowerTriangularMatrix{Complex{NF}},   # output: spectral coefficients
                    map::AbstractMatrix{NF},                    # input: gridded values
                    S::SpectralTransform{NF}
                    ) where {NF<:AbstractFloat}
    
    lmax, mmax = size(alms) .- 1            # maximum degree l, order m of spherical harmonics
    nlon, nlat = size(map)                  # number of longitudes, latitudes in grid space
    nlat_half = (nlat+1) ÷ 2                # half the number of longitudes
    nfreq = nlon÷2 + 1                      # Number of fourier frequencies (real FFTs)

    @unpack cos_colat, sin_colat, lon_offset = S
    @unpack recompute_legendre, Λ, Λs, legendre_weights = S
    @unpack rfft_plan = S
    
    @boundscheck (lmax, mmax) <= size(Λ) .- 1 || throw(BoundsError)
    @boundscheck mmax+1 <= nfreq || throw(BoundsError)

    # preallocate work warrays
    fn = zeros(Complex{NF},nfreq)   # Fourier-transformed northern latitude
    fs = zeros(Complex{NF},nfreq)   # Fourier-transformed southern latitude

    # partial sums are accumulated in alms, force zeros initially.
    fill!(alms,0)   

    @inbounds for ilat in 1:nlat_half   # loop over northern latitudes only due to symmetry
        ilat_s = nlat - ilat + 1        # corresponding southern latitude index

        # Fourier transform in zonal direction
        LinearAlgebra.mul!(fn, rfft_plan, @view(map[:,ilat]))       # Northern latitude
        LinearAlgebra.mul!(fs, rfft_plan, @view(map[:,ilat_s]))     # Southern latitude

        # Legendre transform in meridional direction
        # Recalculate or use precomputed Legendre polynomials Λ
        Λ_ilat = get_legendre_polynomials!(Λ,Λs,ilat,cos_colat[ilat],recompute_legendre)
        legendre_weight = legendre_weights[ilat]                    # weights normalised with π/nlat

        lm = 0                                                      # single index for spherical harmonics
        for m in 1:mmax+1                                           # Σ_{m=0}^{mmax}, but 1-based index
            w = legendre_weight * cis((m-1) * -lon_offset)          # apply Legendre weights on Gaussian
            an = fn[m] * w                                          # latitudes for leakage-free transform
            as = fs[m] * w
            for l in m:lmax+1
                lm += 1                                             # next coefficient of spherical harmonics
                c = isodd(l+m) ? an - as : an + as                  # odd/even wavenumbers
                alms[lm] += c * Λ_ilat[lm]                          # Legendre polynomials in Λ at latitude
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
function spectral(  map::AbstractMatrix{NF};        # gridded field nlon x nlat
                    recompute_legendre::Bool=true,  # saves memory
                    kwargs...                       # additional keyword arguments
                    ) where NF                      # number format NF

    # check grid is compatible with triangular spectral truncation
    nlon, nlat = size(map)
    @boundscheck nlon == 2nlat || throw(BoundsError)

    # Spectral resolution via triangular truncation
    T = triangular_truncation(;nlon)                # largest trunc that satisfies the constraints
    @unpack trunc = T                               # spectral truncation

    # allocate spectral transform struct
    radius = 1  # only needed for argument compatibility
    S = SpectralTransform(NF,nlon,nlat,trunc,radius,recompute_legendre)
    return spectral(map,S;kwargs...)
end

"""
    alms = spectral(map,S)

Forward spectral transform (grid to spectral space) from the gridded field `map` on a regular Gaussian
grid (with Gaussian latitudes) and the SpectralTransform struct `S` into the spectral coefficients of
the Legendre polynomials `alms`. This function allocates `alms` and executes `spectral!(alms,map,S)`."""
function spectral(  map::AbstractMatrix{NF},    # gridded field nlon x nlat
                    S::SpectralTransform{NF};   # spectral transform struct
                    one_more_l::Bool=false      # additional degree l for recursion?
                    ) where NF                  # number format NF

    # check grid is compatible with triangular spectral truncation
    nlon, nlat = size(map)
    @boundscheck nlon == 2nlat || throw(BoundsError)
    @boundscheck nlon == S.nlon || throw(BoundsError)

    alms = LowerTriangularMatrix{Complex{NF}}(undef,S.lmax+1+one_more_l,S.mmax+1)
    return spectral!(alms,map,S)                # in-place version
end