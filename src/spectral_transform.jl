"""
    S = SpectralTransform{NF<:AbstractFloat}(...)

SpectralTransform struct that contains all parameters and preallocated arrays
for the spectral transform."""
struct SpectralTransform{NF<:AbstractFloat}

    # GRID
    grid::Type{<:AbstractGrid}  # grid type used
    nresolution::Int            # resolution parameter of grid
    order::Int                  # 1,2,3 for linear, quadratic, cubic grid

    # SPECTRAL RESOLUTION
    lmax::Int               # Maximum degree l=[0,lmax] of spherical harmonics
    mmax::Int               # Maximum order m=[0,l] of spherical harmonics
    nfreq_max::Int          # Maximum (at Equator) number of Fourier frequencies (real FFT)

    # CORRESPONDING GRID SIZE
    nlon_max::Int           # Maximum number of longitude points (at Equator)
    nlons::Vector{Int}      # Number of longitude points per ring
    nlat::Int               # Number of latitude rings
    nlat_half::Int          # nlat on one hemisphere (incl equator if nlat odd)
    
    # CORRESPONDING GRID VECTORS
    colat::Vector{NF}                   # Gaussian colatitudes (0,π) North to South Pole 
    cos_colat::Vector{NF}               # Cosine of colatitudes
    sin_colat::Vector{NF}               # Sine of colatitudes
    lon_offsets::Matrix{Complex{NF}}    # Offset of first lon per ring from prime meridian

    # NORMALIZATION
    norm_sphere::NF         # normalization of the l=0,m=0 mode

    # FFT plans, one plan for each latitude ring
    rfft_plans::Vector{FFTW.rFFTWPlan{NF}}              # FFT plan for grid to spectral transform
    brfft_plans::Vector{FFTW.rFFTWPlan{Complex{NF}}}    # spectral to grid transform (inverse)

    # LEGENDRE POLYNOMIALS
    recompute_legendre::Bool                # Pre or recompute Legendre polynomials
    Λ::LowerTriangularMatrix{Float64}       # Legendre polynomials for one latitude (requires recomputing)
    Λs::Vector{LowerTriangularMatrix{NF}}   # Legendre polynomials for all latitudes (all precomputed)
    
    # QUADRATURE (integration for the Legendre polynomials, extra normalisation of π/nlat included)
    quadrature_weights::Vector{NF}          # including Δlon in ΔΩ = sinθ Δθ Δϕ (solid angle of grid point)

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
    S = SpectralTransform(NF,grid,trunc,recompute_legendre)

Generator function for a SpectralTransform struct. With `NF` the number format,
`grid` the grid type `<:AbstractGrid` and spectral truncation `trunc` this function sets up
necessary constants for the spetral transform. Also plans the Fourier transforms, retrieves the colatitudes,
and preallocates the Legendre polynomials (if recompute_legendre == false) and quadrature weights."""
function SpectralTransform( ::Type{NF},                     # Number format NF
                            grid::Type{<:AbstractGrid},     # type of spatial grid used
                            trunc::Int,                     # Spectral truncation
                            recompute_legendre::Bool) where NF

    # SPECTRAL RESOLUTION
    lmax = trunc                # Maximum degree l=[0,lmax] of spherical harmonics
    mmax = trunc                # Maximum order m=[0,l] of spherical harmonics
    order = truncation_order(grid)

    # RESOLUTION PARAMETERS
    nresolution = get_resolution(grid,trunc)        # resolution parameter, nlat_half/nside for HEALPixGrid
    nlat_half = get_nlat_half(grid,nresolution)     # contains equator for HEALPix
    nlat = 2nlat_half - nlat_odd(grid)              # one less if grids have odd # of latitude rings
    nlon_max = get_nlon(grid,nresolution)           # number of longitudes around the equator
                                                    # number of longitudes per latitude ring (one hemisphere only)
    nlons = [get_nlon_per_ring(grid,nresolution,i) for i in 1:nlat_half]
    nfreq_max = nlon_max÷2 + 1                      # maximum number of fourier frequencies (real FFTs)

    # LATITUDE VECTORS (based on Gaussian, equi-angle or HEALPix latitudes)
    colat = get_colat(grid,nresolution)             # colatitude in radians                             
    cos_colat = cos.(colat)                         # cos(colat)
    sin_colat = sin.(colat)                         # sin(colat)

    # PLAN THE FFTs
    rfft_plans = [FFTW.plan_rfft(zeros(NF,nlon)) for nlon in nlons]
    brfft_plans = [FFTW.plan_brfft(zeros(Complex{NF},nlon÷2 + 1),nlon) for nlon in nlons]

    # NORMALIZATION
    norm_sphere = 2sqrt(π)  # norm_sphere at l=0,m=0 translates to 1s everywhere in grid space
    Δlons = 2π./nlons       # normalization for forward transform to be baked into the quadrature weights

    # LONGITUDE OFFSETS OF FIRST GRID POINT PER RING (0 for full and octahedral grids)
    _, lons = get_colatlons(grid,nresolution)
    lon0s = [lons[each_index_in_ring(grid,i,nresolution)[1]] for i in 1:nlat]
    lon_offsets = [cispi(m*lon0/π) for m in 0:mmax, lon0 in lon0s]

    # PREALLOCATE LEGENDRE POLYNOMIALS, lmax+2 for one more degree l for meridional gradient recursion
    Λ = zeros(LowerTriangularMatrix,lmax+2,mmax+1)  # Legendre polynomials for one latitude

    # allocate memory in Λs for polynomials at all latitudes or allocate dummy array if precomputed
    # Λs is of size (lmax+2) x (mmax+1) x nlat_half unless recomputed, one more degree l as before
    # for recomputed only Λ is used, not Λs, create dummy array of 0-size instead
    b = ~recompute_legendre                 # true for precomputed
    Λs = [zeros(LowerTriangularMatrix,b*(lmax+2),b*(mmax+1)) for _ in 1:b*nlat_half]

    if recompute_legendre == false          # then precompute all polynomials
        for ilat in 1:nlat_half             # only one hemisphere due to symmetry
            Legendre.λlm!(Λs[ilat], lmax+1, mmax, cos_colat[ilat])
            # underflow_small!(Λs[ilat],sqrt(floatmin(NF)))
        end
    end

    # QUADRATURE WEIGHTS (Gaussian, Clenshaw-Curtis, or Riemann depending on grid)
    quadrature_weights = get_quadrature_weights(grid,nresolution)
    quadrature_weights .*= Δlons

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
    SpectralTransform{NF}(  grid,nresolution,order,
                            lmax,mmax,nfreq_max,
                            nlon_max,nlons,nlat,nlat_half,
                            colat,cos_colat,sin_colat,lon_offsets,
                            norm_sphere,
                            rfft_plans,brfft_plans,
                            recompute_legendre,Λ,Λs,quadrature_weights,
                            ϵlms,grad_x,grad_y1,grad_y2,
                            grad_y_vordiv1,grad_y_vordiv2,vordiv_to_uv_x,
                            vordiv_to_uv1,vordiv_to_uv2,
                            eigenvalues,eigenvalues⁻¹)
end

"""
    S = SpectralTransform(P::Parameters)

Generator function for a SpectralTransform struct pulling in parameters from a Parameters struct."""
function SpectralTransform(P::Parameters)
    @unpack NF, grid, trunc, recompute_legendre = P
    return SpectralTransform(NF,grid,trunc,recompute_legendre)
end

"""
    S = SpectralTransform(  alms::AbstractMatrix{Complex{NF}};
                            recompute_legendre::Bool=true,
                            grid::Type{<:AbstractGrid}=FullGaussianGrid)

Generator function for a `SpectralTransform` struct based on the size of the spectral
coefficients `alms` and the grid `grid`. Recomputes the Legendre polynomials by default."""
function SpectralTransform( alms::AbstractMatrix{Complex{NF}};  # spectral coefficients
                            recompute_legendre::Bool=true,      # saves memory
                            grid::Type{<:AbstractGrid}=FullGaussianGrid,
                            ) where NF                          # number format NF

    _, mmax = size(alms) .- 1                           # -1 for 0-based degree l, order m
    return SpectralTransform(NF,grid,mmax,recompute_legendre)
end

"""
    S = SpectralTransform(  map::AbstractGrid;
                            recompute_legendre::Bool=true)

Generator function for a `SpectralTransform` struct based on the size and grid type of
gridded field `map`. Recomputes the Legendre polynomials by default."""
function SpectralTransform( map::AbstractGrid{NF};          # gridded field
                            recompute_legendre::Bool=true,  # saves memory
                            ) where NF                      # number format NF

    grid = typeof(map)
    trunc = get_truncation(map)
    return SpectralTransform(NF,grid,trunc,recompute_legendre)
end

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
    gridded!(   map::AbstractGrid,
                alms::LowerTriangularMatrix,
                S::SpectralTransform)

Spectral transform (spectral to grid) of the spherical harmonic coefficients `alms` to a gridded field
`map`. The spectral transform is number format-flexible as long as the parametric types of `map`, `alms`, `S`
are identical. The spectral transform is grid-flexible as long as the `typeof(map)<:AbstractGrid`. 
Uses the precalculated arrays, FFT plans and other constants in the SpectralTransform struct `S`."""
function gridded!(  map::AbstractGrid{NF},                      # gridded output
                    alms::LowerTriangularMatrix{Complex{NF}},   # spectral coefficients input
                    S::SpectralTransform{NF}                    # precomputed parameters struct
                    ) where {NF<:AbstractFloat}                 # number format NF

    @unpack nlat, nlons, nlat_half, nfreq_max, order = S
    @unpack cos_colat, lon_offsets = S
    @unpack recompute_legendre, Λ, Λs = S
    @unpack brfft_plans = S

    recompute_legendre && @boundscheck size(alms) == size(Λ) || throw(BoundsError)
    recompute_legendre || @boundscheck size(alms) == size(Λs[1]) || throw(BoundsError)
    lmax, mmax = size(alms) .- 1            # maximum degree l, order m of spherical harmonics

    @boundscheck mmax+1 <= nfreq_max || throw(BoundsError)
    @boundscheck nlat == length(cos_colat) || throw(BoundsError)
    @boundscheck typeof(map) <: S.grid || throw(BoundsError)
    @boundscheck get_nresolution(map) == S.nresolution || throw(BoundsError)

    # preallocate work arrays
    gn = zeros(Complex{NF}, nfreq_max)      # phase factors for northern latitudes
    gs = zeros(Complex{NF}, nfreq_max)      # phase factors for southern latitudes

    Λw = Legendre.Work(Legendre.λlm!, Λ, Legendre.Scalar(zero(NF)))

    @inbounds for ilat_n in 1:nlat_half     # symmetry: loop over northern latitudes only
        ilat_s = nlat - ilat_n + 1          # southern latitude index
        nlon = nlons[ilat_n]                # number of longitudes on this ring
        nfreq  = nlon÷2 + 1                 # linear max Fourier frequency wrt to nlon
        nfreqm = nlon÷(order+1) + 1         # (lin/quad/cub) max frequency to shorten loop over m
        not_equator = ilat_n != ilat_s      # is the latitude ring not on equator?

        # Recalculate or use precomputed Legendre polynomials Λ
        Λ_ilat = recompute_legendre ? 
            Legendre.unsafe_legendre!(Λw, Λ, lmax, mmax, cos_colat[ilat_n]) : Λs[ilat_n]

        # inverse Legendre transform by looping over wavenumbers l,m
        lm = 1                              # single index for non-zero l,m indices
        for m in 1:min(nfreq,mmax+1)        # Σ_{m=0}^{mmax}, but 1-based index
            acc_odd  = zero(Complex{NF})    # accumulator for isodd(l+m)
            acc_even = zero(Complex{NF})    # accumulator for iseven(l+m)

            # integration over l = m:lmax+1
            lm_end = lm + lmax-m+1              # first index lm plus lmax-m+1 (length of column -1)
            even_degrees = iseven(lm+lm_end)    # is there an even number of degrees in column m?

            # anti-symmetry: sign change of odd harmonics on southern hemisphere
            # but put both into one loop for contiguous memory access
            for lm_even in lm:2:lm_end-even_degrees     
                # split into even, i.e. iseven(l+m)
                # acc_even += alms[lm_even] * Λ_ilat[lm_even], but written with muladd
                acc_even = muladd(alms[lm_even],Λ_ilat[lm_even],acc_even)

                # and odd (isodd(l+m)) harmonics
                # acc_odd += alms[lm_odd] * Λ_ilat[lm_odd], but written with muladd
                acc_odd = muladd(alms[lm_even+1],Λ_ilat[lm_even+1],acc_odd)
            end

            # for even number of degrees, one acc_even iteration is skipped, do now
            acc_even = even_degrees ? muladd(alms[lm_end],Λ_ilat[lm_end],acc_even) : acc_even

            # CORRECT FOR LONGITUDE OFFSETTS
            w = lon_offsets[m,ilat_n]           # longitude offset rotation
            acc_n = (acc_even + acc_odd)*w      # accumulators for northern
            acc_s = (acc_even - acc_odd)*w      # and southern hemisphere
            
            gn[m] += acc_n                      # accumulate in phase factors for northern
            gs[m] += acc_s                      # and southern hemisphere

            lm = lm_end + 1                     # first index of next m column
        end

        # Inverse Fourier transform in zonal direction
        brfft_plan = brfft_plans[ilat_n]        # FFT planned wrt nlon on ring
        js = each_index_in_ring(map,ilat_n)     # in-ring indices northern ring
        LinearAlgebra.mul!(view(map.v,js),brfft_plan,view(gn,1:nfreq))  # perform FFT

        # southern latitude, don't call redundant 2nd fft if ring is on equator 
        js = each_index_in_ring(map,ilat_s)     # in-ring indices southern ring
        not_equator && LinearAlgebra.mul!(view(map.v,js),brfft_plan,view(gs,1:nfreq))  # perform FFT

        fill!(gn, zero(Complex{NF}))            # set phase factors back to zero
        fill!(gs, zero(Complex{NF}))
    end

    return map
end

"""
    spectral!(  alms::LowerTriangularMatrix,
                map::AbstractGrid,
                S::SpectralTransform)

Spectral transform (grid to spectral space) from the gridded field `map` on a `grid<:AbstractGrid` to
a `LowerTriangularMatrix` of spherical harmonic coefficients `alms`. Uses FFT in the zonal direction,
and a Legendre Transform in the meridional direction exploiting symmetries. The spectral transform is
number format-flexible as long as the parametric types of `map`, `alms`, `S` are identical.
The spectral transform is grid-flexible as long as the `typeof(map)<:AbstractGrid`. 
Uses the precalculated arrays, FFT plans and other constants in the SpectralTransform struct `S`."""
function spectral!( alms::LowerTriangularMatrix{Complex{NF}},   # output: spectral coefficients
                    map::AbstractGrid{NF},                      # input: gridded values
                    S::SpectralTransform{NF}
                    ) where {NF<:AbstractFloat}
    
    @unpack nlat, nlat_half, nlons, nfreq_max, order, cos_colat = S
    @unpack recompute_legendre, Λ, Λs, quadrature_weights = S
    @unpack rfft_plans, lon_offsets = S
    
    recompute_legendre && @boundscheck size(alms) == size(Λ) || throw(BoundsError)
    recompute_legendre || @boundscheck size(alms) == size(Λs[1]) || throw(BoundsError)
    lmax, mmax = size(alms) .- 1    # maximum degree l, order m of spherical harmonics

    @boundscheck mmax+1 <= nfreq_max || throw(BoundsError)
    @boundscheck nlat == length(cos_colat) || throw(BoundsError)
    @boundscheck typeof(map) <: S.grid || throw(BoundsError)
    @boundscheck get_nresolution(map) == S.nresolution || throw(BoundsError)

    # preallocate work warrays
    fn = zeros(Complex{NF},nfreq_max)       # Fourier-transformed northern latitude
    fs = zeros(Complex{NF},nfreq_max)       # Fourier-transformed southern latitude

    # partial sums are accumulated in alms, force zeros initially.
    fill!(alms,0)

    Λw = Legendre.Work(Legendre.λlm!, Λ, Legendre.Scalar(zero(NF)))

    @inbounds for ilat_n in 1:nlat_half     # symmetry: loop over northern latitudes only
        ilat_s = nlat - ilat_n + 1          # corresponding southern latitude index
        nlon = nlons[ilat_n]                # number of longitudes on this ring
        nfreq  = nlon÷2 + 1                 # linear max Fourier frequency wrt to nlon
        nfreqm = nlon÷(order+1) + 1         # (lin/quad/cub) max frequency to shorten loop over m
        not_equator = ilat_n != ilat_s      # is the latitude ring not on equator?

        # Fourier transform in zonal direction
        rfft_plan = rfft_plans[ilat_n]          # FFT planned wrt nlon on ring
        js = each_index_in_ring(map,ilat_n)     # in-ring indices northern ring
        LinearAlgebra.mul!(view(fn,1:nfreq),rfft_plan,view(map.v,js))   # Northern latitude

        js = each_index_in_ring(map,ilat_s)     # in-ring indices northern ring
                                                # Southern latitude (don't call FFT on Equator)
                                                # then fill fs with zeros and no changes needed further down
        not_equator ? LinearAlgebra.mul!(view(fs,1:nfreq),rfft_plan,view(map.v,js)) : fill!(fs,0)

        # Legendre transform in meridional direction
        # Recalculate or use precomputed Legendre polynomials Λ
        Λ_ilat = recompute_legendre ?
            Legendre.unsafe_legendre!(Λw, Λ, lmax, mmax, cos_colat[ilat_n]) : Λs[ilat_n]
        quadrature_weight = quadrature_weights[ilat_n]  # weights normalised with π/nlat

        lm = 1                                          # single index for spherical harmonics
        for m in 1:min(nfreq,mmax+1)                    # Σ_{m=0}^{mmax}, but 1-based index

            # QUADRATURE WEIGHTS and LONGITUDE OFFSET
            w = lon_offsets[m,ilat_n]                   # longitude offset rotation
            quadrature_weight *= conj(w)                # complex conjugate for back transform
            an = fn[m] * quadrature_weight              # weighted northern latitude
            as = fs[m] * quadrature_weight              # weighted southern latitude
            
            # LEGENDRE TRANSFORM
            a_even = an + as                            # sign flip due to anti-symmetry with
            a_odd = an - as                             # odd polynomials 

            # integration over l = m:lmax+1
            lm_end = lm + lmax-m+1                      # first index lm plus lmax-m+1 (length of column -1)
            even_degrees = iseven(lm+lm_end)            # is there an even number of degrees in column m?
            
            # anti-symmetry: sign change of odd harmonics on southern hemisphere
            # but put both into one loop for contiguous memory access
            for lm_even in lm:2:lm_end-even_degrees
                # split into even, i.e. iseven(l+m)
                # alms[lm_even] += a_even * Λ_ilat[lm_even], but written with muladd
                alms[lm_even] = muladd(a_even,Λ_ilat[lm_even],alms[lm_even])
                
                # and odd (isodd(l+m)) haxwxwrmonics
                # alms[lm_odd] += a_odd * Λ_ilat[lm_odd], but written with muladd
                alms[lm_even+1] = muladd(a_odd,Λ_ilat[lm_even+1],alms[lm_even+1])
            end

            # for even number of degrees, one even iteration is skipped, do now
            alms[lm_end] = even_degrees ? muladd(a_even,Λ_ilat[lm_end],alms[lm_end]) : alms[lm_end]

            lm = lm_end + 1                             # first index of next m column
        end
    end

    return alms
end

"""
    map = gridded(  alms::AbstractMatrix;
                    recompute_legendre::Bool=true,
                    grid::Type{<:AbstractGrid}=FullGaussianGrid)

Spectral transform (spectral to grid space) from spherical coefficients `alms` to a newly allocated gridded
field `map`. Based on the size of `alms` the grid type `grid`, the spatial resolution is retrieved based
on the truncation defined for `grid`. SpectralTransform struct `S` is allocated to execute `gridded(alms,S)`."""
function gridded(   alms::AbstractMatrix{Complex{NF}};  # spectral coefficients
                    recompute_legendre::Bool=true,      # saves memory
                    grid::Type{<:AbstractGrid}=FullGaussianGrid,
                    ) where NF                          # number format NF

    _, mmax = size(alms) .- 1                           # -1 for 0-based degree l, order m
    S = SpectralTransform(NF,grid,mmax,recompute_legendre)
    return gridded(alms,S)
end

"""
    map = gridded(  alms::AbstractMatrix,
                    S::SpectralTransform)

Spectral transform (spectral to grid space) from spherical coefficients `alms` to a newly allocated gridded
field `map` with precalculated properties based on the SpectralTransform struct `S`. `alms` is converted to
a `LowerTriangularMatrix` to execute the in-place `gridded!`."""
function gridded(   alms::AbstractMatrix{Complex{NF}},  # spectral coefficients
                    S::SpectralTransform{NF}            # struct for spectral transform parameters
                    ) where NF                          # number format NF
 
    map = zeros(S.grid{NF},S.nresolution)       # preallocate output
    almsᴸ = LowerTriangularMatrix(alms)         # drop the upper triangle entries
    gridded!(map,almsᴸ,S)                       # now execute the in-place version
    return map
end

"""
    alms = spectral(    map::AbstractMatrix;
                        grid::Type{<:AbstractGrid}=FullGaussianGrid,
                        kwargs...)

Converts `map` to `grid(map)` to execute spectral(map::AbstractGrid;kwargs...)."""
function spectral(  map::AbstractMatrix;            # gridded field
                    grid::Type{<:AbstractGrid}=FullGaussianGrid,
                    kwargs...)

    return spectral(grid(map),kwargs...)
end

"""
    alms = spectral(    map::AbstractMatrix;
                        grid::Type{<:AbstractGrid}=FullGaussianGrid,
                        kwargs...)

Converts `map` to `grid(map)` to execute spectral(map::AbstractGrid;kwargs...)."""
function spectral(  map::AbstractGrid{NF};          # gridded field
                    recompute_legendre::Bool=true,  # saves memory
                    ) where NF                      # number format NF

    grid = typeof(map)
    trunc = get_truncation(map)
    S = SpectralTransform(NF,grid,trunc,recompute_legendre)
    return spectral(map,S)
end

"""
    alms = spectral(    map::AbstractMatrix;
                        grid::Type{<:AbstractGrid}=FullGaussianGrid,
                        kwargs...)

Spectral transform (grid to spectral) `map` to `grid(map)` to execute spectral(map::AbstractGrid;kwargs...)."""
function spectral(  map::AbstractGrid{NF},      # gridded field
                    S::SpectralTransform{NF},   # spectral transform struct
                    ) where NF                  # number format NF

    # always use one more l for consistency with vector quantities
    alms = LowerTriangularMatrix{Complex{NF}}(undef,S.lmax+2,S.mmax+1)
    return spectral!(alms,map,S)                # in-place version
end