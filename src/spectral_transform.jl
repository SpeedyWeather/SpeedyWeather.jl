"""
    S = SpectralTransform{NF<:AbstractFloat}(...)

SpectralTransform struct that contains all parameters and preallocated arrays
for the spectral transform."""
struct SpectralTransform{NF<:AbstractFloat}

    # GRID
    Grid::Type{<:AbstractGrid}  # grid type used
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
    rfft_plans::Vector{AbstractFFTs.Plan}     # FFT plan for grid to spectral transform
    brfft_plans::Vector{AbstractFFTs.Plan}    # spectral to grid transform (inverse)

    # LEGENDRE POLYNOMIALS
    recompute_legendre::Bool                # Pre or recompute Legendre polynomials
    Λ::LowerTriangularMatrix{NF}            # Legendre polynomials for one latitude (requires recomputing)
    Λs::Vector{LowerTriangularMatrix{NF}}   # Legendre polynomials for all latitudes (all precomputed)
    
    # SOLID ANGLES ΔΩ FOR QUADRATURE (integration for the Legendre polynomials, extra normalisation of π/nlat included)
    solid_angles::Vector{NF}                # = ΔΩ = sinθ Δθ Δϕ (solid angle of grid point)

    # RECURSION FACTORS
    ϵlms::LowerTriangularMatrix{NF}         # precomputed for meridional gradients grad_y1, grad_y2

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
    S = SpectralTransform(NF,Grid,trunc,recompute_legendre)

Generator function for a SpectralTransform struct. With `NF` the number format,
`Grid` the grid type `<:AbstractGrid` and spectral truncation `trunc` this function sets up
necessary constants for the spetral transform. Also plans the Fourier transforms, retrieves the colatitudes,
and preallocates the Legendre polynomials (if recompute_legendre == false) and quadrature weights."""
function SpectralTransform( ::Type{NF},                     # Number format NF
                            Grid::Type{<:AbstractGrid},     # type of spatial grid used
                            trunc::Int,                     # Spectral truncation
                            recompute_legendre::Bool) where NF

    # SPECTRAL RESOLUTION
    lmax = trunc                # Maximum degree l=[0,lmax] of spherical harmonics
    mmax = trunc                # Maximum order m=[0,l] of spherical harmonics
    order = truncation_order(Grid)

    # RESOLUTION PARAMETERS
    nresolution = get_resolution(Grid,trunc)        # resolution parameter nlat_half
    nlat_half = nresolution                         # number of latitude rings on one hemisphere incl equator
    nlat = get_nlat(Grid,nlat_half)                 # 2nlat_half but one less if grids have odd # of lat rings
    nlon_max = get_nlon_max(Grid,nlat_half)         # number of longitudes around the equator
                                                    # number of longitudes per latitude ring (one hemisphere only)
    nlons = [get_nlon_per_ring(Grid,nlat_half,j) for j in 1:nlat_half]
    nfreq_max = nlon_max÷2 + 1                      # maximum number of fourier frequencies (real FFTs)

    # LATITUDE VECTORS (based on Gaussian, equi-angle or HEALPix latitudes)
    colat = get_colat(Grid,nlat_half)               # colatitude in radians                             
    cos_colat = cos.(colat)                         # cos(colat)
    sin_colat = sin.(colat)                         # sin(colat)

    # PLAN THE FFTs
    FFT_package = NF <: Union{Float32,Float64} ? FFTW : GenericFFT
    rfft_plans = [FFT_package.plan_rfft(zeros(NF,nlon)) for nlon in nlons]
    brfft_plans = [FFT_package.plan_brfft(zeros(Complex{NF},nlon÷2 + 1),nlon) for nlon in nlons]
    
    # NORMALIZATION
    norm_sphere = 2sqrt(π)  # norm_sphere at l=0,m=0 translates to 1s everywhere in grid space

    # LONGITUDE OFFSETS OF FIRST GRID POINT PER RING (0 for full and octahedral grids)
    _, lons = get_colatlons(Grid,nlat_half)
    lon1s = [lons[each_index_in_ring(Grid,j,nlat_half)[1]] for j in 1:nlat_half]
    lon_offsets = [cispi(m*lon1/π) for m in 0:mmax, lon1 in lon1s]
    Grid <: AbstractHEALPixGrid && fill!(lon_offsets,1)     # no rotation for HEALPix at the moment
    
    # PREALLOCATE LEGENDRE POLYNOMIALS, lmax+2 for one more degree l for meridional gradient recursion
    Λ = zeros(LowerTriangularMatrix{NF},lmax+2,mmax+1)  # Legendre polynomials for one latitude

    # allocate memory in Λs for polynomials at all latitudes or allocate dummy array if precomputed
    # Λs is of size (lmax+2) x (mmax+1) x nlat_half unless recomputed, one more degree l as before
    # for recomputed only Λ is used, not Λs, create dummy array of 0-size instead
    b = ~recompute_legendre                 # true for precomputed
    Λs = [zeros(LowerTriangularMatrix{NF},b*(lmax+2),b*(mmax+1)) for _ in 1:nlat_half]

    if recompute_legendre == false          # then precompute all polynomials
        Λtemp = zeros(NF,lmax+2,mmax+1)     # preallocate matrix
        for j in 1:nlat_half                # only one hemisphere due to symmetry
            Legendre.λlm!(Λtemp, lmax+1, mmax, Float64(cos_colat[j]))   # always precalculate in Float64
            # underflow!(Λtemp,sqrt(floatmin(NF)))
            copyto!(Λs[j],Λtemp)
        end
    end

    # SOLID ANGLES WITH QUADRATURE WEIGHTS (Gaussian, Clenshaw-Curtis, or Riemann depending on grid)
    # solid angles are ΔΩ = sinθ Δθ Δϕ for every grid point with
    # sin(θ)dθ are the quadrature weights approximate the integration over latitudes
    # and sum up to 2 over all latitudes as ∫sin(θ)dθ = 2 over 0...π.
    # Δϕ = 2π/nlon is the azimuth every grid point covers
    solid_angles = get_solid_angles(Grid,nlat_half)

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
    SpectralTransform{NF}(  Grid,nresolution,order,
                            lmax,mmax,nfreq_max,
                            nlon_max,nlons,nlat,nlat_half,
                            colat,cos_colat,sin_colat,lon_offsets,
                            norm_sphere,
                            rfft_plans,brfft_plans,
                            recompute_legendre,Λ,Λs,solid_angles,
                            ϵlms,grad_x,grad_y1,grad_y2,
                            grad_y_vordiv1,grad_y_vordiv2,vordiv_to_uv_x,
                            vordiv_to_uv1,vordiv_to_uv2,
                            eigenvalues,eigenvalues⁻¹)
end

"""
    S2 = spectral_transform_for_full_grid(S::SpectralTransform)

Create a spectral transform struct `S2` similar to the input `S`, but for the corresponding
full grid of the grid in `S`. The FFT is replanned and lon_offsets are set to 1 (i.e. no rotation).
Solid angles for the Legendre transform are recomputed, but all other arrays fields for S, S2
point to the same place in memory, e.g. the Legendre polynomials aren't recomputed or stored twice."""
function spectral_transform_for_full_grid(S::SpectralTransform{NF}) where NF

    FullGrid = full_grid(S.Grid)    # corresponding full grid
    
    # unpack everything that does not have to be recomputed for the full grid
    @unpack nresolution, order, lmax, mmax, nfreq_max, nlon_max, nlat, nlat_half = S
    @unpack colat, cos_colat, sin_colat, norm_sphere, recompute_legendre, Λ, Λs = S
    @unpack ϵlms, grad_x, grad_y1, grad_y2, grad_y_vordiv1, grad_y_vordiv2 = S
    @unpack vordiv_to_uv_x, vordiv_to_uv1, vordiv_to_uv2, eigenvalues, eigenvalues⁻¹ = S

    # recalculate what changes on the full grid: FFT and offsets (always 1)
    nlons = [get_nlon_per_ring(FullGrid,nlat_half,j) for j in 1:nlat_half]
    FFT_package = NF <: Union{Float32,Float64} ? FFTW : GenericFFT
    rfft_plans = [FFT_package.plan_rfft(zeros(NF,nlon)) for nlon in nlons]
    brfft_plans = [FFT_package.plan_brfft(zeros(Complex{NF},nlon÷2 + 1),nlon) for nlon in nlons]
    
    lon_offsets = copy(S.lon_offsets)
    fill!(lon_offsets,1)            # =*1, i.e. no longitude offset rotation on full grid

    solid_angles = get_solid_angles(FullGrid,nlat_half)

    return SpectralTransform{NF}(   FullGrid,nresolution,order,
                                    lmax,mmax,nfreq_max,
                                    nlon_max,nlons,nlat,nlat_half,
                                    colat,cos_colat,sin_colat,lon_offsets,
                                    norm_sphere,
                                    rfft_plans,brfft_plans,
                                    recompute_legendre,Λ,Λs,solid_angles,
                                    ϵlms,grad_x,grad_y1,grad_y2,
                                    grad_y_vordiv1,grad_y_vordiv2,vordiv_to_uv_x,
                                    vordiv_to_uv1,vordiv_to_uv2,
                                    eigenvalues,eigenvalues⁻¹)
end

"""
    S = SpectralTransform(P::Parameters)

Generator function for a SpectralTransform struct pulling in parameters from a Parameters struct."""
function SpectralTransform(P::Parameters)
    @unpack NF, Grid, trunc, recompute_legendre = P
    return SpectralTransform(NF,Grid,trunc,recompute_legendre)
end

"""
    S = SpectralTransform()

As-empty-as-possible constructor for an instance of SpectralTransform."""
function SpectralTransform()
    return SpectralTransform(Float64,FullGaussianGrid,0,true)
end

"""
    S = SpectralTransform(  alms::AbstractMatrix{Complex{NF}};
                            recompute_legendre::Bool=true,
                            Grid::Type{<:AbstractGrid}=FullGaussianGrid)

Generator function for a `SpectralTransform` struct based on the size of the spectral
coefficients `alms` and the grid `Grid`. Recomputes the Legendre polynomials by default."""
function SpectralTransform( alms::AbstractMatrix{Complex{NF}};  # spectral coefficients
                            recompute_legendre::Bool=true,      # saves memory
                            Grid::Type{<:AbstractGrid}=FullGaussianGrid,
                            ) where NF                          # number format NF

    _, mmax = size(alms) .- 1                           # -1 for 0-based degree l, order m
    return SpectralTransform(NF,Grid,mmax,recompute_legendre)
end

"""
    S = SpectralTransform(  map::AbstractGrid;
                            recompute_legendre::Bool=true)

Generator function for a `SpectralTransform` struct based on the size and grid type of
gridded field `map`. Recomputes the Legendre polynomials by default."""
function SpectralTransform( map::AbstractGrid{NF};          # gridded field
                            recompute_legendre::Bool=true,  # saves memory
                            ) where NF                      # number format NF

    Grid = typeof(map)
    trunc = get_truncation(map)
    return SpectralTransform(NF,Grid,trunc,recompute_legendre)
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
                    S::SpectralTransform{NF};                   # precomputed parameters struct
                    unscale_coslat::Bool=false                  # unscale with cos(lat) on the fly?
                    ) where {NF<:AbstractFloat}                 # number format NF

    @unpack nlat, nlons, nlat_half, nfreq_max, order = S
    @unpack cos_colat, sin_colat, lon_offsets = S
    @unpack recompute_legendre, Λ, Λs = S
    @unpack brfft_plans = S

    recompute_legendre && @boundscheck size(alms) == size(Λ) || throw(BoundsError)
    recompute_legendre || @boundscheck size(alms) == size(Λs[1]) || throw(BoundsError)
    lmax, mmax = size(alms) .- 1            # maximum degree l, order m of spherical harmonics

    @boundscheck mmax+1 <= nfreq_max || throw(BoundsError)
    @boundscheck nlat == length(cos_colat) || throw(BoundsError)
    @boundscheck typeof(map) <: S.Grid || throw(BoundsError)
    @boundscheck get_nresolution(map) == S.nresolution || throw(BoundsError)

    # preallocate work arrays
    gn = zeros(Complex{NF}, nfreq_max)      # phase factors for northern latitudes
    gs = zeros(Complex{NF}, nfreq_max)      # phase factors for southern latitudes

    Λw = Legendre.Work(Legendre.λlm!, Λ, Legendre.Scalar(zero(Float64)))

    @inbounds for j_north in 1:nlat_half    # symmetry: loop over northern latitudes only
        j_south = nlat - j_north + 1        # southern latitude index
        nlon = nlons[j_north]               # number of longitudes on this ring
        nfreq  = nlon÷2 + 1                 # linear max Fourier frequency wrt to nlon
        nfreqm = nlon÷(order+1) + 1         # (lin/quad/cub) max frequency to shorten loop over m
        not_equator = j_north != j_south    # is the latitude ring not on equator?

        # Recalculate or use precomputed Legendre polynomials Λ
        recompute_legendre && Legendre.unsafe_legendre!(Λw, Λ, lmax, mmax, Float64(cos_colat[j_north]))
        Λj = recompute_legendre ? Λ : Λs[j_north]

        # inverse Legendre transform by looping over wavenumbers l,m
        lm = 1                              # single index for non-zero l,m indices
        @simd for m in 1:min(nfreq,mmax+1)  # Σ_{m=0}^{mmax}, but 1-based index
            acc_odd  = zero(Complex{NF})    # accumulator for isodd(l+m)
            acc_even = zero(Complex{NF})    # accumulator for iseven(l+m)

            # integration over l = m:lmax+1
            lm_end = lm + lmax-m+1              # first index lm plus lmax-m+1 (length of column -1)
            even_degrees = iseven(lm+lm_end)    # is there an even number of degrees in column m?

            # anti-symmetry: sign change of odd harmonics on southern hemisphere
            # but put both into one loop for contiguous memory access
            for lm_even in lm:2:lm_end-even_degrees     
                # split into even, i.e. iseven(l+m)
                # acc_even += alms[lm_even] * Λj[lm_even], but written with muladd
                acc_even = muladd(alms[lm_even],Λj[lm_even],acc_even)

                # and odd (isodd(l+m)) harmonics
                # acc_odd += alms[lm_odd] * Λj[lm_odd], but written with muladd
                acc_odd = muladd(alms[lm_even+1],Λj[lm_even+1],acc_odd)
            end

            # for even number of degrees, one acc_even iteration is skipped, do now
            acc_even = even_degrees ? muladd(alms[lm_end],Λj[lm_end],acc_even) : acc_even

            acc_n = (acc_even + acc_odd)        # accumulators for northern
            acc_s = (acc_even - acc_odd)        # and southern hemisphere
            
            # CORRECT FOR LONGITUDE OFFSETTS
            o = lon_offsets[m,j_north]          # longitude offset rotation            

            gn[m] = muladd(acc_n,o,gn[m])       # accumulate in phase factors for northern
            gs[m] = muladd(acc_s,o,gs[m])       # and southern hemisphere

            lm = lm_end + 1                     # first index of next m column
        end

        if unscale_coslat
            @inbounds cosθ = sin_colat[j_north] # sin(colat) = cos(lat)
            gn ./= cosθ                         # scale in place          
            gs ./= cosθ
        end

        # INVERSE FOURIER TRANSFORM in zonal direction
        brfft_plan = brfft_plans[j_north]       # FFT planned wrt nlon on ring
        ilons = each_index_in_ring(map,j_north) # in-ring indices northern ring
        LinearAlgebra.mul!(view(map.data,ilons),brfft_plan,view(gn,1:nfreq))  # perform FFT

        # southern latitude, don't call redundant 2nd fft if ring is on equator 
        ilons = each_index_in_ring(map,j_south) # in-ring indices southern ring
        not_equator && LinearAlgebra.mul!(view(map.data,ilons),brfft_plan,view(gs,1:nfreq))  # perform FFT

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
    @unpack recompute_legendre, Λ, Λs, solid_angles = S
    @unpack rfft_plans, lon_offsets = S
    
    recompute_legendre && @boundscheck size(alms) == size(Λ) || throw(BoundsError)
    recompute_legendre || @boundscheck size(alms) == size(Λs[1]) || throw(BoundsError)
    lmax, mmax = size(alms) .- 1    # maximum degree l, order m of spherical harmonics

    @boundscheck mmax+1 <= nfreq_max || throw(BoundsError)
    @boundscheck nlat == length(cos_colat) || throw(BoundsError)
    @boundscheck typeof(map) <: S.Grid || throw(BoundsError)
    @boundscheck get_nresolution(map) == S.nresolution || throw(BoundsError)

    # preallocate work warrays
    fn = zeros(Complex{NF},nfreq_max)       # Fourier-transformed northern latitude
    fs = zeros(Complex{NF},nfreq_max)       # Fourier-transformed southern latitude

    # partial sums are accumulated in alms, force zeros initially.
    fill!(alms,0)

    Λw = Legendre.Work(Legendre.λlm!, Λ, Legendre.Scalar(zero(Float64)))

    @inbounds for j_north in 1:nlat_half    # symmetry: loop over northern latitudes only
        j_south = nlat - j_north + 1        # corresponding southern latitude index
        nlon = nlons[j_north]               # number of longitudes on this ring
        nfreq  = nlon÷2 + 1                 # linear max Fourier frequency wrt to nlon
        nfreqm = nlon÷(order+1) + 1         # (lin/quad/cub) max frequency to shorten loop over m
        not_equator = j_north != j_south    # is the latitude ring not on equator?

        # FOURIER TRANSFORM in zonal direction
        rfft_plan = rfft_plans[j_north]         # FFT planned wrt nlon on ring
        ilons = each_index_in_ring(map,j_north) # in-ring indices northern ring
        LinearAlgebra.mul!(view(fn,1:nfreq),rfft_plan,view(map.data,ilons))   # Northern latitude

        ilons = each_index_in_ring(map,j_south) # in-ring indices southern ring
                                                # Southern latitude (don't call FFT on Equator)
                                                # then fill fs with zeros and no changes needed further down
        not_equator ? LinearAlgebra.mul!(view(fs,1:nfreq),rfft_plan,view(map.data,ilons)) : fill!(fs,0)

        # LEGENDRE TRANSFORM in meridional direction
        # Recalculate or use precomputed Legendre polynomials Λ
        recompute_legendre && Legendre.unsafe_legendre!(Λw, Λ, lmax, mmax, Float64(cos_colat[j_north]))
        Λj = recompute_legendre ? Λ : Λs[j_north]
        
        # SOLID ANGLES including quadrature weights (sinθ Δθ) and azimuth (Δϕ) on ring j
        ΔΩ = solid_angles[j_north]                      # = sinθ Δθ Δϕ, solid angle for a grid point

        lm = 1                                          # single index for spherical harmonics
        @simd for m in 1:min(nfreq,mmax+1)              # Σ_{m=0}^{mmax}, but 1-based index

            an, as = fn[m], fs[m]

            # SOLID ANGLE QUADRATURE WEIGHTS and LONGITUDE OFFSET
            o = lon_offsets[m,j_north]                  # longitude offset rotation
            ΔΩ *= conj(o)                               # complex conjugate for rotation back to prime meridian

            # LEGENDRE TRANSFORM
            a_even = (an + as)*ΔΩ                       # sign flip due to anti-symmetry with
            a_odd = (an - as)*ΔΩ                        # odd polynomials 

            # integration over l = m:lmax+1
            lm_end = lm + lmax-m+1                      # first index lm plus lmax-m+1 (length of column -1)
            even_degrees = iseven(lm+lm_end)            # is there an even number of degrees in column m?
            
            # anti-symmetry: sign change of odd harmonics on southern hemisphere
            # but put both into one loop for contiguous memory access
            for lm_even in lm:2:lm_end-even_degrees
                # lm_odd = lm_even+1
                # split into even, i.e. iseven(l+m)
                # alms[lm_even] += a_even * Λj[lm_even]#, but written with muladd
                alms[lm_even] = muladd(a_even,Λj[lm_even],alms[lm_even])
                
                # and odd (isodd(l+m)) harmonics
                # alms[lm_odd] += a_odd * Λj[lm_odd]#, but written with muladd
                alms[lm_even+1] = muladd(a_odd,Λj[lm_even+1],alms[lm_even+1])
            end

            # for even number of degrees, one even iteration is skipped, do now
            alms[lm_end] = even_degrees ? muladd(a_even,Λj[lm_end],alms[lm_end]) : alms[lm_end]

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
                    Grid::Type{<:AbstractGrid}=FullGaussianGrid,
                    ) where NF                          # number format NF

    _, mmax = size(alms) .- 1                           # -1 for 0-based degree l, order m
    S = SpectralTransform(NF,Grid,mmax,recompute_legendre)
    return gridded(alms,S)
end

"""
    map = gridded(  alms::AbstractMatrix,
                    S::SpectralTransform)

Spectral transform (spectral to grid space) from spherical coefficients `alms` to a newly allocated gridded
field `map` with precalculated properties based on the SpectralTransform struct `S`. `alms` is converted to
a `LowerTriangularMatrix` to execute the in-place `gridded!`."""
function gridded(   alms::AbstractMatrix,       # spectral coefficients
                    S::SpectralTransform{NF},   # struct for spectral transform parameters
                    ) where NF                  # number format NF
 
    map = zeros(S.Grid{NF},S.nresolution)               # preallocate output
    almsᴸ = LowerTriangularMatrix{Complex{NF}}(alms)    # drop the upper triangle and convert to NF
    gridded!(map,almsᴸ,S)                               # now execute the in-place version
    return map
end

"""
    alms = spectral(    map::AbstractMatrix;
                        Grid::Type{<:AbstractGrid}=FullGaussianGrid,
                        kwargs...)

Converts `map` to `grid(map)` to execute spectral(map::AbstractGrid;kwargs...)."""
function spectral(  map::AbstractMatrix;            # gridded field
                    Grid::Type{<:AbstractGrid}=FullGaussianGrid,
                    kwargs...)

    return spectral(Grid(map);kwargs...)
end

"""
    alms = spectral(    map::AbstractGrid;
                        Grid::Type{<:AbstractGrid}=FullGaussianGrid,
                        kwargs...)

Converts `map` to `Grid(map)` to execute spectral(map::AbstractGrid;kwargs...)."""
function spectral(  map::AbstractGrid{NF};          # gridded field
                    recompute_legendre::Bool=true,  # saves memory
                    ) where NF                      # number format NF

    Grid = typeof(map)
    trunc = get_truncation(map)
    S = SpectralTransform(NF,Grid,trunc,recompute_legendre)
    return spectral(map,S)
end

"""
    alms = spectral(    map::AbstractMatrix;
                        Grid::Type{<:AbstractGrid}=FullGaussianGrid,
                        kwargs...)

Spectral transform (grid to spectral) `map` to `grid(map)` to execute spectral(map::AbstractGrid;kwargs...)."""
function spectral(  map::AbstractGrid{NF},      # gridded field
                    S::SpectralTransform{NF},   # spectral transform struct
                    ) where NF                  # number format NF

    # always use one more l for consistency with vector quantities
    alms = LowerTriangularMatrix{Complex{NF}}(undef,S.lmax+2,S.mmax+1)
    return spectral!(alms,map,S)                # in-place version
end