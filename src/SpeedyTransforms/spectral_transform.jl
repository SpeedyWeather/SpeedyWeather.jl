"""
SpectralTransform struct that contains all parameters and precomputed arrays
to perform a spectral transform. Fields are
$(TYPEDFIELDS)"""
struct SpectralTransform{
    NF,
    ArrayType,
    VectorType,                 # <: ArrayType{NF, 1},
    VectorComplexType,          # <: ArrayType{Complex{NF}, 1},
    VectorIntType,              # <: ArrayType{Int, 1},
    MatrixComplexType,          # <: ArrayType{Complex{NF}, 2},
    LowerTriangularMatrixType,  # <: LowerTriangularArray{NF, 1, ArrayType{NF}},
    LowerTriangularArrayType,   # <: LowerTriangularArray{NF, 2, ArrayType{NF}},
}
    # GRID
    Grid::Type{<:AbstractGridArray} # grid type used
    nlat_half::Int                  # resolution parameter of grid (# of latitudes on one hemisphere, Eq incl)

    # SPECTRAL RESOLUTION
    lmax::Int               # Maximum degree l=[0, lmax] of spherical harmonics
    mmax::Int               # Maximum order m=[0, l] of spherical harmonics
    nfreq_max::Int          # Maximum (at Equator) number of Fourier frequencies (real FFT)
    m_truncs::VectorIntType # Maximum order m to retain per latitude ring

    # CORRESPONDING GRID SIZE
    nlon_max::Int           # Maximum number of longitude points (at Equator)
    nlons::VectorIntType    # Number of longitude points per ring
    nlat::Int               # Number of latitude rings

    # CORRESPONDING GRID VECTORS
    colat::VectorType               # Gaussian colatitudes (0, π) North to South Pole 
    cos_colat::VectorType           # Cosine of colatitudes
    sin_colat::VectorType           # Sine of colatitudes
    lon_offsets::MatrixComplexType  # Offset of first lon per ring from prime meridian

    # NORMALIZATION
    norm_sphere::NF                         # normalization of the l=0, m=0 mode

    # FFT plans, one plan for each latitude ring
    rfft_plans::Vector{AbstractFFTs.Plan}   # FFT plan for grid to spectral transform
    brfft_plans::Vector{AbstractFFTs.Plan}  # spectral to grid transform (inverse)

    # LEGENDRE POLYNOMIALS
    recompute_legendre::Bool                # Pre or recompute Legendre polynomials
    Λ::AssociatedLegendrePolMatrix{NF}      # Legendre polynomials for one latitude (requires recomputing)
    Λs::Vector{LowerTriangularMatrixType}   # Legendre polynomials for all latitudes (all precomputed)
    
    # SOLID ANGLES ΔΩ FOR QUADRATURE
    # (integration for the Legendre polynomials, extra normalisation of π/nlat included)
    # vector is pole to pole although only northern hemisphere required
    solid_angles::VectorType                # = ΔΩ = sinθ Δθ Δϕ (solid angle of grid point)

    # RECURSION FACTORS
    ϵlms::LowerTriangularMatrixType     # precomputed for meridional gradients grad_y1, grad_y2

    # GRADIENT MATRICES (on unit sphere, no 1/radius-scaling included)
    grad_x ::VectorComplexType          # = i*m but precomputed
    grad_y1::LowerTriangularMatrixType  # precomputed meridional gradient factors, term 1
    grad_y2::LowerTriangularMatrixType  # term 2

    # GRADIENT MATRICES FOR U, V -> Vorticity, Divergence
    grad_y_vordiv1::LowerTriangularMatrixType
    grad_y_vordiv2::LowerTriangularMatrixType

    # GRADIENT MATRICES FOR Vorticity, Divergence -> U, V
    vordiv_to_uv_x::LowerTriangularMatrixType
    vordiv_to_uv1::LowerTriangularMatrixType
    vordiv_to_uv2::LowerTriangularMatrixType

    # EIGENVALUES (on unit sphere, no 1/radius²-scaling included)
    eigenvalues  ::VectorType           # = -l*(l+1), degree l of spherical harmonic
    eigenvalues⁻¹::VectorType           # = -1/(l*(l+1))
end

"""
$(TYPEDSIGNATURES)
Generator function for a SpectralTransform struct. With `NF` the number format,
`Grid` the grid type `<:AbstractGrid` and spectral truncation `lmax, mmax` this function sets up
necessary constants for the spetral transform. Also plans the Fourier transforms, retrieves the colatitudes,
and preallocates the Legendre polynomials (if recompute_legendre == false) and quadrature weights."""
function SpectralTransform(
    ::Type{NF},                                 # Number format NF
    Grid::Type{<:AbstractGridArray},            # type of spatial grid used
    lmax::Int,                                  # Spectral truncation: degrees
    mmax::Int;                                  # Spectral truncation: orders
    ArrayType::Type{<:AbstractArray} = Array,   # Array type used for spectral coefficients
    recompute_legendre::Bool = true,            # re or precompute legendre polynomials?
    legendre_shortcut::Symbol = :linear,        # shorten Legendre loop over order m
    dealiasing::Real=DEFAULT_DEALIASING
) where NF

    Grid = RingGrids.nonparametric_type(Grid)   # always use nonparametric super type

    # RESOLUTION PARAMETERS
    nlat_half = get_nlat_half(mmax, dealiasing) # resolution parameter nlat_half,
                                                # number of latitude rings on one hemisphere incl equator
    nlat = get_nlat(Grid, nlat_half)            # 2nlat_half but one less if grids have odd # of lat rings
    nlon_max = get_nlon_max(Grid, nlat_half)    # number of longitudes around the equator
                                                # number of longitudes per latitude ring (one hemisphere only)
    nlons = [RingGrids.get_nlon_per_ring(Grid, nlat_half, j) for j in 1:nlat_half]
    nfreq_max = nlon_max÷2 + 1                      # maximum number of fourier frequencies (real FFTs)

    # LATITUDE VECTORS (based on Gaussian, equi-angle or HEALPix latitudes)
    colat = RingGrids.get_colat(Grid, nlat_half)    # colatitude in radians                             
    cos_colat = cos.(colat)                         # cos(colat)
    sin_colat = sin.(colat)                         # sin(colat)

    # LEGENDRE SHORTCUT OVER ORDER M (to have linear/quad/cubic truncation per latitude ring)
    if legendre_shortcut in (:linear, :quadratic, :cubic)
        s = (linear=1, quadratic=2, cubic=3)[legendre_shortcut]
        m_truncs = [nlons[j]÷(s+1) + 1 for j in 1:nlat_half]
    elseif legendre_shortcut == :linquad_coslat²
        m_truncs = [ceil(Int, nlons[j]/(2 + sin_colat[j]^2)) for j in 1:nlat_half]
    elseif legendre_shortcut == :lincub_coslat
        m_truncs = [ceil(Int, nlons[j]/(2 + 2sin_colat[j])) for j in 1:nlat_half]
    else
        throw(error("legendre_shortcut=$legendre_shortcut not defined."))
    end
    m_truncs = min.(m_truncs, mmax+1)    # only to mmax in any case (otherwise BoundsError)

    # PLAN THE FFTs
    FFT_package = NF <: Union{Float32, Float64} ? FFTW : GenericFFT
    rfft_plans = [FFT_package.plan_rfft(zeros(NF, nlon)) for nlon in nlons]
    brfft_plans = [FFT_package.plan_brfft(zeros(Complex{NF}, nlon÷2 + 1), nlon) for nlon in nlons]
    
    # NORMALIZATION
    norm_sphere = 2sqrt(π)  # norm_sphere at l=0, m=0 translates to 1s everywhere in grid space

    # LONGITUDE OFFSETS OF FIRST GRID POINT PER RING (0 for full and octahedral grids)
    _, lons = RingGrids.get_colatlons(Grid, nlat_half)
    rings = eachring(Grid, nlat_half)                       # compute ring indices
    lon1s = [lons[rings[j].start] for j in 1:nlat_half]     # pick lons at first index for each ring
    lon_offsets = [cispi(m*lon1/π) for m in 0:mmax, lon1 in lon1s]

    # PREALLOCATE LEGENDRE POLYNOMIALS, +1 for 1-based indexing
    # Legendre polynomials for one latitude
    Λ = AssociatedLegendrePolArray{NF, 2, 1, Vector{NF}}(zeros(LowerTriangularMatrix{NF}, lmax+1, mmax+1))

    # allocate memory in Λs for polynomials at all latitudes or allocate dummy array if precomputed
    # Λs is of size (lmax+1) x (mmax+1) x nlat_half unless recomputed
    # for recomputed only Λ is used, not Λs, create dummy array of 0-size instead
    b = ~recompute_legendre                 # true for precomputed
    Λs = [zeros(LowerTriangularMatrix{NF}, b*(lmax+1), b*(mmax+1)) for _ in 1:nlat_half]

    if recompute_legendre == false          # then precompute all polynomials
        Λtemp = zeros(NF, lmax+1, mmax+1)   # preallocate matrix
        for j in 1:nlat_half                # only one hemisphere due to symmetry
            Legendre.λlm!(Λtemp, lmax, mmax, Float64(cos_colat[j]))   # always precalculate in Float64
            # underflow!(Λtemp, sqrt(floatmin(NF)))
            copyto!(Λs[j], Λtemp)
        end
    end

    # SOLID ANGLES WITH QUADRATURE WEIGHTS (Gaussian, Clenshaw-Curtis, or Riemann depending on grid)
    # solid angles are ΔΩ = sinθ Δθ Δϕ for every grid point with
    # sin(θ)dθ are the quadrature weights approximate the integration over latitudes
    # and sum up to 2 over all latitudes as ∫sin(θ)dθ = 2 over 0...π.
    # Δϕ = 2π/nlon is the azimuth every grid point covers
    solid_angles = get_solid_angles(Grid, nlat_half)

    # RECURSION FACTORS
    ϵlms = get_recursion_factors(lmax+1, mmax)

    # GRADIENTS (on unit sphere, hence 1/radius-scaling is omitted)
    grad_x = [im*m for m in 0:mmax]         # zonal gradient (precomputed currently not used)

    # meridional gradient for scalars (coslat scaling included)
    grad_y1 = zeros(LowerTriangularMatrix, lmax+1, mmax+1)      # term 1, mul with harmonic l-1, m
    grad_y2 = zeros(LowerTriangularMatrix, lmax+1, mmax+1)      # term 2, mul with harmonic l+1, m

    for m in 0:mmax                         # 0-based degree l, order m
        for l in m:lmax           
            grad_y1[l+1, m+1] = -(l-1)*ϵlms[l+1, m+1]
            grad_y2[l+1, m+1] = (l+2)*ϵlms[l+2, m+1]
        end
    end

    # meridional gradient used to get from u, v/coslat to vorticity and divergence
    grad_y_vordiv1 = zeros(LowerTriangularMatrix, lmax+1, mmax+1)   # term 1, mul with harmonic l-1, m
    grad_y_vordiv2 = zeros(LowerTriangularMatrix, lmax+1, mmax+1)   # term 2, mul with harmonic l+1, m

    for m in 0:mmax                         # 0-based degree l, order m
        for l in m:lmax          
            grad_y_vordiv1[l+1, m+1] = (l+1)*ϵlms[l+1, m+1]
            grad_y_vordiv2[l+1, m+1] = l*ϵlms[l+2, m+1]
        end
    end

    # zonal integration (sort of) to get from vorticity and divergence to u, v*coslat
    vordiv_to_uv_x = LowerTriangularMatrix([-m/(l*(l+1)) for l in 0:lmax, m in 0:mmax])
    vordiv_to_uv_x[1, 1] = 0

    # meridional integration (sort of) to get from vorticity and divergence to u, v*coslat
    vordiv_to_uv1 = zeros(LowerTriangularMatrix, lmax+1, mmax+1)    # term 1, to be mul with harmonic l-1, m
    vordiv_to_uv2 = zeros(LowerTriangularMatrix, lmax+1, mmax+1)    # term 2, to be mul with harmonic l+1, m

    for m in 0:mmax                         # 0-based degree l, order m
        for l in m:lmax           
            vordiv_to_uv1[l+1, m+1] = ϵlms[l+1, m+1]/l
            vordiv_to_uv2[l+1, m+1] = ϵlms[l+2, m+1]/(l+1)
        end
    end

    vordiv_to_uv1[1, 1] = 0                 # remove NaN from 0/0

    # EIGENVALUES (on unit sphere, hence 1/radius²-scaling is omitted)
    eigenvalues = [-l*(l+1) for l in 0:mmax+1]
    eigenvalues⁻¹ = inv.(eigenvalues)
    eigenvalues⁻¹[1] = 0                    # set the integration constant to 0

    # guarantee a nonparametric type to construct lower triangular types correctly
    ArrayType_ = RingGrids.nonparametric_type(ArrayType)
    return SpectralTransform{
        NF,
        ArrayType_,
        ArrayType_{NF, 1},
        ArrayType_{Complex{NF}, 1},
        ArrayType_{Int, 1},
        ArrayType_{Complex{NF}, 2},
        LowerTriangularArray{NF, 1, ArrayType_{NF, 1}},
        LowerTriangularArray{NF, 2, ArrayType_{NF, 2}},
    }(
        Grid, nlat_half,
        lmax, mmax, nfreq_max, m_truncs,
        nlon_max, nlons, nlat,
        colat, cos_colat, sin_colat, lon_offsets,
        norm_sphere,
        rfft_plans, brfft_plans,
        recompute_legendre, Λ, Λs, solid_angles,
        ϵlms, grad_x, grad_y1, grad_y2,
        grad_y_vordiv1, grad_y_vordiv2, vordiv_to_uv_x,
        vordiv_to_uv1, vordiv_to_uv2,
        eigenvalues, eigenvalues⁻¹
    )
end

"""
$(TYPEDSIGNATURES)
Generator function for a `SpectralTransform` struct based on the size of the spectral
coefficients `alms` and the grid `Grid`. Recomputes the Legendre polynomials by default."""
function SpectralTransform(
    alms::LowerTriangularArray{NF, N, ArrayType};   # spectral coefficients
    recompute_legendre::Bool = true,                # saves memory
    Grid::Type{<:AbstractGrid} = DEFAULT_GRID,
    dealiasing::Real = DEFAULT_DEALIASING,
) where {NF, N, ArrayType}                          # number format NF (can be complex)
    lmax, mmax = size(alms, ZeroBased, as=Matrix)   # 0-based degree l, order m
    return SpectralTransform(real(NF), Grid, lmax, mmax; ArrayType, recompute_legendre, dealiasing)
end

"""
$(TYPEDSIGNATURES)
Generator function for a `SpectralTransform` struct based on the size and grid type of
gridded field `grids`. Recomputes the Legendre polynomials by default."""
function SpectralTransform(
    grids::AbstractGridArray{NF, N, ArrayType};     # gridded field
    recompute_legendre::Bool=true,                  # saves memory as transform is garbage collected anyway
    one_more_degree::Bool=false,                    # returns a square LowerTriangularMatrix by default
    dealiasing::Real = DEFAULT_DEALIASING,
) where {NF, N, ArrayType}                          # number format NF
    Grid = RingGrids.nonparametric_type(typeof(grids))
    trunc = get_truncation(grids, dealiasing)
    return SpectralTransform(NF, Grid, trunc+one_more_degree, trunc; ArrayType, recompute_legendre, dealiasing)
end

# CHECK MATCHING SIZES
# TODO use dispatch to return false if array types don't match?
function ismatching(S::SpectralTransform, L::LowerTriangularArray)
    return (S.lmax, S.mmax) == size(L, ZeroBased, as=Matrix)[1:2]
end

function ismatching(S::SpectralTransform, grid::AbstractGridArray)
    match = S.Grid == RingGrids.nonparametric_type(typeof(grid)) && S.nlat_half == grid.nlat_half
    return match
end

# make `ismatching` commutative
ismatching(L::LowerTriangularArray, S::SpectralTransform) = ismatching(S, L)
ismatching(G::AbstractGridArray,    S::SpectralTransform) = ismatching(S, G)

function Base.DimensionMismatch(S::SpectralTransform, L::LowerTriangularArray)
    s = "SpectralTransform for $(S.lmax+1)x$(S.mmax+1)xN LowerTriangularArrays and $(Base.dims2string(size(L, as=Matrix))) "*
        "LowerTriangularArray do not match."
    return DimensionMismatch(s)
end

function Base.DimensionMismatch(S::SpectralTransform{NF1}, G::AbstractGridArray{NF2}) where {NF1, NF2}
    s = "SpectralTransform{$NF1}($(S.Grid), nlat_half=$(S.nlat_half)) and "*
        "$(RingGrids.nonparametric_type(G)){$NF2} with nlat_half=$(G.nlat_half) do not match."
    return DimensionMismatch(s)
end

"""
$(TYPEDSIGNATURES)
Recursion factors `ϵ` as a function of degree `l` and order `m` (0-based) of the spherical harmonics.
ϵ(l, m) = sqrt((l^2-m^2)/(4*l^2-1)) and then converted to number format NF."""
function ϵlm(::Type{NF}, l::Int, m::Int) where NF
    return convert(NF, sqrt((l^2-m^2)/(4*l^2-1)))
end

"""
$(TYPEDSIGNATURES)
Recursion factors `ϵ` as a function of degree `l` and order `m` (0-based) of the spherical harmonics.
ϵ(l, m) = sqrt((l^2-m^2)/(4*l^2-1)) with default number format Float64."""
ϵlm(l::Int, m::Int) = ϵlm(Float64, l, m)

"""
$(TYPEDSIGNATURES)     
Returns a matrix of recursion factors `ϵ` up to degree `lmax` and order `mmax` of number format `NF`."""
function get_recursion_factors( ::Type{NF}, # number format NF
                                lmax::Int,  # max degree l of spherical harmonics (0-based here)
                                mmax::Int   # max order m of spherical harmonics
                                ) where {NF<:AbstractFloat}

    # +1 for 0-based to 1-based
    ϵlms = zeros(LowerTriangularMatrix{NF}, lmax+1, mmax+1)      
    for m in 1:mmax+1                   # loop over 1-based l, m
        for l in m:lmax+1
            ϵlms[l, m] = ϵlm(NF, l-1, m-1) # convert to 0-based l, m for function call
        end
    end
    return ϵlms
end

# if number format not provided use Float64
get_recursion_factors(lmax::Int, mmax::Int) = get_recursion_factors(Float64, lmax, mmax)

"""$(TYPEDSIGNATURES)
Spectral transform (spectral to grid space) from n-dimensional array `specs` of spherical harmonic
coefficients to an n-dimensional array `grids` of ring grids. Uses FFT in the zonal direction,
and a Legendre Transform in the meridional direction exploiting symmetries. The spectral transform is
number format-flexible but `grids` and the spectral transform `S` have to have the same number format.
Uses the precalculated arrays, FFT plans and other constants in the SpectralTransform struct `S`.
The spectral transform is grid-flexible as long as the `typeof(grids)<:AbstractGridArray` and `S.Grid`
matches."""
function transform!(                # SPECTRAL TO GRID
    grids::AbstractGridArray,       # gridded output
    specs::LowerTriangularArray,    # spectral coefficients input
    S::SpectralTransform{NF};       # precomputed transform
    unscale_coslat::Bool=false,     # unscale with cos(lat) on the fly?
) where NF                          # number format NF

    (; nlat, nlons, nlat_half, nfreq_max ) = S
    (; cos_colat, sin_colat, lon_offsets ) = S
    (; recompute_legendre, Λ, Λs, m_truncs ) = S
    (; brfft_plans ) = S

    @boundscheck ismatching(S, grids) || throw(DimensionMismatch(S, grids))
    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))

    lmax = specs.m - 1            # 0-based maximum degree l of spherical harmonics
    mmax = specs.n - 1            # 0-based maximum order m of spherical harmonics

    # preallocate work arrays
    gn = zeros(Complex{NF}, nfreq_max)      # phase factors for northern latitudes
    gs = zeros(Complex{NF}, nfreq_max)      # phase factors for southern latitudes

    rings = eachring(grids)                 # precomputed ring indices

    Λw = Legendre.Work(Legendre.λlm!, Λ, Legendre.Scalar(zero(Float64)))

    # loop over all specs/grids (e.g. vertical dimension)
    # k, k_grid only differ when specs/grids have a singleton dimension
    @inbounds for (k, k_grid) in zip(eachmatrix(specs), eachgrid(grids))
        for j_north in 1:nlat_half              # symmetry: loop over northern latitudes only
            j_south = nlat - j_north + 1        # southern latitude index
            nlon = nlons[j_north]               # number of longitudes on this ring
            nfreq  = nlon÷2 + 1                 # linear max Fourier frequency wrt to nlon
            m_trunc = m_truncs[j_north]         # (lin/quad/cub) max frequency to shorten loop over m
            not_equator = j_north != j_south    # is the latitude ring not on equator?

            # Recalculate or use precomputed Legendre polynomials Λ
            recompute_legendre && Legendre.unsafe_legendre!(Λw, Λ, lmax, mmax, Float64(cos_colat[j_north]))
            Λj = recompute_legendre ? Λ : Λs[j_north]

            # INVERSE LEGENDRE TRANSFORM by looping over wavenumbers l, m
            lm = 1                              # single index for non-zero l, m indices
            for m in 1:m_trunc                  # Σ_{m=0}^{mmax}, but 1-based index, shortened to m_trunc
                acc_odd  = zero(Complex{NF})    # accumulator for isodd(l+m)
                acc_even = zero(Complex{NF})    # accumulator for iseven(l+m)

                # integration over l = m:lmax+1
                lm_end = lm + lmax-m+1              # first index lm plus lmax-m+1 (length of column -1)
                even_degrees = iseven(lm+lm_end)    # is there an even number of degrees in column m?

                # anti-symmetry: sign change of odd harmonics on southern hemisphere
                # but put both into one loop for contiguous memory access
                for lm_even in lm:2:lm_end-even_degrees     
                    # split into even, i.e. iseven(l+m)
                    # acc_even += specs[lm_even] * Λj[lm_even], but written with muladd
                    acc_even = muladd(specs[lm_even, k], Λj[lm_even], acc_even)

                    # and odd (isodd(l+m)) harmonics
                    # acc_odd += specs[lm_odd] * Λj[lm_odd], but written with muladd
                    acc_odd = muladd(specs[lm_even+1, k], Λj[lm_even+1], acc_odd)
                end

                # for even number of degrees, one acc_even iteration is skipped, do now
                acc_even = even_degrees ? muladd(specs[lm_end, k], Λj[lm_end], acc_even) : acc_even

                acc_n = (acc_even + acc_odd)        # accumulators for northern
                acc_s = (acc_even - acc_odd)        # and southern hemisphere
                
                # CORRECT FOR LONGITUDE OFFSETTS
                o = lon_offsets[m, j_north]         # longitude offset rotation            

                gn[m] = muladd(acc_n, o, gn[m])     # accumulate in phase factors for northern
                gs[m] = muladd(acc_s, o, gs[m])     # and southern hemisphere

                lm = lm_end + 1                     # first index of next m column
            end

            if unscale_coslat
                @inbounds cosθ = sin_colat[j_north] # sin(colat) = cos(lat)
                gn ./= cosθ                         # scale in place          
                gs ./= cosθ
            end

            # INVERSE FOURIER TRANSFORM in zonal direction
            brfft_plan = brfft_plans[j_north]       # FFT planned wrt nlon on ring
            ilons = rings[j_north]                  # in-ring indices northern ring
            LinearAlgebra.mul!(view(grids.data, ilons, k_grid), brfft_plan, view(gn, 1:nfreq))  # perform FFT

            # southern latitude, don't call redundant 2nd fft if ring is on equator 
            ilons = rings[j_south]                  # in-ring indices southern ring
            not_equator && LinearAlgebra.mul!(view(grids.data, ilons, k_grid), brfft_plan, view(gs, 1:nfreq))  # perform FFT

            fill!(gn, zero(Complex{NF}))            # set phase factors back to zero
            fill!(gs, zero(Complex{NF}))
        end
    end

    return grids
end

"""$(TYPEDSIGNATURES)
Spectral transform (grid to spectral space) from n-dimensional array of `grids` to an n-dimensional
array `specs` of spherical harmonic coefficients. Uses FFT in the zonal direction,
and a Legendre Transform in the meridional direction exploiting symmetries. The spectral transform is
number format-flexible but `grids` and the spectral transform `S` have to have the same number format.
Uses the precalculated arrays, FFT plans and other constants in the SpectralTransform struct `S`.
The spectral transform is grid-flexible as long as the `typeof(grids)<:AbstractGridArray` and `S.Grid`
matches."""
function transform!(                    # grid -> spectral
    specs::LowerTriangularArray,        # output: spectral coefficients
    grids::AbstractGridArray{NF},       # input: gridded values
    S::SpectralTransform{NF}            # precomputed spectral transform
) where NF                              # number format
    
    (; nlat, nlat_half, nlons, nfreq_max, cos_colat ) = S
    (; recompute_legendre, Λ, Λs, solid_angles ) = S
    (; rfft_plans, lon_offsets, m_truncs ) = S
    
    @boundscheck ismatching(S, grids) || throw(DimensionMismatch(S, grids))
    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))

    lmax = specs.m - 1            # 0-based maximum degree l of spherical harmonics
    mmax = specs.n - 1            # 0-based maximum order m of spherical harmonics

    # preallocate work arrays
    fn = zeros(Complex{NF}, nfreq_max)      # Fourier-transformed northern latitude
    fs = zeros(Complex{NF}, nfreq_max)      # Fourier-transformed southern latitude

    rings = eachring(grids)                 # precompute ring indices

    # partial sums are accumulated in specs, force zeros initially.
    fill!(specs, 0)

    Λw = Legendre.Work(Legendre.λlm!, Λ, Legendre.Scalar(zero(Float64)))

    # loop over all specs/grids (e.g. vertical dimension)
    # k, k_grid only differ when specs/grids have a singleton dimension
    @inbounds for (k, k_grid) in zip(eachmatrix(specs), eachgrid(grids))
        for j_north in 1:nlat_half              # symmetry: loop over northern latitudes only
            j_south = nlat - j_north + 1        # corresponding southern latitude index
            nlon = nlons[j_north]               # number of longitudes on this ring
            nfreq  = nlon÷2 + 1                 # linear max Fourier frequency wrt to nlon
            m_trunc = m_truncs[j_north]         # (lin/quad/cub) max frequency to shorten loop over m
            not_equator = j_north != j_south    # is the latitude ring not on equator?

            # FOURIER TRANSFORM in zonal direction
            rfft_plan = rfft_plans[j_north]         # FFT planned wrt nlon on ring
            ilons = rings[j_north]                  # in-ring indices northern ring
            LinearAlgebra.mul!(view(fn, 1:nfreq), rfft_plan, view(grids.data, ilons, k_grid))   # Northern latitude

            ilons = rings[j_south]                  # in-ring indices southern ring
                                                    # Southern latitude (don't call FFT on Equator)
                                                    # then fill fs with zeros and no changes needed further down
            not_equator ? LinearAlgebra.mul!(view(fs, 1:nfreq), rfft_plan, view(grids.data, ilons, k_grid)) : fill!(fs, 0)

            # LEGENDRE TRANSFORM in meridional direction
            # Recalculate or use precomputed Legendre polynomials Λ
            recompute_legendre && Legendre.unsafe_legendre!(Λw, Λ, lmax, mmax, Float64(cos_colat[j_north]))
            Λj = recompute_legendre ? Λ : Λs[j_north]
            
            # SOLID ANGLES including quadrature weights (sinθ Δθ) and azimuth (Δϕ) on ring j
            ΔΩ = solid_angles[j_north]                      # = sinθ Δθ Δϕ, solid angle for a grid point

            lm = 1                                          # single index for spherical harmonics
            for m in 1:m_trunc                              # Σ_{m=0}^{mmax}, but 1-based index

                an, as = fn[m], fs[m]

                # SOLID ANGLE QUADRATURE WEIGHTS and LONGITUDE OFFSET
                o = lon_offsets[m, j_north]                 # longitude offset rotation
                ΔΩ_rotated = ΔΩ*conj(o)                     # complex conjugate for rotation back to prime meridian

                # LEGENDRE TRANSFORM
                a_even = (an + as)*ΔΩ_rotated               # sign flip due to anti-symmetry with
                a_odd = (an - as)*ΔΩ_rotated                # odd polynomials 

                # integration over l = m:lmax+1
                lm_end = lm + lmax-m+1                      # first index lm plus lmax-m+1 (length of column -1)
                even_degrees = iseven(lm+lm_end)            # is there an even number of degrees in column m?
                
                # anti-symmetry: sign change of odd harmonics on southern hemisphere
                # but put both into one loop for contiguous memory access
                for lm_even in lm:2:lm_end-even_degrees
                    # lm_odd = lm_even+1
                    # split into even, i.e. iseven(l+m)
                    # specs[lm_even] += a_even * Λj[lm_even]#, but written with muladd
                    specs[lm_even, k] = muladd(a_even, Λj[lm_even], specs[lm_even, k])
                    
                    # and odd (isodd(l+m)) harmonics
                    # specs[lm_odd] += a_odd * Λj[lm_odd]#, but written with muladd
                    specs[lm_even+1, k] = muladd(a_odd, Λj[lm_even+1], specs[lm_even+1, k])
                end

                # for even number of degrees, one even iteration is skipped, do now
                specs[lm_end, k] = even_degrees ? muladd(a_even, Λj[lm_end], specs[lm_end, k]) : specs[lm_end, k]

                lm = lm_end + 1                             # first index of next m column
            end
        end
    end

    return specs
end

# CONVENIENCE/ALLOCATING VERSIONS

"""$(TYPEDSIGNATURES)
Spherical harmonic transform from `grids` to a newly allocated `specs::LowerTriangularArray`
using the precomputed spectral transform `S`."""
function transform(             # GRID TO SPECTRAL
    grids::AbstractGridArray,   # input grid
    S::SpectralTransform{NF},   # precomputed spectral transform
) where NF
    ks = size(grids)[2:end]     # the non-horizontal dimensions
    specs = zeros(LowerTriangularArray{Complex{NF}}, S.lmax+1, S.mmax+1, ks...)
    transform!(specs, grids, S)
    return specs
end

"""$(TYPEDSIGNATURES)
Spherical harmonic transform from `specs` to a newly allocated `grids::AbstractGridArray`
using the precomputed spectral transform `S`."""
function transform(                 # SPECTRAL TO GRID
    specs::LowerTriangularArray,    # input spectral coefficients
    S::SpectralTransform{NF};       # precomputed spectral transform
    kwargs...                       # pass on unscale_coslat=true/false(default)
) where NF
    ks = size(specs)[2:end]     # the non-horizontal dimensions
    grids = zeros(S.Grid{NF}, S.nlat_half, ks...)
    transform!(grids, specs, S; kwargs...)
    return grids
end

"""
$(TYPEDSIGNATURES)
Spectral transform (spectral to grid space) from spherical coefficients `alms` to a newly allocated gridded
field `map`. Based on the size of `alms` the grid type `grid`, the spatial resolution is retrieved based
on the truncation defined for `grid`. SpectralTransform struct `S` is allocated to execute `transform(alms, S)`."""
function transform(                     # SPECTRAL TO GRID
    specs::LowerTriangularArray{NF};    # spectral coefficients input
    recompute_legendre::Bool = true,    # saves memory
    Grid::Type{<:AbstractGrid} = DEFAULT_GRID,
    dealiasing::Real = DEFAULT_DEALIASING,
    kwargs...                           # pass on unscale_coslat=true/false(default)
) where NF                              # number format NF
    lmax, mmax = size(specs, ZeroBased, as=Matrix)
    S = SpectralTransform(real(NF), Grid, lmax, mmax; recompute_legendre, dealiasing)
    return transform(specs, S; kwargs...)
end

"""
$(TYPEDSIGNATURES)
Spectral transform (grid to spectral space) from `grids` to a newly allocated `LowerTriangularArray`.
Based on the size of `grids` and the keyword `dealiasing` the spectral resolution trunc is
retrieved. SpectralTransform struct `S` is allocated to execute `transform(grids, S)`."""
function transform(
    grids::AbstractGridArray{NF};       # gridded fields
    recompute_legendre::Bool = true,    # saves memory
    one_more_degree::Bool = false,      # for lmax+2 x mmax+1 output size
    dealiasing::Real = DEFAULT_DEALIASING,
) where NF                              # number format NF

    Grid = RingGrids.nonparametric_type(typeof(grids))
    trunc = get_truncation(grids.nlat_half, dealiasing)
    S = SpectralTransform(NF, Grid, trunc+one_more_degree, trunc; recompute_legendre, dealiasing)
    return transform(grids, S)
end