const DEFAULT_NLAYERS = 1
const DEFAULT_GRID = FullGaussianGrid

abstract type AbstractSpectralTransform end

"""SpectralTransform struct that contains all parameters and precomputed arrays
to perform a spectral transform. Fields are
$(TYPEDFIELDS)"""
struct SpectralTransform{
    NF,
    ArrayType,                  # non-parametric array type
    VectorType,                 # <: ArrayType{NF, 1},
    VectorComplexType,          # <: ArrayType{Complex{NF}, 1},
    MatrixComplexType,          # <: ArrayType{Complex{NF}, 2},
    ArrayComplexType,           # <: ArrayType{Complex{NF}, 3},
    LowerTriangularMatrixType,  # <: LowerTriangularArray{NF, 1, ArrayType{NF}},
    LowerTriangularArrayType,   # <: LowerTriangularArray{NF, 2, ArrayType{NF}},
} <: AbstractSpectralTransform
    # GRID
    Grid::Type{<:AbstractGridArray} # grid type used
    nlat_half::Int                  # resolution parameter of grid (# of latitudes on one hemisphere, Eq incl)
    nlayers::Int                    # number of layers in the vertical (for scratch memory size)

    # SPECTRAL RESOLUTION
    lmax::Int                       # Maximum degree l=[0, lmax] of spherical harmonics
    mmax::Int                       # Maximum order m=[0, l] of spherical harmonics
    nfreq_max::Int                  # Maximum (at Equator) number of Fourier frequencies (real FFT)
    LegendreShortcut::Type{<:AbstractLegendreShortcut} # Legendre shortcut for truncation of m loop
    mmax_truncation::Vector{Int}   # Maximum order m to retain per latitude ring

    # CORRESPONDING GRID SIZE
    nlon_max::Int                   # Maximum number of longitude points (at Equator)
    nlons::Vector{Int}              # Number of longitude points per ring
    nlat::Int                       # Number of latitude rings

    # CORRESPONDING GRID VECTORS
    coslat::VectorType              # Cosine of latitudes, north to south
    coslat⁻¹::VectorType            # inverse of coslat inv.(coslat)
    lon_offsets::MatrixComplexType  # Offset of first lon per ring from prime meridian

    # NORMALIZATION
    norm_sphere::NF                         # normalization of the l=0, m=0 mode

    # FFT plans, one plan for each latitude ring, batched in the vertical
    rfft_plans::Vector{AbstractFFTs.Plan}   # FFT plan for grid to spectral transform
    brfft_plans::Vector{AbstractFFTs.Plan}  # spectral to grid transform (inverse)

    # FFT plans, but unbatched
    rfft_plans_1D::Vector{AbstractFFTs.Plan}
    brfft_plans_1D::Vector{AbstractFFTs.Plan}

    # LEGENDRE POLYNOMIALS, for all latitudes, precomputed
    legendre_polynomials::LowerTriangularArrayType
    
    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    # state is undetermined, only read after writing to it
    scratch_memory::ScratchMemory{NF, ArrayComplexType} 
    scratch_memory_grid::VectorType                 # scratch memory with 1-stride for FFT output
    scratch_memory_spec::VectorComplexType
    scratch_memory_column_north::VectorComplexType  # scratch memory for vertically batched Legendre transform
    scratch_memory_column_south::VectorComplexType  # scratch memory for vertically batched Legendre transform

    jm_index_size::Int                              # number of indices per layer in kjm_indices
    kjm_indices::ArrayType                          # precomputed kjm loop indices map

    # SOLID ANGLES ΔΩ FOR QUADRATURE
    # (integration for the Legendre polynomials, extra normalisation of π/nlat included)
    # vector is pole to pole although only northern hemisphere required
    solid_angles::VectorType                # = ΔΩ = sinθ Δθ Δϕ (solid angle of grid point)

    # GRADIENT MATRICES (on unit sphere, no 1/radius-scaling included)
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

# eltype of a transform is the number format used within
Base.eltype(S::SpectralTransform{NF}) where NF = NF

"""
$(TYPEDSIGNATURES)
Generator function for a SpectralTransform struct. With `NF` the number format,
`Grid` the grid type `<:AbstractGrid` and spectral truncation `lmax, mmax` this function sets up
necessary constants for the spetral transform. Also plans the Fourier transforms, retrieves the colatitudes,
and preallocates the Legendre polynomials and quadrature weights."""
function SpectralTransform(
    ::Type{NF},                                     # Number format NF
    lmax::Integer,                                  # Spectral truncation: degrees
    mmax::Integer,                                  # Spectral truncation: orders
    nlat_half::Integer;                             # grid resolution, latitude rings on one hemisphere incl equator
    Grid::Type{<:AbstractGridArray} = DEFAULT_GRID, # type of spatial grid used
    ArrayType::Type{<:AbstractArray} = Array,       # Array type used for spectral coefficients (can be parametric)
    nlayers::Integer = DEFAULT_NLAYERS,             # number of layers in the vertical (for scratch memory size)
    LegendreShortcut::Type{<:AbstractLegendreShortcut} = LegendreShortcutLinear,   # shorten Legendre loop over order m
) where NF

    Grid = RingGrids.nonparametric_type(Grid)               # always use nonparametric concrete type
    ArrayType_ = RingGrids.nonparametric_type(ArrayType)    # drop parameters of ArrayType

    # guarantee a nonparametric type to construct lower triangular types correctly
    ArrayType_ = RingGrids.nonparametric_type(ArrayType)

    # RESOLUTION PARAMETERS
    nlat = get_nlat(Grid, nlat_half)            # 2nlat_half but one less if grids have odd # of lat rings
    nlon_max = get_nlon_max(Grid, nlat_half)    # number of longitudes around the equator
                                                # number of longitudes per latitude ring (one hemisphere only)
    nlons = [RingGrids.get_nlon_per_ring(Grid, nlat_half, j) for j in 1:nlat_half]
    nfreq_max = nlon_max÷2 + 1                      # maximum number of fourier frequencies (real FFTs)

    # LATITUDE VECTORS (based on Gaussian, equi-angle or HEALPix latitudes)
    latd = RingGrids.get_latd(Grid, nlat_half)      # latitude in degrees (90˚Nto -90˚N)
    colat = RingGrids.get_colat(Grid, nlat_half)    # colatitude in radians
    cos_colat = cos.(colat)                         # cos(colat)
    coslat = cosd.(latd)                            # cos(lat)
    coslat⁻¹ = inv.(coslat)                         # 1/cos(lat)

    # LEGENDRE SHORTCUT OVER ORDER M (0-based), truncate the loop over order m 
    mmax_truncation = [LegendreShortcut(nlons[j], latd[j]) for j in 1:nlat_half]
    mmax_truncation = min.(mmax_truncation, mmax)   # only to mmax in any case (otherwise BoundsError)

    # NORMALIZATION
    norm_sphere = 2sqrt(π)  # norm_sphere at l=0, m=0 translates to 1s everywhere in grid space
    
    # LONGITUDE OFFSETS OF FIRST GRID POINT PER RING (0 for full and octahedral grids)
    lons, _ = RingGrids.get_lonlats(Grid, nlat_half)
    rings = eachring(Grid, nlat_half)                       # compute ring indices
    lon1s = [lons[rings[j].start] for j in 1:nlat_half]     # pick lons at first index for each ring
    lon_offsets = [cispi(m*lon1/π) for m in 0:mmax, lon1 in lon1s]
    
    # PRECOMPUTE LEGENDRE POLYNOMIALS, +1 for 1-based indexing
    legendre_polynomials = zeros(LowerTriangularMatrix{NF}, lmax+1, mmax+1, nlat_half)
    legendre_polynomials_j = zeros(NF, lmax+1, mmax+1)      # temporary for one latitude
    for j in 1:nlat_half                                    # only one hemisphere due to symmetry
        Legendre.λlm!(legendre_polynomials_j, lmax, mmax, cos_colat[j])             # precompute
        legendre_polynomials[:, j] = LowerTriangularArray(legendre_polynomials_j)  # store
    end
    
    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    scratch_memory = ScratchMemory(NF, ArrayType_, nfreq_max, nlayers, nlat_half)

    # SCRATCH MEMORY TO 1-STRIDE DATA FOR FFTs
    scratch_memory_grid  = ArrayType_(zeros(NF, nlon_max*nlayers))
    scratch_memory_spec  = ArrayType_(zeros(Complex{NF}, nfreq_max*nlayers))

    # SCRATCH MEMORY COLUMNS FOR VERTICALLY BATCHED LEGENDRE TRANSFORM
    scratch_memory_column_north = ArrayType_(zeros(Complex{NF}, nlayers))
    scratch_memory_column_south = ArrayType_(zeros(Complex{NF}, nlayers))

    rfft_plans = Vector{AbstractFFTs.Plan}(undef, nlat_half)
    brfft_plans = Vector{AbstractFFTs.Plan}(undef, nlat_half)
    rfft_plans_1D = Vector{AbstractFFTs.Plan}(undef, nlat_half)
    brfft_plans_1D = Vector{AbstractFFTs.Plan}(undef, nlat_half)

    fake_grid_data = adapt(ArrayType_, zeros(Grid{NF}, nlat_half, nlayers))

    # PLAN THE FFTs
    plan_FFTs!(
        rfft_plans, brfft_plans, rfft_plans_1D, brfft_plans_1D,
        fake_grid_data, scratch_memory.north, rings, nlons
    )
    
    # PRECOMPUTE KJM INDICES FOR LEGENDRE TRANSFORM (0-based)
    # For GPU it's quicker to precompute the indices for the loops in the 
    # legendre transform and store them in a 3D array rather than computing them 
    # on the fly. We also store the jm_index_size for the loop so we can 
    # truncate to fewer layers if needed. 
    jm_index_size = sum(mmax_truncation .+ 1)
    kjm_indices = zeros(Int, jm_index_size * nlayers, 3)
    i = 0
    for k in 1:nlayers
        for (j, mmax_j) in enumerate(mmax_truncation) 
            for m in 1:mmax_j+1
                i += 1
                kjm_indices[i, :] .= [k, j, m]
            end
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

    return SpectralTransform{
        NF,
        ArrayType_,
        ArrayType_{NF, 1},
        ArrayType_{Complex{NF}, 1},
        ArrayType_{Complex{NF}, 2},
        ArrayType_{Complex{NF}, 3},
        LowerTriangularArray{NF, 1, ArrayType_{NF, 1}},
        LowerTriangularArray{NF, 2, ArrayType_{NF, 2}},
    }(
        Grid, nlat_half, nlayers,
        lmax, mmax, nfreq_max, 
        LegendreShortcut, mmax_truncation,
        nlon_max, nlons, nlat,
        coslat, coslat⁻¹, lon_offsets,
        norm_sphere,
        rfft_plans, brfft_plans, rfft_plans_1D, brfft_plans_1D,
        legendre_polynomials,
        scratch_memory, 
        scratch_memory_grid, scratch_memory_spec,
        scratch_memory_column_north, scratch_memory_column_south,
        jm_index_size, kjm_indices,
        solid_angles, grad_y1, grad_y2,
        grad_y_vordiv1, grad_y_vordiv2, vordiv_to_uv_x,
        vordiv_to_uv1, vordiv_to_uv2,
        eigenvalues, eigenvalues⁻¹
    )
end

"""$(TYPEDSIGNATURES)
Generator function for a `SpectralTransform` struct based on the size of the spectral
coefficients `specs`. Use keyword arguments `nlat_half`, `Grid` or `deliasing` (if `nlat_half`
not provided) to define the grid."""
function SpectralTransform(
    specs::LowerTriangularArray{NF, N, ArrayType};  # spectral coefficients
    nlat_half::Integer = 0,                         # resolution parameter nlat_half
    dealiasing::Real = DEFAULT_DEALIASING,          # dealiasing factor
    kwargs...
) where {NF, N, ArrayType}                          # number format NF (can be complex)
    lmax, mmax = size(specs, ZeroBased, as=Matrix)  # 0-based degree l, order m

    # get nlat_half from dealiasing if not provided
    nlat_half = nlat_half > 0 ? nlat_half : get_nlat_half(mmax, dealiasing)
    nlayers = size(specs, 2)
    return SpectralTransform(real(NF), lmax, mmax, nlat_half; ArrayType, nlayers, kwargs...)
end

"""$(TYPEDSIGNATURES)
Generator function for a `SpectralTransform` struct based on the size and grid type of `grids`.
Use keyword arugments `trunc`, `dealiasing` (ignored if `trunc` is used) or `one_more_degree`
to define the spectral truncation."""
function SpectralTransform(
    grids::AbstractGridArray{NF, N, ArrayType};     # gridded field
    trunc::Integer = 0,                             # spectral truncation
    dealiasing::Real = DEFAULT_DEALIASING,          # dealiasing factor
    one_more_degree::Bool = false,                  # returns a square LowerTriangularMatrix by default
    kwargs...
) where {NF, N, ArrayType}                          # number format NF
    Grid = RingGrids.nonparametric_type(typeof(grids))
    trunc = trunc > 0 ? trunc : get_truncation(grids, dealiasing)
    nlayers = size(grids, 2)
    return SpectralTransform(NF, trunc+one_more_degree, trunc, grids.nlat_half; Grid, ArrayType, nlayers, kwargs...)
end

"""$(TYPEDSIGNATURES)
Generator function for a `SpectralTransform` struct to transform between `grids` and `specs`."""
function SpectralTransform(
    grids::AbstractGridArray{NF1, N, ArrayType1},
    specs::LowerTriangularArray{NF2, N, ArrayType2};
    kwargs...
) where {NF1, NF2, N, ArrayType1, ArrayType2}           # number formats 1 and 2
    
    # infer types for SpectralTransform
    NF = promote_type(real(eltype(grids)), real(eltype(specs)))
    Grid = RingGrids.nonparametric_type(typeof(grids))
    _ArrayType1 = RingGrids.nonparametric_type(ArrayType1)
    _ArrayType2 = RingGrids.nonparametric_type(ArrayType2)
    @assert _ArrayType1 == _ArrayType2 "ArrayTypes of grids ($_ArrayType1) and specs ($_ArrayType2) do not match."

    # get resolution
    lmax, mmax = size(specs, ZeroBased, as=Matrix)      # 0-based degree l, order m
    nlayers = size(grids, 2)
    @assert nlayers == size(specs, 2) "Number of layers in grids ($nlayers) and lower triangular matrices ($(size(specs, 2))) do not match."
    return SpectralTransform(NF, lmax, mmax, grids.nlat_half; Grid, ArrayType=ArrayType1, nlayers, kwargs...)
end

# make commutative
SpectralTransform(specs::LowerTriangularArray, grids::AbstractGridArray) = SpectralTransform(grids, specs)

# CHECK MATCHING SIZES
"""$(TYPEDSIGNATURES)
Spectral transform `S` and lower triangular matrix `L` match if the
spectral dimensions `(lmax, mmax)` match and the number of vertical layers is
equal or larger in the transform (constraints due to allocated scratch memory size)."""
function ismatching(S::SpectralTransform, L::LowerTriangularArray)
    resolution_math = (S.lmax, S.mmax) == size(L, ZeroBased, as=Matrix)[1:2]
    vertical_match = length(axes(L, 2)) <= S.nlayers
    return resolution_math && vertical_match
end

"""$(TYPEDSIGNATURES)
Spectral transform `S` and `grid` match if the resolution `nlat_half` and the
type of the grid match and the number of vertical layers is equal or larger in
the transform (constraints due to allocated scratch memory size)."""
function ismatching(S::SpectralTransform, grid::AbstractGridArray)
    type_match = S.Grid == RingGrids.nonparametric_type(typeof(grid))
    resolution_match = S.nlat_half == grid.nlat_half
    vertical_match = size(grid, 2) <= S.nlayers
    return type_match && resolution_match && vertical_match
end

# make `ismatching` commutative
ismatching(L::LowerTriangularArray, S::SpectralTransform) = ismatching(S, L)
ismatching(G::AbstractGridArray,    S::SpectralTransform) = ismatching(S, G)

function Base.DimensionMismatch(S::SpectralTransform, L::LowerTriangularArray)
    s = "SpectralTransform for $(S.lmax+1)x$(S.mmax+1)x$(S.nlayers) LowerTriangularArrays "*
        "with $(Base.dims2string(size(L, as=Matrix))) "*
        "LowerTriangularArray do not match."
    return DimensionMismatch(s)
end

function Base.DimensionMismatch(S::SpectralTransform, G::AbstractGridArray)
    sz = (RingGrids.get_npoints2D(S.Grid, S.nlat_half), S.nlayers)
    s = "SpectralTransform for $(Base.dims2string(sz)) $(S.Grid) with "*
        "$(Base.dims2string(size(G))) $(RingGrids.nonparametric_type(G)) do not match."
    return DimensionMismatch(s)
end

"""
$(TYPEDSIGNATURES)
Recursion factors `ϵ` as a function of degree `l` and order `m` (0-based) of the spherical harmonics.
ϵ(l, m) = sqrt((l^2-m^2)/(4*l^2-1))."""
recursion_factor(l::Int, m::Int) = sqrt((l^2-m^2)/(4*l^2-1))

"""
$(TYPEDSIGNATURES)     
Returns a matrix of recursion factors `ϵ` up to degree `lmax` and order `mmax` of number format `NF`."""
function get_recursion_factors( ::Type{NF}, # number format NF
                                lmax::Int,  # max degree l of spherical harmonics (0-based here)
                                mmax::Int   # max order m of spherical harmonics
                                ) where NF

    # +1 for 0-based to 1-based
    ϵlms = zeros(LowerTriangularMatrix{NF}, lmax+1, mmax+1)      
    for m in 1:mmax+1                                   # loop over 1-based l, m
        for l in m:lmax+1
            ϵlms[l, m] = recursion_factor(l-1, m-1)     # convert to 0-based l, m for function arguments
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
transform!(                             # SPECTRAL TO GRID
    grids::AbstractGridArray,           # gridded output
    specs::LowerTriangularArray,        # spectral coefficients input
    S::SpectralTransform;               # precomputed transform
    unscale_coslat::Bool = false,       # unscale with cos(lat) on the fly?
) = transform!(grids, specs, S.scratch_memory, S; unscale_coslat)

function transform!(                                # SPECTRAL TO GRID
    grids::AbstractGridArray,                       # gridded output
    specs::LowerTriangularArray,                    # spectral coefficients input
    scratch_memory::ScratchMemory,                  # explicit scratch memory to use
    S::SpectralTransform;                           # precomputed transform
    unscale_coslat::Bool = false,                   # unscale with cos(lat) on the fly?
)
    # catch incorrect sizes early
    @boundscheck ismatching(S, grids) || throw(DimensionMismatch(S, grids))
    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))

    # use scratch memory for Legendre but not yet Fourier-transformed data
    g_north = scratch_memory.north    # phase factors for northern latitudes
    g_south = scratch_memory.south    # phase factors for southern latitudes

    # INVERSE LEGENDRE TRANSFORM in meridional direction
    _legendre!(g_north, g_south, specs, S; unscale_coslat)

    # INVERSE FOURIER TRANSFORM in zonal direction
    _fourier!(grids, g_north, g_south, S)

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
transform!(                             # GRID TO SPECTRAL
    specs::LowerTriangularArray,        # output: spectral coefficients
    grids::AbstractGridArray,           # input: gridded values
    S::SpectralTransform,               # precomputed spectral transform
) = transform!(specs, grids, S.scratch_memory, S)
    
function transform!(                               # GRID TO SPECTRAL
    specs::LowerTriangularArray,                   # output: spectral coefficients
    grids::AbstractGridArray,                      # input: gridded values
    scratch_memory::ScratchMemory,                 # explicit scratch memory to use
    S::SpectralTransform,                          # precomputed spectral transform
)
    # catch incorrect sizes early
    @boundscheck ismatching(S, grids) || throw(DimensionMismatch(S, grids))
    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))

    # use scratch memory for Fourier but not yet Legendre-transformed data
    f_north = scratch_memory.north    # phase factors for northern latitudes
    f_south = scratch_memory.south    # phase factors for southern latitudes

    # FOURIER TRANSFORM in zonal direction
    _fourier!(f_north, f_south, grids, S)    
    
    # LEGENDRE TRANSFORM in meridional direction
    _legendre!(specs, f_north, f_south, S)

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
function transform(
    specs::LowerTriangularArray;                # SPECTRAL TO GRID
    unscale_coslat::Bool = false,               # separate from kwargs as argument for transform!
    kwargs...                                   # arguments for SpectralTransform constructor
)
    S = SpectralTransform(specs; kwargs...)     # precompute transform

    return transform(specs, S; unscale_coslat)  # do the transform
end

"""
$(TYPEDSIGNATURES)
Spectral transform (grid to spectral space) from `grids` to a newly allocated `LowerTriangularArray`.
Based on the size of `grids` and the keyword `dealiasing` the spectral resolution trunc is
retrieved. SpectralTransform struct `S` is allocated to execute `transform(grids, S)`."""
function transform(grids::AbstractGridArray; kwargs...) # GRID TO SPECTRAL
    S = SpectralTransform(grids; kwargs...)             # precompute transform
    return transform(grids, S)                          # do the transform
end