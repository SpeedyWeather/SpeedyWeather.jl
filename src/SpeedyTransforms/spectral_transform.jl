const DEFAULT_NLAYERS = 1
const DEFAULT_GRID = FullGaussianGrid
const DEFAULT_NF = Float64
const DEFAULT_ARRAYTYPE = Array

abstract type AbstractSpectralTransform end

"""SpectralTransform struct that contains all parameters and precomputed arrays
to perform a spectral transform. Fields are
$(TYPEDFIELDS)"""
struct SpectralTransform{
    NF,
    ArrayType,                  # non-parametric array type
    SpectrumType,               # <: AbstractSpectrum
    GridType,                   # <: AbstractGrid
    VectorType,                 # <: ArrayType{NF, 1},
    VectorComplexType,          # <: ArrayType{Complex{NF}, 1},
    MatrixComplexType,          # <: ArrayType{Complex{NF}, 2},
    ArrayComplexType,           # <: ArrayType{Complex{NF}, 3},
    LowerTriangularMatrixType,  # <: LowerTriangularArray{NF, 1, ArrayType{NF}, SpectrumType},
    LowerTriangularArrayType,   # <: LowerTriangularArray{NF, 2, ArrayType{NF}, SpectrumType},
} <: AbstractSpectralTransform

    # SPECTRAL RESOLUTION
    spectrum::SpectrumType          # spectral trunction 
    nfreq_max::Int                  # Maximum (at Equator) number of Fourier frequencies (real FFT)
    LegendreShortcut::Type{<:AbstractLegendreShortcut} # Legendre shortcut for truncation of m loop
    mmax_truncation::Vector{Int}    # Maximum order m to retain per latitude ring
    
    # GRID
    grid::GridType              # grid used, including nlat_half for resolution, indices for rings, etc.
    nlayers::Int                # number of layers in the vertical (for scratch memory size)

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
    scratch_memory_north::ArrayComplexType
    scratch_memory_south::ArrayComplexType
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
array_type(S::SpectralTransform{NF, A}) where {NF, A} = A

"""
$(TYPEDSIGNATURES)
Generator function for a SpectralTransform struct. With `NF` the number format,
`Grid` the grid type `<:AbstractGrid` and spectral truncation `lmax, mmax` this function sets up
necessary constants for the spetral transform. Also plans the Fourier transforms, retrieves the colatitudes,
and preallocates the Legendre polynomials and quadrature weights."""
function SpectralTransform(
    spectrum::AbstractSpectrum,                     # Spectral truncation
    grid::AbstractGrid;                             # grid used and resolution, e.g. FullGaussianGrid
    NF::Type{<:Real} = DEFAULT_NF,                  # Number format NF
    ArrayType::Type{<:AbstractArray} = DEFAULT_ARRAYTYPE,   # Array type used for spectral coefficients (can be parametric)
    nlayers::Integer = DEFAULT_NLAYERS,             # number of layers in the vertical (for scratch memory size)
    LegendreShortcut::Type{<:AbstractLegendreShortcut} = LegendreShortcutLinear,   # shorten Legendre loop over order m
)

    (; nlat_half) = grid    # number of latitude rings on one hemisphere incl equator
    ArrayType_ = RingGrids.nonparametric_type(ArrayType)    # drop parameters of ArrayType
    (; lmax, mmax) = spectrum                               # 1-based spectral truncation order and degree
    
    # RESOLUTION PARAMETERS
    nlat = get_nlat(grid)           # 2nlat_half but one less if grids have odd # of lat rings
    nlon_max = get_nlon_max(grid)   # number of longitudes around the equator
                                    # number of longitudes per latitude ring (one hemisphere only)
    nlons = [RingGrids.get_nlon_per_ring(grid, j) for j in 1:nlat_half]
    nfreq_max = nlon_max÷2 + 1              # maximum number of fourier frequencies (real FFTs)

    # LATITUDE VECTORS (based on Gaussian, equi-angle or HEALPix latitudes)
    latd = RingGrids.get_latd(grid)         # latitude in degrees (90˚Nto -90˚N)
    colat = RingGrids.get_colat(grid)       # colatitude in radians
    cos_colat = cos.(colat)                 # cos(colat)
    coslat = cosd.(latd)                    # cos(lat)
    coslat⁻¹ = inv.(coslat)                 # 1/cos(lat)

    # LEGENDRE SHORTCUT OVER ORDER M (0-based), truncate the loop over order m 
    mmax_truncation = [LegendreShortcut(nlons[j], latd[j]) for j in 1:nlat_half]
    mmax_truncation = min.(mmax_truncation, mmax-1)   # only to mmax in any case (otherwise BoundsError)

    # NORMALIZATION
    norm_sphere = 2sqrt(π)  # norm_sphere at l=0, m=0 translates to 1s everywhere in grid space
    
    # LONGITUDE OFFSETS OF FIRST GRID POINT PER RING (0 for full and octahedral grids)
    lons, _ = RingGrids.get_lonlats(grid)
    rings = eachring(grid)                                  # compute ring indices
    lon1s = [lons[rings[j].start] for j in 1:nlat_half]     # pick lons at first index for each ring
    lon_offsets = [cispi(m*lon1/π) for m in 0:mmax-1, lon1 in lon1s]
    
    # PRECOMPUTE LEGENDRE POLYNOMIALS
    legendre_polynomials = zeros(LowerTriangularArray{NF}, spectrum, nlat_half)
    legendre_polynomials_j = zeros(NF, lmax, mmax)          # temporary for one latitude
    for j in 1:nlat_half                                    # only one hemisphere due to symmetry
        Legendre.λlm!(legendre_polynomials_j, lmax - 1, mmax - 1, cos_colat[j])     # precompute l, m 0-based
        legendre_polynomials[:, j] = LowerTriangularArray(legendre_polynomials_j)   # store
    end
    
    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    # convert all to the arraytype (e.g. moves array to GPU)
    scratch_memory_north = ArrayType_(zeros(Complex{NF}, nfreq_max, nlayers, nlat_half))
    scratch_memory_south = ArrayType_(zeros(Complex{NF}, nfreq_max, nlayers, nlat_half))

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

    fake_grid_data = adapt(ArrayType_, zeros(NF, grid, nlayers))

    # PLAN THE FFTs
    plan_FFTs!(
        rfft_plans, brfft_plans, rfft_plans_1D, brfft_plans_1D,
        fake_grid_data, scratch_memory_north, rings, nlons
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
    solid_angles = get_solid_angles(grid)

    # RECURSION FACTORS
    ϵlms = get_recursion_factors(spectrum)

    # GRADIENTS (on unit sphere, hence 1/radius-scaling is omitted)
    # meridional gradient for scalars (coslat scaling included)
    grad_y1 = zeros(LowerTriangularMatrix, spectrum)      # term 1, mul with harmonic l-1, m
    grad_y2 = zeros(LowerTriangularMatrix, spectrum)      # term 2, mul with harmonic l+1, m

    for m in 1:mmax                         # 1-based degree l, order m
        for l in m:lmax           
            grad_y1[l, m] = -(l-2)*ϵlms[l, m]
            grad_y2[l, m] = (l+1)*ϵlms[l+1, m]
        end
    end

    # meridional gradient used to get from u, v/coslat to vorticity and divergence
    grad_y_vordiv1 = zeros(LowerTriangularMatrix, spectrum)   # term 1, mul with harmonic l-1, m
    grad_y_vordiv2 = zeros(LowerTriangularMatrix, spectrum)   # term 2, mul with harmonic l+1, m

    for m in 1:mmax                         # 1-based degree l, order m
        for l in m:lmax          
            grad_y_vordiv1[l, m] = l*ϵlms[l, m]
            grad_y_vordiv2[l, m] = (l-1)*ϵlms[l+1, m]
        end
    end

    # zonal integration (sort of) to get from vorticity and divergence to u, v*coslat
    vordiv_to_uv_x = LowerTriangularMatrix([-m/(l*(l+1)) for l in 0:(lmax-1), m in 0:(mmax-1)])
    vordiv_to_uv_x[1, 1] = 0

    # meridional integration (sort of) to get from vorticity and divergence to u, v*coslat
    vordiv_to_uv1 = zeros(LowerTriangularMatrix, spectrum)    # term 1, to be mul with harmonic l-1, m
    vordiv_to_uv2 = zeros(LowerTriangularMatrix, spectrum)    # term 2, to be mul with harmonic l+1, m

    for m in 1:mmax                         # 1-based degree l, order m
        for l in m:lmax           
            vordiv_to_uv1[l, m] = ϵlms[l, m]/(l-1)
            vordiv_to_uv2[l, m] = ϵlms[l+1, m]/l
        end
    end

    vordiv_to_uv1[1, 1] = 0                 # remove NaN from 0/0

    # EIGENVALUES (on unit sphere, hence 1/radius²-scaling is omitted)
    eigenvalues = [-l*(l+1) for l in 0:mmax]
    eigenvalues⁻¹ = inv.(eigenvalues)
    eigenvalues⁻¹[1] = 0                    # set the integration constant to 0

    return SpectralTransform{
        NF,
        ArrayType_,
        typeof(spectrum),
        typeof(grid),
        ArrayType_{NF, 1},
        ArrayType_{Complex{NF}, 1},
        ArrayType_{Complex{NF}, 2},
        ArrayType_{Complex{NF}, 3},
        LowerTriangularArray{NF, 1, ArrayType_{NF, 1}, typeof(spectrum)},
        LowerTriangularArray{NF, 2, ArrayType_{NF, 2}, typeof(spectrum)},
    }(
        spectrum, nfreq_max, 
        LegendreShortcut, mmax_truncation,
        grid, nlayers,
        nlon_max, nlons, nlat,
        coslat, coslat⁻¹, lon_offsets,
        norm_sphere,
        rfft_plans, brfft_plans, rfft_plans_1D, brfft_plans_1D,
        legendre_polynomials,
        scratch_memory_north, scratch_memory_south,
        scratch_memory_grid, scratch_memory_spec,
        scratch_memory_column_north, scratch_memory_column_south,
        jm_index_size, kjm_indices,
        solid_angles, grad_y1, grad_y2,
        grad_y_vordiv1, grad_y_vordiv2, vordiv_to_uv_x,
        vordiv_to_uv1, vordiv_to_uv2,
        eigenvalues, eigenvalues⁻¹
    )
end

SpectralTransform(spectrum::AbstractSpectrum, grid::AbstractGrid; kwargs...) =
    SpectralTransform(DEFAULT_NF, spectrum, grid; kwargs...)

"""$(TYPEDSIGNATURES)
Generator function for a `SpectralTransform` struct based on the spectral truncation
given by integers `lmax` and `mmax.""" 
SpectralTransform(
    ::Type{NF},                                     # Number format NF
    lmax::Integer,                                  # spectral trunction degree
    mmax::Integer,                                  # spectral truncation order
    vargs...; kwargs...                             # grid resolution, latitude rings on one hemisphere incl equator
) where NF = SpectralTransform(NF, Spectrum(lmax, mmax), vargs...; kwargs...)
    
"""$(TYPEDSIGNATURES)
Generator function for a `SpectralTransform` struct based on the size of the spectral
coefficients `specs`. Use keyword arguments `nlat_half`, `Grid` or `deliasing` (if `nlat_half`
not provided) to define the grid."""
function SpectralTransform(
    specs::LowerTriangularArray;                    # spectral coefficients
    Grid::Type{<:AbstractGrid} = DEFAULT_GRID,      # grid type, e.g. FullGaussianGrid
    nlat_half::Integer = 0,                         # resolution parameter nlat_half
    dealiasing::Real = DEFAULT_DEALIASING,          # dealiasing factor
    kwargs...
)
    (; spectrum) = specs
    NF = eltype(specs)                     # number format of the spectral coefficients
    ArrayType = LowerTriangularArrays.array_type(specs)

    # get nlat_half from dealiasing if not provided
    nlat_half = nlat_half > 0 ? nlat_half : get_nlat_half(spectrum.mmax - 1, dealiasing)
    grid = Grid(nlat_half)                  # create grid with nlat_half
    nlayers = size(specs, 2)
    return SpectralTransform(real(NF), spectrum, grid; ArrayType, nlayers, kwargs...)
end

"""$(TYPEDSIGNATURES)
Generator function for a `SpectralTransform` struct based on the size and grid type of `grids`.
Use keyword arugments `trunc`, `dealiasing` (ignored if `trunc` is used) or `one_more_degree`
to define the spectral truncation. `trunc` is assumed to be zero-indexed (i.e. starting with l=0)"""
function SpectralTransform(
    field::AbstractField;                       # gridded field
    trunc::Integer = 0,                         # spectral truncation (0-indexed)
    dealiasing::Real = DEFAULT_DEALIASING,      # dealiasing factor
    one_more_degree::Bool = false,              # returns a square LowerTriangularMatrix by default
    kwargs...
)
    NF = eltype(field)                          # number format of the spectral coefficients
    ArrayType = RingGrids.array_type(field)

    # get trunc from dealiasing if not provided
    trunc = trunc > 0 ? trunc : get_truncation(field.grid, dealiasing)
    spectrum = Spectrum(trunc+1+one_more_degree, trunc+1)
    nlayers = size(field, 2)
    return SpectralTransform(NF, spectrum, field.grid; ArrayType, nlayers, kwargs...)
end

"""$(TYPEDSIGNATURES)
Generator function for a `SpectralTransform` struct to transform between `field` and `specs`."""
function SpectralTransform(
    field::AbstractField,
    specs::LowerTriangularArray;
    kwargs...
) 
    # infer types for SpectralTransform
    NF = promote_type(real(eltype(field)), real(eltype(specs)))
    ArrayType = RingGrids.nonparametric_type(RingGrids.array_type(field))
    ArrayType2 = LowerTriangularArrays.nonparametric_type(LowerTriangularArrays.array_type(specs))
    @assert ArrayType == ArrayType2 "ArrayTypes of field ($_ArrayType1) and specs ($_ArrayType2) do not match."

    # get resolution
    (; spectrum) = specs
    (; grid) = field 
    nlayers = size(field, 2)
    @assert nlayers == size(specs, 2) "Number of layers in field ($nlayers) and specs ($(size(specs, 2))) do not match."
    return SpectralTransform(NF, spectrum, grid; ArrayType, nlayers, kwargs...)
end

# make commutative
SpectralTransform(specs::LowerTriangularArray, field::AbstractField; kwargs...) =
    SpectralTransform(field, specs; kwargs...)

# CHECK MATCHING SIZES
"""$(TYPEDSIGNATURES)
Spectral transform `S` and lower triangular matrix `L` match if the
spectral dimensions `(lmax, mmax)` match and the number of vertical layers is
equal or larger in the transform (constraints due to allocated scratch memory size)."""
function ismatching(S::SpectralTransform, L::LowerTriangularArray)
    resolution_math = resolution(S.spectrum) == size(L, OneBased, as=Matrix)[1:2]
    vertical_match = length(axes(L, 2)) <= S.nlayers
    return resolution_math && vertical_match
end

"""$(TYPEDSIGNATURES)
Spectral transform `S` and `grid` match if the resolution `nlat_half` and the
type of the grid match and the number of vertical layers is equal or larger in
the transform (constraints due to allocated scratch memory size)."""
function ismatching(S::SpectralTransform, field::AbstractField)
    grid_match = S.grid == field.grid   # grid type and resolution match
    vertical_match = size(field, 2) <= S.nlayers
    return grid_match && vertical_match
end

# make `ismatching` commutative
ismatching(L::LowerTriangularArray, S::SpectralTransform) = ismatching(S, L)
ismatching(F::AbstractField,        S::SpectralTransform) = ismatching(S, F)

function Base.DimensionMismatch(S::SpectralTransform, L::LowerTriangularArray)
    s = "SpectralTransform for $(S.spectrum.lmax)x$(S.spectrum.mmax)x$(S.nlayers) LowerTriangularArrays "*
        "with $(Base.dims2string(size(L, as=Matrix))) "*
        "LowerTriangularArray do not match."
    return DimensionMismatch(s)
end

function Base.DimensionMismatch(S::SpectralTransform, F::AbstractField)
    sz = (RingGrids.get_npoints2D(S.Grid, S.nlat_half), S.nlayers)
    s = "SpectralTransform for $(Base.dims2string(sz)) $(S.Grid) with "*
        "$(Base.dims2string(size(F))) $F do not match."
    return DimensionMismatch(s)
end

"""
$(TYPEDSIGNATURES)
Recursion factors `ϵ` as a function of degree `l` and order `m` (0-based) of the spherical harmonics.
ϵ(l, m) = sqrt((l^2-m^2)/(4*l^2-1))."""
recursion_factor(l::Int, m::Int) = sqrt((l^2-m^2)/(4*l^2-1))

"""
$(TYPEDSIGNATURES)     
Returns a matrix of recursion factors `ϵ` up to degree `lmax` and order `mmax` (1-based) of the `spectrum` in number format `NF`."""
function get_recursion_factors( ::Type{NF}, # number format NF
                                spectrum::Spectrum,
                                ) where NF
    (; lmax, mmax) = spectrum 

    ϵlms = zeros(LowerTriangularMatrix{NF}, lmax+1, mmax)      
    for m in 1:mmax                                     # loop over 1-based l, m
        for l in m:lmax+1
            ϵlms[l, m] = recursion_factor(l-1, m-1)     # convert to 0-based l, m for function arguments
        end
    end
    return ϵlms
end

# if number format not provided use Float64
get_recursion_factors(spectrum::Spectrum) = get_recursion_factors(Float64, spectrum)

"""$(TYPEDSIGNATURES)
Spectral transform (spectral to grid space) from n-dimensional array `specs` of spherical harmonic
coefficients to an n-dimensional array `grids` of ring grids. Uses FFT in the zonal direction,
and a Legendre Transform in the meridional direction exploiting symmetries. The spectral transform is
number format-flexible but `grids` and the spectral transform `S` have to have the same number format.
Uses the precalculated arrays, FFT plans and other constants in the SpectralTransform struct `S`.
The spectral transform is grid-flexible as long as `field::AbstractField` and `S.Grid` match."""
function transform!(                    # SPECTRAL TO GRID
    field::AbstractField,               # gridded output
    specs::LowerTriangularArray,        # spectral coefficients input
    S::SpectralTransform;               # precomputed transform
    unscale_coslat::Bool = false,       # unscale with cos(lat) on the fly?
)
    # catch incorrect sizes early
    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))

    # use scratch memory for Legendre but not yet Fourier-transformed data
    g_north = S.scratch_memory_north    # phase factors for northern latitudes
    g_south = S.scratch_memory_south    # phase factors for southern latitudes

    # INVERSE LEGENDRE TRANSFORM in meridional direction
    _legendre!(g_north, g_south, specs, S; unscale_coslat)

    # INVERSE FOURIER TRANSFORM in zonal direction
    _fourier!(field, g_north, g_south, S)

    return field
end

"""$(TYPEDSIGNATURES)
Spectral transform (grid to spectral space) from n-dimensional array of `grids` to an n-dimensional
array `specs` of spherical harmonic coefficients. Uses FFT in the zonal direction,
and a Legendre Transform in the meridional direction exploiting symmetries. The spectral transform is
number format-flexible but `grids` and the spectral transform `S` have to have the same number format.
Uses the precalculated arrays, FFT plans and other constants in the SpectralTransform struct `S`.
The spectral transform is grid-flexible as long as `field::AbstractField` and `S.Grid` match."""
function transform!(                # GRID TO SPECTRAL
    specs::LowerTriangularArray,    # output: spectral coefficients
    field::AbstractField,           # input: gridded values
    S::SpectralTransform,           # precomputed spectral transform
)
    # catch incorrect sizes early
    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))

    # use scratch memory for Fourier but not yet Legendre-transformed data
    f_north = S.scratch_memory_north    # phase factors for northern latitudes
    f_south = S.scratch_memory_south    # phase factors for southern latitudes

    # FOURIER TRANSFORM in zonal direction
    _fourier!(f_north, f_south, field, S)    
    
    # LEGENDRE TRANSFORM in meridional direction
    _legendre!(specs, f_north, f_south, S)

    return specs
end

# CONVENIENCE/ALLOCATING VERSIONS

"""$(TYPEDSIGNATURES)
Spherical harmonic transform from `grids` to a newly allocated `specs::LowerTriangularArray`
using the precomputed spectral transform `S`."""
function transform(             # GRID TO SPECTRAL
    field::AbstractField,       # input field
    S::SpectralTransform{NF},   # precomputed spectral transform
) where NF
    ks = size(field)[2:end]     # the non-horizontal dimensions
    specs = zeros(LowerTriangularArray{Complex{NF}}, S.spectrum, ks...)
    transform!(specs, field, S)
    return specs
end


"""$(TYPEDSIGNATURES)
Spherical harmonic transform from `specs` to a newly allocated `field::AbstractField`
using the precomputed spectral transform `S`."""
function transform(                 # SPECTRAL TO GRID
    specs::LowerTriangularArray,    # input spectral coefficients
    S::SpectralTransform;           # precomputed spectral transform
    kwargs...                       # pass on unscale_coslat=true/false(default)
)
    ks = size(specs)[2:end]         # the non-horizontal dimensions
    field = zeros(eltype(S), S.grid, ks...)
    transform!(field, specs, S; kwargs...)
    return field
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
Spectral transform (grid to spectral space) from `field` to a newly allocated `LowerTriangularArray`.
Based on the size of `field` and the keyword `dealiasing` the spectral resolution trunc is
retrieved. SpectralTransform struct `S` is allocated to execute `transform(field, S)`."""
function transform(field::AbstractField; kwargs...) # GRID TO SPECTRAL
    S = SpectralTransform(field; kwargs...)         # precompute transform
    return transform(field, S)                      # do the transform
end