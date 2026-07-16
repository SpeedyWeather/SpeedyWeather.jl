const DEFAULT_NLAYERS = 1
const DEFAULT_GRID = FullGaussianGrid
const DEFAULT_NF = Float32
const DEFAULT_ARRAYTYPE = Array

abstract type AbstractSpectralTransform{NF, AR} end
Architectures.nonparametric_type(S::AbstractSpectralTransform) = nonparametric_type(typeof(S))

"""SpectralTransform struct that contains all parameters and precomputed arrays
to perform a spectral transform. Fields are
$(TYPEDFIELDS)"""
struct SpectralTransform{
        NF,
        AR,                         # <: AbstractArchitecture
        SpectrumType,               # <: AbstractSpectrum
        GridType,                   # <: AbstractGrid
        VectorType,                 # <: ArrayType{NF, 1},
        ArrayTypeIntMatrix,         # <: ArrayType{Int, 2}
        MatrixComplexType,          # <: ArrayType{Complex{NF}, 2},
        LowerTriangularArrayType,   # <: LowerTriangularArray{NF, 2, ArrayType{NF}, ...},
        ScratchType,                # <: ScratchMemory{ArrayComplexType, VectorComplexType},
        GradientType,               # <: NamedTuple for gradients
        IntType,                    # <: Integer
        B,                          # <: Bool
        RFFTSerialType,             # <: Vector{<:AbstractFFTs.Plan}  (K=1, 1-D forward plans)
        BRFFTSerialType,            # <: Vector{<:AbstractFFTs.Plan}  (K=1, 1-D inverse plans)
        RFFTBatchedType,            # <: Dict{Int, <:Vector{<:AbstractFFTs.Plan}}  (K>1, 2-D forward)
        BRFFTBatchedType,           # <: Dict{Int, <:Vector{<:AbstractFFTs.Plan}}  (K>1, 2-D inverse)
    } <: AbstractSpectralTransform{NF, AR}

    # Architecture
    architecture::AR

    # SPECTRAL RESOLUTION
    spectrum::SpectrumType          # spectral truncation
    nfreq_max::IntType              # Maximum (at Equator) number of Fourier frequencies (real FFT)
    LegendreShortcut::Type{<:AbstractLegendreShortcut} # Legendre shortcut for truncation of m loop
    mmax_truncation::Vector{Int}    # Maximum order m to retain per latitude ring

    # GRID
    grid::GridType                  # grid used, including nlat_half for resolution, indices for rings, etc.
    nlayers::IntType                # max number of layers in the vertical (= max planned batch K; for scratch memory size)

    # CORRESPONDING GRID SIZE
    nlon_max::IntType               # Maximum number of longitude points (at Equator)
    nlons::Vector{Int}              # Number of longitude points per ring
    nlat::IntType                   # Number of latitude rings
    rings::Vector{UnitRange{Int}}   # precomputed ring indices

    # CORRESPONDING GRID VECTORS
    coslat::VectorType              # Cosine of latitudes, north to south
    coslat⁻¹::VectorType            # inverse of coslat inv.(coslat)
    lon_offsets::MatrixComplexType  # Offset of first lon per ring from prime meridian

    # NORMALIZATION
    norm_sphere::NF                 # normalization of the l=0, m=0 mode

    # FFT plans, split into concretely-typed serial (K=1, 1-D) and batched (K>1, 2-D) sets so that
    # `plan * view` / `mul!` is a static call (the K=1 and K>1 plans are different concrete FFTW/cuFFT
    # types and cannot share one container). FFTW/cuFFT plans bake K into the plan at construction.
    # The serial plans are the per-layer fallback (used by `_fourier_serial!`), always present. Hot K
    # values (batched dycore transforms, prognostic spec→grid, U/V) are pre-planned in the batched
    # Dicts, keyed by K; everything else falls back to the serial plans in a loop.
    rfft_plan_serial::RFFTSerialType     # grid → spectral (forward), K=1 per-ring 1-D plans
    brfft_plan_serial::BRFFTSerialType   # spectral → grid (inverse), K=1 per-ring 1-D plans
    rfft_plans_batched::RFFTBatchedType  # grid → spectral (forward), K>1 per-ring 2-D plans, keyed by K
    brfft_plans_batched::BRFFTBatchedType # spectral → grid (inverse), K>1 per-ring 2-D plans, keyed by K

    # LEGENDRE POLYNOMIALS, for all latitudes, precomputed
    legendre_polynomials::LowerTriangularArrayType

    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    # state is undetermined, only read after writing to it
    scratch_memory::ScratchType

    jm_index_size::IntType                      # number of indices per layer in kjm_indices
    kjm_indices::ArrayTypeIntMatrix             # precomputed kjm loop indices map for legendre transform

    # SOLID ANGLES ΔΩ FOR QUADRATURE
    # (integration for the Legendre polynomials, extra normalisation of π/nlat included)
    # vector is pole to pole although only northern hemisphere required
    solid_angles::VectorType                    # = ΔΩ = sinθ Δθ Δϕ (solid angle of grid point)

    gradients::GradientType                     # precomputed gradient and integration matrices

    # CUDA GRAPHS
    # toggle for the CUDA-Graphs accelerated batched Fourier transform
    # Set to `false` to fall back to the generic (allocating) per-ring GPU
    # Meaningless for non-CUDA architectures.
    cuda_graphs::B
end

# eltype of a transform is the number format used within
Base.eltype(S::AbstractSpectralTransform{NF}) where {NF} = NF
Architectures.array_type(S::AbstractSpectralTransform{NF, AR}) where {NF, AR} = array_type(AR)
Architectures.nonparametric_type(::Type{<:SpectralTransform}) = SpectralTransform

"""
$(TYPEDSIGNATURES)
Generator function for a SpectralTransform struct. With `NF` the number format,
`Grid` the grid type `<:AbstractGrid` and spectral truncation `lmax, mmax` this function sets up
necessary constants for the spetral transform. Also plans the Fourier transforms, retrieves the colatitudes,
and preallocates the Legendre polynomials and quadrature weights."""
function SpectralTransform(
        spectrum::AbstractSpectrum,                     # Spectral truncation
        grid::AbstractGrid;                             # grid used and resolution, e.g. FullGaussianGrid
        NF::Type{<:Real} = DEFAULT_NF,                                                  # Number format NF
        nlayers::Integer = DEFAULT_NLAYERS,                                             # scratch size — max layer count any single transform call may carry
        transform_batch::AbstractVector{<:Integer} = Int[1, nlayers],                   # list of batch sizes K to pre-plan FFTs for (independent of scratch size)
        LegendreShortcut::Type{<:AbstractLegendreShortcut} = LegendreShortcutLinear,    # shorten Legendre loop over order m
        cuda_graphs::Bool = true,                                            # use CUDA-Graphs accelerated Fourier path (CUDA only)
    )
    # planned_K controls which Ks get pre-built FFT plans. K=1 is always planned (it is the
    # per-layer fallback used by `_fourier_serial!`).
    planned_K = sort!(unique!(Int[1; transform_batch]))
    @assert maximum(planned_K) <= nlayers "transform_batch $planned_K contains an entry larger than nlayers=$nlayers; scratch would be too small."
    (; lmax, mmax, architecture) = spectrum                       # 1-based spectral truncation order and degree

    ArrayType = array_type(architecture)
    ArrayType_ = nonparametric_type(ArrayType)      # drop parameters of ArrayType

    (; nlat_half) = grid                            # number of latitude rings on one hemisphere incl equator

    # RESOLUTION PARAMETERS
    nlat = get_nlat(grid)           # 2nlat_half but one less if grid has odd # of lat rings
    nlon_max = get_nlon_max(grid)   # number of longitudes around the equator
    # number of longitudes per latitude ring (one hemisphere only)
    nlons = [RingGrids.get_nlon_per_ring(grid, j) for j in 1:nlat_half]
    nfreq_max = nlon_max ÷ 2 + 1              # maximum number of fourier frequencies (real FFTs)
    rings = on_architecture(CPU(), eachring(grid))  # precomputed ring indices

    # LATITUDE VECTORS (based on Gaussian, equi-angle or HEALPix latitudes)
    latd = RingGrids.get_latd(grid)         # latitude in degrees (90˚Nto -90˚N)
    colat = RingGrids.get_colat(grid)       # colatitude in radians
    cos_colat = cos.(colat)                 # cos(colat)
    coslat = cosd.(latd)                    # cos(lat)
    coslat⁻¹ = inv.(coslat)                 # 1/cos(lat)

    # LEGENDRE SHORTCUT OVER ORDER M (0-based), truncate the loop over order m
    mmax_truncation = [LegendreShortcut(nlons[j], latd[j]) for j in 1:nlat_half]
    mmax_truncation = min.(mmax_truncation, mmax - 1)   # only to mmax in any case (otherwise BoundsError)

    # NORMALIZATION
    norm_sphere = 2sqrt(π)  # norm_sphere at l=0, m=0 translates to 1s everywhere in grid space

    # LONGITUDE OFFSETS OF FIRST GRID POINT PER RING (0 for full and octahedral grids)
    lons, _ = RingGrids.get_lonlats(grid)
    rings = adapt(Array, eachring(grid))                    # compute ring indices (on CPU)
    lon1s = [lons[rings[j].start] for j in 1:nlat_half]     # pick lons at first index for each ring
    lon_offsets = [cispi(m * lon1 / π) for m in 0:(mmax - 1), lon1 in lon1s]

    # PRECOMPUTE LEGENDRE POLYNOMIALS
    legendre_polynomials = zeros(LowerTriangularArray{NF}, spectrum, nlat_half)
    legendre_polynomials_j = zeros(NF, lmax, mmax)          # temporary for one latitude
    for j in 1:nlat_half                                    # only one hemisphere due to symmetry
        Legendre.λlm!(legendre_polynomials_j, lmax - 1, mmax - 1, cos_colat[j])     # precompute l, m 0-based
        legendre_polynomials[:, j] = LowerTriangularArray(legendre_polynomials_j)   # store
    end

    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    # CPU always chunks unplanned K down to the largest planned batch, so scratch only needs
    # to hold one such chunk. GPU shouldn't chunk and must fit the full input K.
    scratch_memory = ScratchMemory(NF, architecture, grid, _scratch_nlayers(architecture, nlayers, planned_K))

    fake_grid_data = on_architecture(architecture, zeros(NF, grid, nlayers))

    # PLAN THE FFTs: concretely-typed serial (K=1) plan vectors + batched (K>1) plan Dicts
    rfft_plan_serial, brfft_plan_serial, rfft_plans_batched, brfft_plans_batched = plan_FFTs(
        planned_K, fake_grid_data, scratch_memory.north, rings, nlons
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
            for m in 1:(mmax_j + 1)
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

    # PRECOMPUTE GRADIENT AND INTEGRATION MATRICES
    gradients = gradient_arrays(NF, spectrum)

    return SpectralTransform{
        NF,
        typeof(architecture),
        typeof(spectrum),
        typeof(grid),
        array_type(architecture, NF, 1),
        array_type(architecture, Int, 2),
        array_type(architecture, Complex{NF}, 2),
        typeof(legendre_polynomials),
        typeof(scratch_memory),
        typeof(gradients),
        typeof(nlayers),
        typeof(cuda_graphs),
        typeof(rfft_plan_serial),
        typeof(brfft_plan_serial),
        typeof(rfft_plans_batched),
        typeof(brfft_plans_batched),
    }(
        architecture,
        spectrum, nfreq_max,
        LegendreShortcut, mmax_truncation,
        grid, nlayers,
        nlon_max, nlons, nlat, rings,
        coslat, coslat⁻¹, lon_offsets,
        norm_sphere,
        rfft_plan_serial, brfft_plan_serial,
        rfft_plans_batched, brfft_plans_batched,
        legendre_polynomials,
        scratch_memory,
        jm_index_size, kjm_indices,
        solid_angles,
        gradients,
        cuda_graphs,
    )
end

# make commutative
SpectralTransform(grid::AbstractGrid, spectrum::AbstractSpectrum; kwargs...) =
    SpectralTransform(spectrum, grid; kwargs...)

# calculate spectrum if not provided
function SpectralTransform(
        grid::AbstractGrid;
        trunc::Integer = 0,                         # spectral truncation (0-indexed)
        dealiasing::Real = DEFAULT_DEALIASING,      # dealiasing factor
        one_more_degree::Bool = false,              # returns a square LowerTriangularMatrix by default
        kwargs...
    )
    # get trunc from dealiasing if not provided
    trunc = trunc > 0 ? trunc : get_truncation(grid, dealiasing)
    spectrum = Spectrum(trunc + 1 + one_more_degree, trunc + 1; architecture = architecture(grid))
    return SpectralTransform(spectrum, grid; kwargs...)
end

# calculate grid if not provided
function SpectralTransform(
        spectrum::AbstractSpectrum;                     # spectral coefficients
        Grid::Type{<:AbstractGrid} = DEFAULT_GRID,      # grid type, e.g. FullGaussianGrid
        nlat_half::Integer = 0,                         # resolution parameter nlat_half
        dealiasing::Real = DEFAULT_DEALIASING,          # dealiasing factor
        kwargs...
    )
    # get nlat_half from dealiasing if not provided
    nlat_half = nlat_half > 0 ? nlat_half : get_nlat_half(spectrum.mmax - 1, dealiasing)
    grid = Grid(nlat_half, spectrum.architecture)                  # create grid with nlat_half
    return SpectralTransform(spectrum, grid; kwargs...)
end

# use spectrum, NF, and ArrayType from coeffs
function SpectralTransform(
        coeffs::LowerTriangularArray;                   # spectral coefficients
        kwargs...
    )
    (; spectrum) = coeffs
    NF = real(eltype(coeffs))
    nlayers = size(coeffs, 2)
    return SpectralTransform(spectrum; NF, nlayers, kwargs...)
end

# use grid, NF, and ArrayType from field
function SpectralTransform(
        field::AbstractField;                       # gridded field
        kwargs...
    )
    (; grid) = field
    NF = eltype(field)                          # number format of the spectral coefficients
    nlayers = size(field, 2)
    return SpectralTransform(grid; NF, nlayers, kwargs...)
end

"""$(TYPEDSIGNATURES)
Generator function for a `SpectralTransform` struct to transform between `field` and `coeffs`."""
function SpectralTransform(
        field::AbstractField,
        coeffs::LowerTriangularArray;
        kwargs...
    )
    @boundscheck ismatching(architecture(field), architecture(coeffs)) || throw(ArgumentError("Architectures of field and coeffs do not match."))

    # infer types for SpectralTransform
    NF = promote_type(real(eltype(field)), real(eltype(coeffs)))

    # get resolution
    (; spectrum) = coeffs
    (; grid) = field
    nlayers = size(field, 2)
    @boundscheck nlayers == size(coeffs, 2) || throw(DimensionMismatch("Number of layers in field ($nlayers) and coeffs ($(size(coeffs, 2))) do not match."))
    return SpectralTransform(spectrum, grid; NF, nlayers, kwargs...)
end

# make commutative
SpectralTransform(coeffs::LowerTriangularArray, field::AbstractField; kwargs...) =
    SpectralTransform(field, coeffs; kwargs...)

# CHECK MATCHING SIZES
"""$(TYPEDSIGNATURES)
Spectral transform `S` and lower triangular matrix `L` match if the
spectral dimensions `(lmax, mmax)` match and the number of vertical layers is
equal or larger in the transform (constraints due to allocated scratch memory size)."""
function Architectures.ismatching(S::AbstractSpectralTransform, L::LowerTriangularArray; horizontal_only::Bool = false)
    resolution_match = resolution(S.spectrum) == size(L, OneBased, as = Matrix)[1:2]
    vertical_match = horizontal_only ? true : length(axes(L, 2)) <= S.nlayers
    return resolution_match && vertical_match
end

"""$(TYPEDSIGNATURES)
Spectral transform `S` and `grid` match if the resolution `nlat_half` and the
type of the grid match and the number of vertical layers is equal or larger in
the transform (constraints due to allocated scratch memory size)."""
function Architectures.ismatching(S::AbstractSpectralTransform, field::AbstractField; horizontal_only::Bool = false)
    grid_match = S.grid == field.grid   # grid type and resolution match
    vertical_match = horizontal_only ? true : size(field, 2) <= S.nlayers
    return grid_match && vertical_match
end

# make `ismatching` commutative
Architectures.ismatching(L::LowerTriangularArray, S::AbstractSpectralTransform) = ismatching(S, L)
Architectures.ismatching(F::AbstractField, S::AbstractSpectralTransform) = ismatching(S, F)

function Base.DimensionMismatch(S::AbstractSpectralTransform, L::LowerTriangularArray)
    S_ = nonparametric_type(S)
    s = "$S_ for $(S.spectrum.lmax)x$(S.spectrum.mmax)x$(S.nlayers) LowerTriangularArrays " *
        "with $(Base.dims2string(size(L, as = Matrix))) " *
        "LowerTriangularArray do not match."
    return DimensionMismatch(s)
end

function Base.DimensionMismatch(S::AbstractSpectralTransform, F::AbstractField)
    S_ = nonparametric_type(S)
    sz = (RingGrids.get_npoints(S.grid), S.nlayers)
    Grid = nonparametric_type(S.grid)
    Grid2 = nonparametric_type(F.grid)
    s = "$S_ for $(Base.dims2string(sz)) $Grid with " *
        "$(Base.dims2string(size(F))) field on $Grid2 do not match."
    return DimensionMismatch(s)
end

# Layer size of scratch memory (north/south columns in `ScratchMemory`).
# On CPU we typically use `nlayers` only (largest allocated plan and physical layers of the model)
# On GPU we batch transforms much more, and need to use a scratch memory in 
# accordance with the largest batch which is `nlayers` in this case and ISN'T equal to the physical layers of the model
_scratch_nlayers(::AbstractCPU, ::Integer, planned_K::AbstractVector{<:Integer}) = maximum(planned_K)
_scratch_nlayers(::AbstractArchitecture, nlayers::Integer, ::AbstractVector{<:Integer}) = nlayers

# Return true if a call with `K` layers should be split into chunks that are 
# executed sequentially. 
# Architecture-dispatched: 
# CPU chunks to avoid FFTW's alignment-mismatch on x86 (no problem on ARM)
# GPU doesn't need any of that we have seperate plans for each K 
@inline function _needs_chunking(K::Integer, S::SpectralTransform{NF, <:AbstractCPU}) where {NF}
    K > 1 || return false                            # K=1 always handled directly
    haskey(S.rfft_plans_batched, K) && return false  # K is directly planned, no chunking needed
    # K > 1 and not directly planned: chunk through whichever planned K is largest (≤ K_total).
    # Note this includes the degenerate case `keys(rfft_plans) == [1]`: chunks then have K=1
    # and route to the K=1 plan one layer at a time. Scratch (sized to maximum(planned_K))
    # is guaranteed to fit each chunk in all cases.
    return true
end

@inline _needs_chunking(::Integer, ::SpectralTransform) = false

# Largest planned batch K ≤ K_total, excluding the K=1 fallback. Used to pick the chunk size. Typically this will just be the number of layers in the model. 
@inline function _largest_planned_batch(K_total::Integer, S::SpectralTransform)
    best = 1
    for k in keys(S.rfft_plans_batched)
        if k > best && k <= K_total
            best = k
        end
    end
    return best
end

# For CPU only: Split a transform over `K` layers into chunks that each match a pre-built batched FFT plan.
# Each chunk is a sub-view of the input/output along the layer dim; calling `transform!`
# on the sub-views routes through the regular batched path. The trailing remainder
# (size < K_batched) goes through the K=1 serial path one layer at a time.
#
# Why this matters on CPU: FFT plans bake their batch dim K into the plan, and on x86 the
# scratch column-stride (= nfreq_max * sizeof(Complex{NF})) is generally not a multiple
# of the SIMD width, so the K=1 plan can't be applied to columns other than column 1
# without an alignment-mismatch error. Chunking with the sequential execution sidesteps 
# that path for the bulk of the layers and effectively restores the previous behavior before 
# fusion/batching without performance penalties. 
# Greedy chunk sizes: at each position take the largest planned batch that fits the remaining
# layers (1 = the always-planned serial fallback when nothing larger fits). Reproduces the
# previous recursive re-chunking exactly (e.g. K=13, planned {1,3,8} → 8+3+1+1) but iteratively,
# and each chunk calls the NON-chunked core directly
function _transform_chunked!(                       # SPECTRAL TO GRID
        field::AbstractField, coeffs::LowerTriangularArray,
        scratch_memory::ScratchMemory, S::SpectralTransform;
        unscale_coslat::Bool = false,
    )
    K = size(field, 2)
    c = 1
    while c <= K
        len = _largest_planned_batch(K - c + 1, S)  # planned batch (or 1) fitting the remainder
        chunk = c:(c + len - 1)
        field_chunk = wrapped_view(field, :, chunk)
        coeffs_chunk = wrapped_view(coeffs, :, chunk)
        _transform_nonchunked!(field_chunk, coeffs_chunk, scratch_memory, S, unscale_coslat)
        c += len
    end
    return field
end

function _transform_chunked!(                       # GRID TO SPECTRAL
        coeffs::LowerTriangularArray, field::AbstractField,
        scratch_memory::ScratchMemory, S::SpectralTransform,
    )
    K = size(field, 2)
    c = 1
    while c <= K
        len = _largest_planned_batch(K - c + 1, S)  # planned batch (or 1) fitting the remainder
        chunk = c:(c + len - 1)
        field_chunk = wrapped_view(field, :, chunk)
        coeffs_chunk = wrapped_view(coeffs, :, chunk)
        _transform_nonchunked!(coeffs_chunk, field_chunk, scratch_memory, S)
        c += len
    end
    return coeffs
end

"""$(TYPEDSIGNATURES)
Spectral transform (spectral to grid space) from n-dimensional array `coeffs` of spherical harmonic
coefficients to an n-dimensional array `field`. Uses FFT in the zonal direction,
and a Legendre Transform in the meridional direction exploiting symmetries. The spectral transform is
number format-flexible but `field` and the spectral transform `S` have to have the same number format.
Uses the precalculated arrays, FFT plans and other constants in the SpectralTransform struct `S`.
The spectral transform is grid-flexible as long as `field.grid` and `S.grid` match."""
transform!(                                 # SPECTRAL TO GRID
    field::AbstractField,                   # gridded output
    coeffs::LowerTriangularArray,           # spectral coefficients input
    S::AbstractSpectralTransform;           # precomputed transform
    unscale_coslat::Bool = false,           # unscale with cos(lat) on the fly?
) = transform!(field, coeffs, S.scratch_memory, S; unscale_coslat)

function transform!(                        # SPECTRAL TO GRID
        field::AbstractField,               # gridded output
        coeffs::LowerTriangularArray,       # spectral coefficients input
        scratch_memory::ScratchMemory,      # explicit scratch memory to use
        S::SpectralTransform;               # precomputed transform
        unscale_coslat::Bool = false,       # unscale with cos(lat) on the fly?
    )
    # thin forwarder to the positional core (see `_transform_grid!`)
    return _transform_grid!(field, coeffs, scratch_memory, S, unscale_coslat)
end

function _transform_grid!(
        field::AbstractField, coeffs::LowerTriangularArray,
        scratch_memory::ScratchMemory, S::SpectralTransform, unscale_coslat::Bool,
    )
    # On CPU, an unplanned K would route the FFT through `_fourier_serial!`, which on x86
    # crashes with an FFTW alignment mismatch when the scratch column-stride is not a multiple
    # of the SIMD width (consecutive columns have differing alignment). Split into chunks
    # matching the largest available batched plan; each chunk uses the batched FFT path
    # against scratch starting at column 1 (always SIMD-aligned). `_needs_chunking` is
    # statically false on GPU/others, so this branch elides there.
    # Checked **before** the boundscheck below so chunking can handle K > S.nlayers
    K = size(field, 2)
    if _needs_chunking(K, S)
        return _transform_chunked!(field, coeffs, scratch_memory, S; unscale_coslat)
    end
    return _transform_nonchunked!(field, coeffs, scratch_memory, S, unscale_coslat)
end

function _transform_nonchunked!(
        field::AbstractField, coeffs::LowerTriangularArray,
        scratch_memory::ScratchMemory, S::SpectralTransform, unscale_coslat::Bool,
    )
    # catch incorrect sizes early
    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck ismatching(S, coeffs) || throw(DimensionMismatch(S, coeffs))

    # use scratch memory for Legendre but not yet Fourier-transformed data
    g_north = scratch_memory.north    # phase factors for northern latitudes
    g_south = scratch_memory.south    # phase factors for southern latitudes

    # INVERSE LEGENDRE TRANSFORM in meridional direction
    _legendre!(g_north, g_south, coeffs, scratch_memory.column, S; unscale_coslat)

    # INVERSE FOURIER TRANSFORM in zonal direction
    _fourier!(field, g_north, g_south, S)

    return field
end

"""$(TYPEDSIGNATURES)
Spectral transform (grid to spectral space) from n-dimensional array `field` to an n-dimensional
array `coeffs` of spherical harmonic coefficients. Uses FFT in the zonal direction,
and a Legendre Transform in the meridional direction exploiting symmetries. The spectral transform is
number format-flexible but `field` and the spectral transform `S` have to have the same number format.
Uses the precalculated arrays, FFT plans and other constants in the SpectralTransform struct `S`.
The spectral transform is grid-flexible as long as `field.grid` and `S.grid` match."""
transform!(                                             # GRID TO SPECTRAL
    coeffs::LowerTriangularArray,                       # output: spectral coefficients
    field::AbstractField,                               # input: gridded values
    S::AbstractSpectralTransform,                       # precomputed spectral transform
) = transform!(coeffs, field, S.scratch_memory, S)

function transform!(                                    # GRID TO SPECTRAL
        coeffs::LowerTriangularArray,                   # output: spectral coefficients
        field::AbstractField,                           # input: gridded values
        scratch_memory::ScratchMemory,                  # explicit scratch memory to use
        S::SpectralTransform,                           # precomputed spectral transform
    )
    # thin forwarder to the positional core (see `_transform_spec!`)
    return _transform_spec!(coeffs, field, scratch_memory, S)
end

# Positional core of the grid-to-spectral `transform!`; analytic-adjoint AD boundary, see
# `_transform_grid!` and SpeedyTransformsEnzymeExt.
function _transform_spec!(
        coeffs::LowerTriangularArray, field::AbstractField,
        scratch_memory::ScratchMemory, S::SpectralTransform,
    )
    # Same chunking-before-boundscheck pattern as the spec→grid `transform!` above —
    # see comment there.
    K = size(field, 2)
    if _needs_chunking(K, S)
        return _transform_chunked!(coeffs, field, scratch_memory, S)
    end
    return _transform_nonchunked!(coeffs, field, scratch_memory, S)
end

# Non-chunked core, factored out for the same recursion-breaking + AD-boundary reasons as
# `_transform_nonchunked!` above.
function _transform_nonchunked!(
        coeffs::LowerTriangularArray, field::AbstractField,
        scratch_memory::ScratchMemory, S::SpectralTransform,
    )
    # catch incorrect sizes early
    @boundscheck ismatching(S, field) || throw(DimensionMismatch(S, field))
    @boundscheck ismatching(S, coeffs) || throw(DimensionMismatch(S, coeffs))

    # use scratch memory for Fourier but not yet Legendre-transformed data
    f_north = scratch_memory.north    # phase factors for northern latitudes
    f_south = scratch_memory.south    # phase factors for southern latitudes

    # FOURIER TRANSFORM in zonal direction
    _fourier!(f_north, f_south, field, S)

    # LEGENDRE TRANSFORM in meridional direction
    _legendre!(coeffs, f_north, f_south, scratch_memory.column, S)

    return coeffs
end

# CONVENIENCE/ALLOCATING VERSIONS

"""$(TYPEDSIGNATURES)
Spherical harmonic transform from `field` to a newly allocated `coeffs::LowerTriangularArray`
using the precomputed spectral transform `S`."""
function transform(                         # GRID TO SPECTRAL
        field::AbstractField,               # input field
        S::AbstractSpectralTransform,       # precomputed spectral transform
    )
    coeffs = similar(field, S.spectrum, Complex{eltype(S)})
    transform!(coeffs, field, S)
    return coeffs
end

"""$(TYPEDSIGNATURES)
Spherical harmonic transform from `coeffs` to a newly allocated `field::AbstractField`
using the precomputed spectral transform `S`."""
function transform(                     # SPECTRAL TO GRID
        coeffs::LowerTriangularArray,   # input spectral coefficients
        S::AbstractSpectralTransform;   # precomputed spectral transform
        kwargs...                       # pass on unscale_coslat=true/false(default)
    )
    field = similar(coeffs, S.grid, eltype(S))
    transform!(field, coeffs, S; kwargs...)
    return field
end

"""
$(TYPEDSIGNATURES)
Spectral transform (spectral to grid space) from spherical coefficients `alms` to a newly allocated gridded
field `map`. Based on the size of `alms` the grid type `grid`, the spatial resolution is retrieved based
on the truncation defined for `grid`. SpectralTransform struct `S` is allocated to execute `transform(alms, S)`."""
function transform(
        coeffs::LowerTriangularArray;               # SPECTRAL TO GRID
        unscale_coslat::Bool = false,               # separate from kwargs as argument for transform!
        kwargs...                                   # arguments for SpectralTransform constructor
    )
    S = SpectralTransform(coeffs; kwargs...)        # precompute transform
    return transform(coeffs, S; unscale_coslat)     # do the transform
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
