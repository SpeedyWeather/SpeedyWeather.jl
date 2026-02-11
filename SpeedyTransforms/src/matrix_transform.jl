"""MatrixSpectralTransform struct that contains all parameters and precomputed matrices
to perform a spectral transform using dense matrices. Fields are
$(TYPEDFIELDS)"""
struct MatrixSpectralTransform{
        NF,
        AR,                         # <: AbstractArchitecture
        SpectrumType,               # <: AbstractSpectrum
        GridType,                   # <: AbstractGrid
        VectorType,                 # <: ArrayType{NF, 1},
        MatrixType,                 # <: ArrayType{NF, 2},
        MatrixComplexType,          # <: ArrayType{Complex{NF}, 2},
        GradientType,               # <: NamedTuple for gradients
        IntType,                    # <: Integer
    } <: AbstractSpectralTransform{NF, AR}

    # Architecture
    architecture::AR

    # SPECTRAL AND GRID RESOLUTION
    spectrum::SpectrumType              # spectral truncation
    grid::GridType                      # grid used, including nlat_half for resolution, indices for rings, etc.
    nlayers::IntType                    # number of layers in vertical

    # CORRESPONDING GRID VECTORS
    coslat::VectorType                  # Cosine of latitudes, north to south
    coslat⁻¹::VectorType                # inverse of coslat inv.(coslat)

    # NORMALIZATION
    norm_sphere::NF                 # normalization of the l=0, m=0 mode

    # THE ACTUAL TRANSFORM MATRICES for forward = LT(FFT(input)) and backward = IFFT(ILT(input))
    forward::MatrixComplexType          # forward transform matrix
    backward::MatrixComplexType         # backward transform matrix

    backward_real::MatrixType           # real part of backward transform matrix
    backward_imag::MatrixType           # imag part of backward transform matrix

    # SCRATCH MEMORY
    scratch_memory::MatrixType
    scratch_memory_old::MatrixComplexType   # state is undetermined, only read after writing to it

    gradients::GradientType             # precomputed gradient and integration matrices
end

Architectures.nonparametric_type(::Type{<:MatrixSpectralTransform}) = MatrixSpectralTransform

"""
$(TYPEDSIGNATURES)
Generator function for a MatrixSpectralTransform struct. With `NF` the number format,
`Grid` the grid type `<:AbstractGrid` and spectral truncation `lmax, mmax` this function sets up
necessary constants for the spectral transform using dense transformation matrices. Also precomputes
the forward and backward transform matrices, retrieves the colatitudes, and preallocates scratch memory."""
function MatrixSpectralTransform(
        spectrum::AbstractSpectrum,                                                     # Spectral truncation
        grid::AbstractGrid;                                                             # grid used and resolution, e.g. FullGaussianGrid
        NF::Type{<:Real} = DEFAULT_NF,                                                  # Number format NF
        nlayers::Integer = DEFAULT_NLAYERS,                                             # number of layers in the vertical (for scratch memory size)
        LegendreShortcut::Type{<:AbstractLegendreShortcut} = LegendreShortcutLinear,    # shorten Legendre loop over order m
    )
    (; architecture) = spectrum                       # 1-based spectral truncation order and degree

    ArrayType = array_type(architecture)
    ArrayType_ = nonparametric_type(ArrayType)      # drop parameters of ArrayType

    # LATITUDE VECTORS (based on Gaussian, equi-angle or HEALPix latitudes)
    latd = RingGrids.get_latd(grid)                         # latitude in degrees (90˚Nto -90˚N)
    coslat = on_architecture(architecture, cosd.(latd))     # cos(lat)
    coslat⁻¹ = on_architecture(architecture, inv.(coslat))  # 1/cos(lat)

    # Create another SpectralTransform to calculate the transform matrices from (do this on the CPU)
    spectrum_cpu = on_architecture(CPU(), spectrum)
    grid_cpu = on_architecture(CPU(), grid)
    S = SpectralTransform(spectrum_cpu, grid_cpu; NF, nlayers, LegendreShortcut)

    npoints = get_npoints(grid)
    nharmonics = LowerTriangularArrays.nonzeros(spectrum)
    forward = zeros(Complex{NF}, nharmonics, npoints)
    backward = zeros(Complex{NF}, npoints, nharmonics)

    field2D = zeros(NF, grid_cpu)
    coeffs2D = zeros(Complex{NF}, spectrum_cpu)

    progress = ProgressMeter.Progress(length(field2D); dt = 2, desc = "Precalculate matrices:")
    forward_matrix!(forward, S, field2D, coeffs2D, progress)
    backward_matrix!(backward, forward)

    forward = on_architecture(architecture, forward)
    backward = on_architecture(architecture, backward)

    backward_real = real(backward)
    backward_imag = imag(backward)

    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    scratch_memory = on_architecture(architecture, zeros(NF, spectrum, nlayers).data)
    scratch_memory_old = on_architecture(architecture, zeros(Complex{NF}, grid, nlayers).data)

    # PRECOMPUTE GRADIENT AND INTEGRATION MATRICES
    gradients = gradient_arrays(NF, spectrum)

    return MatrixSpectralTransform{
        NF,
        typeof(architecture),
        typeof(spectrum),
        typeof(grid),
        typeof(coslat),
        typeof(backward_real),
        typeof(forward),
        typeof(gradients),
        typeof(nlayers),
    }(
        architecture,
        spectrum, grid, nlayers,
        coslat, coslat⁻¹,
        S.norm_sphere,
        forward,
        backward,
        backward_real,
        backward_imag,
        scratch_memory,
        scratch_memory_old,
        gradients,
    )
end

"""$(TYPEDSIGNATURES)
Compute the forward transform matrix `F` from grid to spectral space using the precomputed
spectral transform `S`. The matrix is computed such that multiplying it by flattened grid
data yields spectral coefficients. This function is not yet implemented."""
function forward_matrix!(F, S::AbstractSpectralTransform, field::AbstractField2D, coeffs::LowerTriangularMatrix, progress = nothing)
    for ij in eachindex(field)
        field .= 0                              # unit vector of input
        GPUArrays.@allowscalar field[ij] = 1
        transform!(coeffs, field, S)            # forward transforms of unit vectors
        F[:, ij] .= coeffs.data                 # are the columns of the transformation matrix F
        isnothing(progress) || ProgressMeter.next!(progress)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Compute the backward transform matrix `B` from spectral to grid space using the precomputed
spectral transform `S`. The matrix is computed such that multiplying it by spectral
coefficients yields flattened grid data. This function is not yet implemented."""
function backward_matrix!(B, F)
    B .= LinearAlgebra.pinv(F)              # Moore Penrose pseudo-inverse
    return B
end

"""$(TYPEDSIGNATURES)
Spectral transform (grid to spectral space) from n-dimensional array `field` to an n-dimensional
array `coeffs` of spherical harmonic coefficients. Uses precomputed dense transform matrices to
perform the transformation. The spectral transform is number format-flexible but `field` and the
spectral transform `M` have to have the same number format. The spectral transform is grid-flexible
as long as `field.grid` and `M.grid` match."""
function transform!(                        # GRID TO SPECTRAL
        coeffs::LowerTriangularArray,       # output: spectral coefficients
        field::AbstractField,               # input: gridded values
        scratch_memory,                     # explicit scratch memory (not used only in spectral to grid)
        M::MatrixSpectralTransform,         # precomputed spectral transform
    )

    # catch incorrect sizes early
    @boundscheck ismatching(M, field, horizontal_only = true) || throw(DimensionMismatch(M, field))
    @boundscheck ismatching(M, coeffs, horizontal_only = true) || throw(DimensionMismatch(M, coeffs))
    # TODO: deactivated temporarily because of Reactant issue
    #@boundscheck size(coeffs, 2) == size(field, 2) || throw(DimensionMismatch(field.data, coeffs.data))
    LinearAlgebra.mul!(coeffs.data, M.forward, field.data)
    return coeffs
end

"""$(TYPEDSIGNATURES)
Spectral transform (spectral to grid space) from n-dimensional array `coeffs` of spherical harmonic
coefficients to an n-dimensional array `field`. Uses precomputed dense transform matrices to perform
the transformation. The spectral transform is number format-flexible but `field` and the spectral
transform `M` have to have the same number format. The spectral transform is grid-flexible as long
as `field.grid` and `M.grid` match."""
function transform_old!(                        # SPECTRAL TO GRID
        field::AbstractField,               # gridded output
        coeffs::LowerTriangularArray,       # spectral coefficients input
        scratch_memory,                     # explicit scratch memory to use
        M::MatrixSpectralTransform;         # precomputed transform
        unscale_coslat::Bool = false,       # unscale with cos(lat) on the fly?
    )

    # catch incorrect sizes early
    @boundscheck ismatching(M, field) || throw(DimensionMismatch(M, field))
    @boundscheck ismatching(M, coeffs) || throw(DimensionMismatch(M, coeffs))

    nlayers = size(coeffs, 2)
    if nlayers < size(scratch_memory, 2)
        # use first n layers of scratch memory
        scratch = view(scratch_memory, :, 1:nlayers)

        # multiply into scratch which is complex typed and then take real part into field
        # imaginary part should be zero but destination is used to store intermediate results
        # explicitly convert to real also for NaN + NaN*im results
        LinearAlgebra.mul!(scratch, M.backward, coeffs.data)
        field.data .= real.(scratch)

    else    # don't use view for 3D transforms with all layers
        # TODO if multiplication yields all real then one could write directly into field.data
        # which would be much faster but if the imaginary part is non-zero this throws an error
        # so for now use the scratch memory as intermediate storage and then copy real part
        LinearAlgebra.mul!(scratch_memory, M.backward, coeffs.data)
        field.data .= real.(scratch_memory)
    end

    unscale_coslat && RingGrids._scale_lat!(field, M.coslat⁻¹)
    return field
end

function transform!(                        # SPECTRAL TO GRID
        field::AbstractField,               # gridded output
        coeffs::LowerTriangularArray,       # spectral coefficients input
        scratch_memory,                     # explicit scratch memory to use
        M::MatrixSpectralTransform;         # precomputed transform
        unscale_coslat::Bool = false,       # unscale with cos(lat) on the fly?
    )

    # catch incorrect sizes early
    @boundscheck ismatching(M, field) || throw(DimensionMismatch(M, field))
    @boundscheck ismatching(M, coeffs) || throw(DimensionMismatch(M, coeffs))

    nlayers = size(coeffs, 2)
    scratch = ndims(coeffs) == 1 ? view(scratch_memory, :, 1) : nlayers < size(scratch_memory, 2) ? view(scratch_memory, :, 1:nlayers) : scratch_memory

    # the result is real-valued, therefore we can split the complex multiplication
    # into two real-valued multiplications
    scratch .= real.(coeffs.data)
    LinearAlgebra.mul!(field.data, M.backward_real, scratch)

    scratch .= imag.(coeffs.data)
    LinearAlgebra.mul!(field.data, M.backward_imag, scratch, -1, 1)

    unscale_coslat && RingGrids._scale_lat!(field, M.coslat⁻¹)
    return field
end

# KERNEL-BASED INVERSE TRANSFORM (spectral to grid) for MatrixSpectralTransform
# This is an alternative to the LinearAlgebra.mul!-based transform! above,
# expressing the same operation as a KernelAbstractions @kernel so that
# Reactant can trace and potentially accelerate it.

@kernel inbounds = true function _inverse_matrix_transform_kernel!(
        field_data,             # output: grid point data, (npoints, nlayers)
        backward_real,          # real part of backward matrix, (npoints, nharmonics)
        backward_imag,          # imag part of backward matrix, (npoints, nharmonics)
        coeffs_data,            # input: spectral coefficients, (nharmonics, nlayers)
        nharmonics,             # number of spectral harmonics (first dim of coeffs_data)
    )
    ij, k = @index(Global, NTuple)

    # dot product: field[ij, k] = Σ_lm B_re[ij, lm] * re(c[lm, k]) - B_im[ij, lm] * im(c[lm, k])
    acc = zero(eltype(field_data))
    for lm in 1:nharmonics
        c = coeffs_data[lm, k]
        acc += backward_real[ij, lm] * real(c) - backward_imag[ij, lm] * imag(c)
    end
    field_data[ij, k] = acc
end

"""$(TYPEDSIGNATURES)
Kernel-based spectral-to-grid transform using `MatrixSpectralTransform`.
Expresses the backward transform as a KernelAbstractions `@kernel` over
`(ij, k)` grid point and layer indices, computing the matrix-vector product
without scratch memory. This is an alternative to the `LinearAlgebra.mul!`-based
`transform!` intended for Reactant compatibility testing."""
function transform_kernel!(                     # SPECTRAL TO GRID (kernel version)
        field::AbstractField,                   # gridded output
        coeffs::LowerTriangularArray,           # spectral coefficients input
        M::MatrixSpectralTransform;             # precomputed transform
        unscale_coslat::Bool = false,           # unscale with cos(lat) on the fly?
    )

    # catch incorrect sizes early
    @boundscheck ismatching(M, field) || throw(DimensionMismatch(M, field))
    @boundscheck ismatching(M, coeffs) || throw(DimensionMismatch(M, coeffs))

    npoints = size(M.backward_real, 1)
    nharmonics = size(M.backward_real, 2)
    nlayers = size(coeffs, 2)

    launch!(
        M.architecture,
        RingGridWorkOrder,
        (npoints, nlayers),
        _inverse_matrix_transform_kernel!,
        field.data,
        M.backward_real,
        M.backward_imag,
        coeffs.data,
        nharmonics,
    )

    unscale_coslat && RingGrids._scale_lat!(field, M.coslat⁻¹)
    return field
end
