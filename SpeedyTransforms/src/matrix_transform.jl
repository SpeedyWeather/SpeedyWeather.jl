"""MatrixSpectralTransform struct that contains all parameters and precomputed matrices
to perform a spectral transform using dense matrices. Fields are
$(TYPEDFIELDS)"""
struct MatrixSpectralTransform{
        NF,
        AR,                         # <: AbstractArchitecture
        ArrayType,                  # non-parametric array type
        SpectrumType,               # <: AbstractSpectrum
        GridType,                   # <: AbstractGrid
        VectorType,                 # <: ArrayType{NF, 1},
        MatrixComplexType,          # <: ArrayType{Complex{NF}, 2},
        GradientType,               # <: NamedTuple for gradients
    } <: AbstractSpectralTransform{NF, AR, ArrayType}

    # Architecture
    architecture::AR

    # SPECTRAL AND GRID RESOLUTION
    spectrum::SpectrumType              # spectral truncation
    grid::GridType                      # grid used, including nlat_half for resolution, indices for rings, etc.
    nlayers::Int                        # number of layers in vertical

    # CORRESPONDING GRID VECTORS
    coslat::VectorType                  # Cosine of latitudes, north to south
    coslat⁻¹::VectorType                # inverse of coslat inv.(coslat)

    # NORMALIZATION
    norm_sphere::NF                 # normalization of the l=0, m=0 mode

    # THE ACTUAL TRANSFORM MATRICES for forward = LT(FFT(input)) and backward = IFFT(ILT(input))
    forward::MatrixComplexType          # forward transform matrix
    backward::MatrixComplexType         # backward transform matrix

    # SCRATCH MEMORY
    scratch_memory::MatrixComplexType   # state is undetermined, only read after writing to it

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
        ArrayType::Type{<:AbstractArray} = DEFAULT_ARRAYTYPE,                           # Array type used for spectral coefficients (can be parametric)
        nlayers::Integer = DEFAULT_NLAYERS,                                             # number of layers in the vertical (for scratch memory size)
        LegendreShortcut::Type{<:AbstractLegendreShortcut} = LegendreShortcutLinear,    # shorten Legendre loop over order m
        architecture::AbstractArchitecture = architecture(ArrayType),                   # architecture that kernels are launched on
    )

    ArrayType_ = nonparametric_type(ArrayType)      # drop parameters of ArrayType

    # LATITUDE VECTORS (based on Gaussian, equi-angle or HEALPix latitudes)
    latd = RingGrids.get_latd(grid)                         # latitude in degrees (90˚Nto -90˚N)
    coslat = on_architecture(architecture, cosd.(latd))     # cos(lat)
    coslat⁻¹ = on_architecture(architecture, inv.(coslat))  # 1/cos(lat)

    # Create another SpectralTransform to calculate the transform matrices from (do this on the CPU)
    spectrum_cpu = on_architecture(CPU(), spectrum)
    grid_cpu = on_architecture(CPU(), grid)
    S = SpectralTransform(spectrum_cpu, grid_cpu; NF, ArrayType=Array, nlayers, LegendreShortcut, architecture=CPU())

    npoints = get_npoints(grid)
    nharmonics = LowerTriangularArrays.nonzeros(spectrum)
    forward = zeros(Complex{NF}, nharmonics, npoints)
    backward = zeros(Complex{NF}, npoints, nharmonics)

    field2D = zeros(NF, grid_cpu)
    coeffs2D = zeros(Complex{NF}, spectrum_cpu)

    progress = ProgressMeter.Progress(length(field2D) + length(coeffs2D); dt=2, desc="Precalculate matrices:")
    forward_matrix!(forward, S, field2D, coeffs2D, progress)
    backward_matrix!(backward, S, field2D, coeffs2D, progress)

    forward = on_architecture(architecture, forward)
    backward = on_architecture(architecture, backward)

    # SCRATCH MEMORY FOR FOURIER NOT YET LEGENDRE TRANSFORMED AND VICE VERSA
    scratch_memory = on_architecture(architecture, zeros(Complex{NF}, grid, nlayers).data)

    # PRECOMPUTE GRADIENT AND INTEGRATION MATRICES
    gradients = gradient_arrays(NF, spectrum)

    return MatrixSpectralTransform{
        NF,
        typeof(architecture),
        ArrayType_,
        typeof(spectrum),
        typeof(grid),
        typeof(coslat),
        typeof(forward),
        typeof(gradients),
    }(
        architecture,
        spectrum, grid, nlayers,
        coslat, coslat⁻¹,
        S.norm_sphere,
        forward,
        backward,
        scratch_memory,
        gradients,
    )
end

"""$(TYPEDSIGNATURES)
Compute the forward transform matrix `F` from grid to spectral space using the precomputed
spectral transform `S`. The matrix is computed such that multiplying it by flattened grid
data yields spectral coefficients. This function is not yet implemented."""
function forward_matrix!(F, S::AbstractSpectralTransform, field::AbstractField2D, coeffs::LowerTriangularMatrix, progress=nothing)
    for ij in eachindex(field)
        field .= 0                              # unit vector of input
        GPUArrays.@allowscalar field[ij] = 1
        transform!(coeffs, field, S)            # forward transforms of unit vectors
        F[:, ij] .= coeffs.data                 # are the columns of the transformation matrix F
        isnothing(progress) || ProgressMeter.next!(progress)
    end
end

"""$(TYPEDSIGNATURES)
Compute the backward transform matrix `B` from spectral to grid space using the precomputed
spectral transform `S`. The matrix is computed such that multiplying it by spectral
coefficients yields flattened grid data. This function is not yet implemented."""
function backward_matrix!(B, S::AbstractSpectralTransform, field::AbstractField2D, coeffs::LowerTriangularMatrix, progress=nothing)
    for ij in eachindex(coeffs)
        coeffs .= 0                             # unit vector of input
        GPUArrays.@allowscalar coeffs[ij] = 1   # TODO what's the unit vector that should be used here?
        transform!(field, coeffs, S)            # backward transforms of unit vectors
        B[:, ij] .= field.data                  # are the columns of the transformation matrix B
        isnothing(progress) || ProgressMeter.next!(progress)
    end
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
    @boundscheck ismatching(M, field, horizontal_only=true) || throw(DimensionMismatch(M, field))
    @boundscheck ismatching(M, coeffs, horizontal_only=true) || throw(DimensionMismatch(M, coeffs))
    @boundscheck size(coeffs, 2) == size(field, 2) || throw(DimensionMismatch(field.data, coeffs.data))
    LinearAlgebra.mul!(coeffs.data, M.forward, field.data)
    return coeffs
end

"""$(TYPEDSIGNATURES)
Spectral transform (spectral to grid space) from n-dimensional array `coeffs` of spherical harmonic
coefficients to an n-dimensional array `field`. Uses precomputed dense transform matrices to perform
the transformation. The spectral transform is number format-flexible but `field` and the spectral
transform `M` have to have the same number format. The spectral transform is grid-flexible as long
as `field.grid` and `M.grid` match."""
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