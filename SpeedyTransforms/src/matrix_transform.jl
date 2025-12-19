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
    latd = RingGrids.get_latd(grid)         # latitude in degrees (90˚Nto -90˚N)
    coslat = cosd.(latd)                    # cos(lat)
    coslat⁻¹ = inv.(coslat)                 # 1/cos(lat)

    # Create another SpectralTransform to calculate the transform matrices from
    S = SpectralTransform(spectrum, grid; NF, ArrayType, nlayers, LegendreShortcut, architecture)

    npoints = get_npoints(grid)
    nharmonics = LowerTriangularArrays.nonzeros(spectrum)
    forward = on_architecture(architecture, zeros(Complex{NF}, nharmonics, npoints))
    backward = on_architecture(architecture, zeros(Complex{NF}, npoints, nharmonics))

    field = on_architecture(architecture, zeros(NF, grid))
    coeffs = on_architecture(architecture, zeros(Complex{NF}, spectrum))

    forward_matrix!(forward, S, field, coeffs)
    backward_matrix!(backward, S, field, coeffs)

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

# TODO Algorithm to compute the actual transform matrices
forward_matrix!(F, S::AbstractSpectralTransform, field::AbstractField, coeffs::LowerTriangularArray) = nothing
backward_matrix!(B, S::AbstractSpectralTransform, field::AbstractField, coeffs::LowerTriangularArray) = nothing

function transform!(                        # GRID TO SPECTRAL
        coeffs::LowerTriangularArray,       # output: spectral coefficients
        field::AbstractField,               # input: gridded values
        scratch_memory,                     # explicit scratch memory (not used only in spectral to grid)
        M::MatrixSpectralTransform,         # precomputed spectral transform
    )
    LinearAlgebra.mul!(coeffs.data, M.forward, field.data)
    return coeffs
end

function transform!(                        # SPECTRAL TO GRID
        field::AbstractField,               # gridded output
        coeffs::LowerTriangularArray,       # spectral coefficients input
        scratch_memory,                     # explicit scratch memory to use
        M::MatrixSpectralTransform;         # precomputed transform
        unscale_coslat::Bool = false,       # unscale with cos(lat) on the fly?
    )
    # use first n layers of scratch memory
    scratch = view(scratch_memory, :, 1:size(coeffs, 2))
    
    # multiply into scratch which is complex typed and then take real part into field
    # imaginary part should be zero but destination is used to store intermediate results
    # explicitly convert to real also for NaN + NaN*im results
    LinearAlgebra.mul!(scratch, M.backward, coeffs.data)
    field.data .= real.(scratch)
    unscale_coslat && RingGrids._scale_lat!(field, M.coslat⁻¹)
    return field
end