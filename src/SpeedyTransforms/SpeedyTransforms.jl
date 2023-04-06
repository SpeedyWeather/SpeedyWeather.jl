module SpeedyTransforms

    import Parameters: @unpack

    # NUMERICS
    import AssociatedLegendrePolynomials as Legendre
    import AbstractFFTs
    import FFTW
    import GenericFFT
    import LinearAlgebra
    import Primes

    # SPEEDYWEATHER MODULES
    using ..LowerTriangularMatrices
    using ..RingGrids

    const DEFAULT_GRID = FullGaussianGrid

    # TRANSFORM
    export  SpectralTransform,
            gridded,
            gridded!,
            spectral,
            spectral!

    # ALIASING
    export  get_nlat_half

    # GRADIENTS
    export  curl!,
            divergence!,
            UV_from_vor!,
            UV_from_vordiv!,
            ∇²!,
            ∇⁻²!,
            ∇!

    # TRUNCATION
    export  spectral_truncation,
            spectral_truncation!,
            spectral_interpolation

    include("aliasing.jl")
    include("spectral_transform.jl")
    include("spectral_gradients.jl")
    include("spectral_truncation.jl")
end