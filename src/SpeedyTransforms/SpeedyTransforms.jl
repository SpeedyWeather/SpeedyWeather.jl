module SpeedyTransforms

    import Parameters: @unpack

    # NUMERICS
    import AssociatedLegendrePolynomials as Legendre
    import AbstractFFTs
    import FFTW
    import GenericFFT
    import LinearAlgebra

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

    include("spectral_transform.jl")
    include("spectral_gradients.jl")
    include("spectral_truncation.jl")
end