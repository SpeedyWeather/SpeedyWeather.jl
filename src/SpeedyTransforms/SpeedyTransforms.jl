module SpeedyTransforms

    import Parameters: @unpack

    # NUMERICS
    import FastGaussQuadrature
    import AssociatedLegendrePolynomials as Legendre
    import AbstractFFTs
    import FFTW
    import GenericFFT
    import Primes
    import LinearAlgebra

    # SPEEDYWEATHER MODULES
    using SpeedyWeather.LowerTriangularMatrices
    using SpeedyWeather.RingGrids
    import SpeedyWeather: AbstractParameters    # TODO remove this dependency

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

    include("utility_functions.jl")

    include("spectral_transform.jl")
    include("spectral_gradients.jl")
    include("spectral_truncation.jl")
end