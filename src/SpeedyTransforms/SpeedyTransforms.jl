module SpeedyTransforms

using DocStringExtensions, Printf

# NUMERICS
import AssociatedLegendrePolynomials as Legendre
import AbstractFFTs
import FFTW
import GenericFFT
import LinearAlgebra
import Primes
import Adapt: adapt

# SPEEDYWEATHER MODULES
using ..LowerTriangularMatrices
using ..RingGrids

# TRANSFORM
export  SpectralTransform,
        transform!,
        transform

# ALIASING
export  get_nlat_half

# GRADIENTS
export  curl,
        divergence,
        curl!,
        divergence!,
        UV_from_vor!,
        UV_from_vordiv!,
        ∇²!, ∇⁻²!, ∇!,
        ∇², ∇⁻², ∇

# TRUNCATION
export  spectral_truncation,
        spectral_truncation!,
        spectral_interpolation

# ANALYSIS
export  power_spectrum

include("aliasing.jl")
include("legendrepolarray.jl")
include("legendre_shortcuts.jl")
include("spectral_transform.jl")
include("fourier.jl")
include("legendre.jl")
include("spectral_gradients.jl")
include("spectral_truncation.jl")
include("spectrum.jl")
include("show.jl")

end