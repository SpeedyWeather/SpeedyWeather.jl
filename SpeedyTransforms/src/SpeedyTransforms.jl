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
import KernelAbstractions: @kernel, @index, @Const, synchronize
import Atomix

# SPEEDYWEATHER MODULES
include("../../Architectures/src/Architectures.jl")
using .Architectures

include("../../Utils/src/Utils.jl")
using .Utils

include("../../RingGrids/src/RingGrids.jl")
using .RingGrids

include("../../LowerTriangularArrays/src/LowerTriangularArrays.jl")
using .LowerTriangularArrays

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
include("legendre_shortcuts.jl")
include("scratch_memory.jl")
include("spectral_transform.jl")
include("fourier.jl")
include("legendre.jl")
include("legendre_ka.jl")
include("spectral_gradients.jl")
include("spectral_truncation.jl")
include("power_spectrum.jl")
include("show.jl")

end