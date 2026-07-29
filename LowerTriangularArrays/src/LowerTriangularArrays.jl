module LowerTriangularArrays

# DOCUMENTATION
using DocStringExtensions
using StyledStrings

# NUMERICS
import LinearAlgebra: Transpose, tril!

# GPU
import Adapt: Adapt, adapt
import GPUArrays
import KernelAbstractions
import KernelAbstractions: @kernel, @index

# SPEEDYWEATHER SUBMODULES
import SpeedyWeatherInternals.Architectures: Architectures, AbstractArchitecture, on_architecture,
    array_type, ismatching, CPU, GPU, architecture, nonparametric_type
export CPU, GPU, on_architecture, architecture                # export device functions
export Architectures

import SpeedyWeatherInternals.KernelLaunching: launch!, SpectralWorkOrder
export KernelLaunching

import SpeedyWeatherInternals.ArrayDimensions: ArrayDimensions, AbstractArrayDimensions, hastime, hasvertical,
    Dimensions2D, Dimensions3D, Dimensions4D,
    DimensionsWithTime, DimensionsWithVertical, DimensionsWithTimeAndVertical
export ArrayDimensions

# CONSTANTS
const DEFAULT_NF = Float32
const DEFAULT_ARCHITECTURE = CPU

# ABSTRACT TYPES AND MAIN TYPES
export AbstractSpectrum, Spectrum, resolution, truncation

export LowerTriangularMatrix, LowerTriangularArray
export LowerTriangularArrayWithTime, LowerTriangularArrayWithVertical, LowerTriangularArrayWithTimeAndVertical
export eachharmonic, eachmatrix, eachorder, orders
export OneBased, ZeroBased
export lta_view
export zero_last_degree!

include("spectrum.jl")
include("lower_triangular_array.jl")
include("rotate_reverse.jl")
include("truncation.jl")

# Trait-based dispatch on dimensions
"""Type alias for all LowerTriangularArrays with a time dimension"""
const LowerTriangularArrayWithTime = LowerTriangularArray{T, N, ArrayType, S, Dims} where {T, N, ArrayType, S, Dims <: DimensionsWithTime}
"""Type alias for all LowerTriangularArrays with a vertical dimension"""
const LowerTriangularArrayWithVertical = LowerTriangularArray{T, N, ArrayType, S, Dims} where {T, N, ArrayType, S, Dims <: DimensionsWithVertical}
"""Type alias for all LowerTriangularArrays with both time and vertical dimensions"""
const LowerTriangularArrayWithTimeAndVertical = LowerTriangularArray{T, N, ArrayType, S, Dims} where {T, N, ArrayType, S, Dims <: DimensionsWithTimeAndVertical}

end
