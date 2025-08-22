module Utils

using DocStringExtensions

# ARCHITECTURE
using Architectures
import Architectures: device

# miscellaneous utility functions
export isincreasing, isdecreasing, clip_negatives!, underflow!
export flipsign!, nans, print_fields
include("utility_functions.jl")

# kernel
export configure_kernel, launch!
export AbstractWorkOrder, SpectralWorkOrder, RingGridWorkOrder, SpectralInnerWorkOrder
export DiagonalWorkOrder, Array3DWorkOrder, LinearWorkOrder
include("kernel_launching.jl")

end 