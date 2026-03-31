module Utils

using DocStringExtensions
using ..Architectures

import Dates

# miscellaneous utility functions (uses Dates for readable_secs)
export isincreasing, isdecreasing, clip_negatives!, underflow!
export flipsign!, nans, print_fields, readable_secs
export _jit, @maybe_jit
include("utility_functions.jl")

# kernel
export configure_kernel, launch!
export AbstractWorkOrder, SpectralWorkOrder, RingGridWorkOrder, SpectralInnerWorkOrder
export DiagonalWorkOrder, Array3DWorkOrder, LinearWorkOrder
include("kernel_launching.jl")

end
