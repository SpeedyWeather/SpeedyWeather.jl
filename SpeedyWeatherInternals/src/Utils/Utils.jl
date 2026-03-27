module Utils

using DocStringExtensions
using StyledStrings
using ..Architectures

import Dates

# miscellaneous utility functions (uses Dates for readable_secs)
export isincreasing, isdecreasing
export print_fields, readable_secs
export _jit, @maybe_jit
include("utility_functions.jl")

# kernel
export configure_kernel, launch!
export AbstractWorkOrder, SpectralWorkOrder, RingGridWorkOrder, SpectralInnerWorkOrder
export DiagonalWorkOrder, Array3DWorkOrder, LinearWorkOrder
include("kernel_launching.jl")

end
