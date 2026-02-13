module Utils

using DocStringExtensions
using ..Architectures

# custom dates with parametric value types — standalone replacement for Dates
# available as `using Utils.TracableDates` (drop-in for `using Dates`)
include("custom_dates.jl")
export TracableDates
import .TracableDates as Dates

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
