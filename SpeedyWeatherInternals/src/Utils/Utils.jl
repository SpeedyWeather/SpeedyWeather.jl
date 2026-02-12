module Utils

using DocStringExtensions
using Dates
using ..Architectures

# miscellaneous utility functions
export isincreasing, isdecreasing, clip_negatives!, underflow!
export flipsign!, nans, print_fields, readable_secs
export _jit, @maybe_jit
include("utility_functions.jl")

# kernel
export configure_kernel, launch!
export AbstractWorkOrder, SpectralWorkOrder, RingGridWorkOrder, SpectralInnerWorkOrder
export DiagonalWorkOrder, Array3DWorkOrder, LinearWorkOrder
include("kernel_launching.jl")

# custom dates with parametric value types for Reactant compatibility
export SpeedyPeriod, SpeedyFixedPeriod
export SpeedyMillisecond, SpeedySecond, SpeedyMinute, SpeedyHour, SpeedyDay, SpeedyWeek, SpeedyMonth, SpeedyYear
export SpeedyUTInstant, SpeedyDateTime
export secondofday
include("custom_dates.jl")

end
