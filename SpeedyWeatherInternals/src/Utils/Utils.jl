module Utils

using DocStringExtensions
using StyledStrings
using ..Architectures

import Dates        # Dates for readable_secs

# miscellaneous utility functions
export isincreasing, isdecreasing
export print_fields, readable_secs
export _jit, @maybe_jit
export set_scalar!, componentwise_ifelse
include("utility_functions.jl")

end
