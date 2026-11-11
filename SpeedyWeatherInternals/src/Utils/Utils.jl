module Utils

using DocStringExtensions

# DocStringExtensions v0.10 can show default values of (keyword) arguments in signatures,
# v0.9 can't (its `TypedMethodSignatures` has no `defaults` keyword) so fall back to the plain one
const TYPEDSIGNATURES = if hasmethod(DocStringExtensions.TypedMethodSignatures, Tuple{}, (:defaults,))
    DocStringExtensions.TypedMethodSignatures(defaults = true)
else
    DocStringExtensions.TYPEDSIGNATURES
end

using StyledStrings
using ..Architectures

import Dates        # Dates for readable_secs

# miscellaneous utility functions
export isincreasing, isdecreasing
export print_fields, readable_secs
export _jit, @maybe_jit
export set_scalar!
include("utility_functions.jl")

end
