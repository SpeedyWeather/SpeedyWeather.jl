module Utils

using DocStringExtensions

# ARRAYS
import ComponentArrays: ComponentArray, ComponentVector, labels, label2index, getaxes

# ARCHITECTURE
include("../../Architectures/src/Architectures.jl")
using .Architectures
import .Architectures: device

# UTILITIES
import MacroTools
import ModelParameters: ModelParameters, AbstractParam
import ConstructionBase: constructorof, getproperties, setproperties

# DOMAINS
import DomainSets: Domain, RealLine, NonnegativeRealLine, PositiveRealLine, NegativeRealLine, UnitInterval
using DomainSets.IntervalSets
export Unbounded, Positive, Nonnegative, Negative
export ComponentVector

# Domain aliases
const Unbounded = RealLine()
const Positive = PositiveRealLine()
const Nonnegative = NonnegativeRealLine()
const Negative = NegativeRealLine()

# miscellaneous utility functions
export isincreasing, isdecreasing, clip_negatives!, underflow!
export flipsign!, nans, print_fields
include("utility_functions.jl")

# kernel
export configure_kernel, launch!
export AbstractWorkOrder, SpectralWorkOrder, RingGridWorkOrder, SpectralInnerWorkOrder
export DiagonalWorkOrder, Array3DWorkOrder, LinearWorkOrder
include("kernel_launching.jl")

# Parameter utilities
export SpeedyParam, SpeedyParams
export @parameterized, parameters, parameterof, reconstruct, stripparams, attributes, bounds, description, value
include("parameters.jl")

end 