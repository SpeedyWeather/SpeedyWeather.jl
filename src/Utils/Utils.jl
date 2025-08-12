module Utils 

using DocStringExtensions

import Dates: Dates, DateTime, Period, Millisecond, Second, Minute, Hour, Day, Week, Month, Year

# ARRAYS
import ComponentArrays: ComponentArray, ComponentVector, labels, label2index, getaxes

# UTILITIES
import MacroTools

import ModelParameters: ModelParameters, AbstractParam

import ConstructionBase: constructorof, getproperties, setproperties

# DOMAINS
import DomainSets: Domain, RealLine, NonnegativeRealLine, PositiveRealLine, NegativeRealLine, UnitInterval

using DomainSets.IntervalSets

export Unbounded, Positive, Nonnegative, Negative

export ComponentVector

export isincreasing, isdecreasing, clip_negatives!, underflow!
export flipsign!, nans, print_fields

export DateTime, Millisecond, Second, Minute, Hour, Day, Week, Month, Year, Century, Millenium

# Domain aliases
const Unbounded = RealLine()
const Positive = PositiveRealLine()
const Nonnegative = NonnegativeRealLine()
const Negative = NegativeRealLine()

include("utility_functions.jl")

export configure_kernel, launch!
export AbstractWorkOrder, SpectralWorkOrder, RingGridWorkOrder, SpectralInnerWorkOrder
export DiagonalWorkOrder, Array3DWorkOrder, LinearWorkOrder

include("kernel_launching.jl")

# Parameter utilities
export SpeedyParam, SpeedyParams
export @parameterized, parameters, parameterof, reconstruct, stripparams, attributes, bounds, description, value

include("parameters.jl")

end 