# MODELS
abstract type ModelSetup end
abstract type Barotropic <: ModelSetup end
abstract type ShallowWater <: ModelSetup end
abstract type PrimitiveEquation <: ModelSetup end
abstract type PrimitiveDryCore <: PrimitiveEquation end
abstract type PrimitiveWetCore <: PrimitiveEquation end

# PARAMETERES, CONSTANTS, COEFFICIENTS
abstract type AbstractParameters{M} end
abstract type Coefficients end


abstract type InitialConditions end
# subtypes defined in initial_conditions.jl

abstract type AbstractOrography{NF} end
abstract type AbstractBoundaries{NF} end
# subtypes defined in boundaries.jl

abstract type AbstractColumnVariables{NF} end
abstract type AbstractConstants{NF} end
