# MODELS
abstract type ModelSetup end
abstract type Barotropic <: ModelSetup end
abstract type ShallowWater <: ModelSetup end
abstract type PrimitiveEquation <: ModelSetup end
abstract type PrimitiveDryCore <: PrimitiveEquation end
abstract type PrimitiveWetCore <: PrimitiveEquation end

# PARAMETERS (to be chosen by user)
abstract type AbstractParameters{M} end

# COEFFICIENTS (bundled parameters)
abstract type Coefficients end

# CONSTANTS (to be defined from parameters & coefficients,
# face the dynamical core/parameterizations and not the user)
abstract type AbstractParameterizationConstants{NF} end
abstract type AbstractDynamicsConstants{NF} end

# INITIAL CONDITIONS AND OROGRAPHY/BOUNDARIES
abstract type InitialConditions end         # subtypes defined in initial_conditions.jl
abstract type AbstractOrography end         # subtypes defined in boundaries.jl
abstract type AbstractBoundaries{NF} end

# ATMOSPHERIC COLUMN FOR PARAMETERIZATIONS
abstract type AbstractColumnVariables{NF} end

# PARAMETERIZATIONS
abstract type BoundaryLayer end

# INPUT/OUTPUT
abstract type AbstractFeedback end
abstract type AbstractOutput end

# IMPLICIT PRECOMPUTED TERMS
abstract type AbstractImplicit{NF} end
