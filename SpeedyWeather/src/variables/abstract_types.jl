abstract type AbstractVariableDim end                           # Dimensions of variables
abstract type AbstractVariable{D <: AbstractVariableDim} end    # Variable type with dimensions D
abstract type AbstractVariables end                             # Container for all variable arrays in the model