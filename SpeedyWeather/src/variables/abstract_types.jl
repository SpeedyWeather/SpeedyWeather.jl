abstract type AbstractTimeStepper <: AbstractModelComponent end

abstract type AbstractVariableDim end                           # Dimensions of variables
abstract type AbstractVariable{D <: AbstractVariableDim} end    # Variable type with dimensions D
abstract type AbstractVariables end                             # Container for all variable arrays in the model

# TODO remove. Define here to pass parsing while testing
struct PrognosticVariables end
struct DiagnosticVariables end
struct Tendencies end
