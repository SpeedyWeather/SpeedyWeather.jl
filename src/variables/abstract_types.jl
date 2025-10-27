abstract type AbstractTimeStepper <: AbstractModelComponent end

abstract type AbstractVariableDims end
abstract type AbstractVariable{VD<:AbstractVariableDims} end

abstract type AbstractVariables end
abstract type AbstractPrognosticVariables <: AbstractVariables end
abstract type AbstractDiagnosticVariables <: AbstractVariables end