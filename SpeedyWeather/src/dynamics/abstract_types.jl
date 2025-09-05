abstract type AbstractTimeStepper <: AbstractModelComponent end

abstract type AbstractVariables end
abstract type AbstractPrognosticVariables <: AbstractVariables end
abstract type AbstractDiagnosticVariables <: AbstractVariables end