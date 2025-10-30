abstract type AbstractTimeStepper <: AbstractModelComponent end

abstract type AbstractVariableDims end

abstract type AbstractVariables end
abstract type AbstractPrognosticVariables <: AbstractVariables end
abstract type AbstractDiagnosticVariables <: AbstractVariables end

"""
    $TYPEDEF

Abstract type for all types of variables declared by [`variables`](@ref). 
Currently only [`PrognosticVariable`](@ref) and [`DiagnosticVariable`](@ref) 
are possible.
"""
abstract type AbstractVariable{VD<:AbstractVariableDims} end
