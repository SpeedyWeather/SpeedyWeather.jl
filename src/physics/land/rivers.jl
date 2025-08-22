abstract type AbstractRivers <: AbstractParameterization end

# no rivers
initialize!(::PrognosticVariables, ::DiagnosticVariables, ::Nothing, ::PrimitiveEquation) = nothing
timestep!(progn, diagn, rivers::Nothing, model) = nothing