abstract type AbstractRivers <: AbstractParameterization end

export NoRivers
struct NoRivers <: AbstractRivers end
NoRivers(SG::SpectralGrid) = NoRivers()
initialize!(rivers::NoRivers, model::PrimitiveEquation) = nothing
initialize!(::PrognosticVariables, ::DiagnosticVariables, ::NoRivers, ::PrimitiveEquation) = nothing
timestep!(progn::PrognosticVariables, diagn::DiagnosticVariables, rivers::NoRivers, model::PrimitiveEquation) = nothing