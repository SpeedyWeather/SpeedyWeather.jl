export Barotropic, ShallowWater, PrimitiveEquation, PrimitiveDry, PrimitiveWet, AbstractModel

abstract type AbstractSimulation{Model} end
abstract type AbstractModel end
abstract type Barotropic <: AbstractModel end
abstract type ShallowWater <: AbstractModel end
abstract type PrimitiveEquation <: AbstractModel end
abstract type PrimitiveDry <: PrimitiveEquation end
abstract type PrimitiveWet <: PrimitiveEquation end

abstract type AbstractModelComponent end

# any model component set to nothing needs no initialization or finalize!
initialize!(::Nothing, ::AbstractModel) = nothing

# allow for components that initialize variables to be nothing
initialize!(p, d, ::Nothing, ::AbstractModel) = nothing
# or sub parts of the prognostic variables; this makes land=nothing, ocean=nothing possible
initialize!(psub, p, d, ::Nothing, ::AbstractModel) = nothing
timestep!(p, d, ::Nothing, ::AbstractModel) = nothing
finalize!(::Nothing, ::AbstractModel) = nothing

# fallback for model components: nothing to initialize
initialize!(::AbstractModelComponent, ::AbstractModel) = nothing
finalize!(::AbstractModelComponent, ::AbstractModel) = nothing

# print all fields with type <: Number
function Base.show(io::IO, P::AbstractModelComponent)
    println(io, "$(typeof(P)) <: $(supertype(typeof(P)))")
    keys = propertynames(P)
    return print_fields(io, P, keys)
end

"""$(TYPEDSIGNATURES)
Returns true if the model `M` has a prognostic variable `var_name`, false otherwise.
The default fallback is that all variables are included. """
has(Model::Type{<:AbstractModel}, var_name::Symbol) = var_name in prognostic_variables(Model)
has(model::AbstractModel, var_name) = has(typeof(model), var_name)
prognostic_variables(model::AbstractModel) = prognostic_variables(typeof(model))

# model class is the abstract supertype
model_class(::Type{<:Barotropic}) = Barotropic
model_class(::Type{<:ShallowWater}) = ShallowWater
model_class(::Type{<:PrimitiveDry}) = PrimitiveDry
model_class(::Type{<:PrimitiveWet}) = PrimitiveWet
model_class(model::AbstractModel) = model_class(typeof(model))

# model type is the parameter-free type of a model
model_type(::Type{<:Barotropic}) = BarotropicModel
model_type(::Type{<:ShallowWater}) = ShallowWaterModel
model_type(::Type{<:PrimitiveDry}) = PrimitiveDryModel
model_type(::Type{<:PrimitiveWet}) = PrimitiveWetModel
model_type(model::AbstractModel) = model_type(typeof(model))

initialize!(model::AbstractModel, ps::Union{ComponentVector, SpeedyParams}; kwargs...) =
    initialize!(reconstruct(model, ps); kwargs...)

function Base.show(io::IO, M::AbstractModel)
    properties = propertynames(M)
    n = length(properties)
    s = "$(model_type(M)) <: $(model_class(M))"
    n == 0 ? print(io, s) : println(io, s)
    for (i, key) in enumerate(properties)
        val = getfield(M, key)
        s = i == n ? "└" : "├"  # choose ending └ for last property
        p = i == n ? print : println
        a = "$s $key: $(typeof(val))"
        a = textwidth(a) > 100 ? string(a[1:97], "...") : a  # truncate long strings
        p(io, a)
    end
    return
end

# Functions to get parameters and parameterization to
# a) initialize variables
"""$(TYPEDSIGNATURES)
Extract the model components with parameters needed for the parameterizations
as NamedTuple. These are the GPU-compatible components of the model."""
@inline function get_model_parameters(model::PrimitiveEquation)
    return NamedTuple{model.model_parameters}(map(field -> getfield(model, field), model.model_parameters))
end

"""$(TYPEDSIGNATURES)
Extract the parameterizations from the model as NamedTuple.
These are the GPU-compatible components of the model."""
@generated function get_parameterizations(model::ModelType) where {ModelType <: PrimitiveEquation}
    # Extract parameterization symbols from the type
    params_type = fieldtype(ModelType, :params)
    param_names = params_type.parameters[1]  # Extract tuple from Val{tuple}
    
    # Generate literal field accesses for type stability
    return :(NamedTuple{$param_names}(tuple($([:(model.$name) for name in param_names]...))))
end

# TODO: better name?
"""$(TYPEDSIGNATURES)
Extract the extra parameterizations from the model that are not part of the 
column-based parameterizations, but define variables such as land and ocean."""
@inline function get_extra_parameterizations(model::PrimitiveEquation)
    return NamedTuple{model.extra_parameterizations}(map(field -> getfield(model, field), model.extra_parameterizations))
end

@inline get_parameterizations(model::Barotropic) = NamedTuple()
@inline get_extra_parameterizations(model::Barotropic) = NamedTuple()

@inline get_parameterizations(model::ShallowWater) = NamedTuple()
@inline get_extra_parameterizations(model::ShallowWater) = NamedTuple()
