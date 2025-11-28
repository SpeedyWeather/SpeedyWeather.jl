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
    print_fields(io, P, keys)
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
end

# Functions to get parameters and parameterization to 
# a) initialize variables 
"""$(TYPEDSIGNATURES)
Extract the model components with parameters needed for the parameterizations
as NamedTuple. These are the GPU-compatible components of the model."""
function get_model_parameters(model::PrimitiveEquation)
    names = map(field -> _get_param_name(model, field), model.model_parameters)
    values = map(field -> _get_param(model, field), model.model_parameters)
    # also include the model class for dispatching inside kernels
    return NamedTuple{names}(values)
end

"""$(TYPEDSIGNATURES)
Extract the parameterizations from the model as NamedTuple.
These are the GPU-compatible components of the model."""
function get_parameterizations(model::PrimitiveEquation)
    names = map(field -> _get_param_name(model, field), model.parameterizations)
    values = map(field -> _get_param(model, field), model.parameterizations)
    return NamedTuple{names}(values)
end

# TODO: better name? 
"""$(TYPEDSIGNATURES)
Extract the extra parameterizations from the model that are not part of the 
column-based parameterizations, but define variables such as land and ocean."""
function get_extra_parameterizations(model::PrimitiveEquation)
    names = map(field -> _get_param_name(model, field), model.extra_parameterizations)
    values = map(field -> _get_param(model, field), model.extra_parameterizations)
    return NamedTuple{names}(values)
end

get_parameterizations(model::Barotropic) = NamedTuple()
get_extra_parameterizations(model::Barotropic) = NamedTuple()

get_parameterizations(model::ShallowWater) = NamedTuple()
get_extra_parameterizations(model::ShallowWater) = NamedTuple()

# helper function to extract parameterizations from model tuples
_get_param(model, field::Symbol) = getfield(model, field)
_get_param(model, obj::Pair{Symbol}) = obj.second
_get_param(model, obj) = error("Unknown parameterization type: $(typeof(obj)), needs to be <:Symbol or <:Pair{Symbol, obj}")

_get_param_name(model, field::Symbol) = field
_get_param_name(model, obj::Pair{Symbol}) = obj.first
_get_param_name(model, obj) = error("Unknown parameterization type: $(typeof(obj)), needs to be <:Symbol or <:Pair{Symbol, obj}")

function Adapt.adapt_structure(to, model::ModelType) where ModelType <: PrimitiveEquation
    return ModelType(NamedTuple{fieldnames(model)}(
        field in merge(model.model_parameters, model.parameterizations) ? adapt_structure(to, getfield(model, field)) : nothing 
        for field in fieldnames(ModelType)
    )...)
end 