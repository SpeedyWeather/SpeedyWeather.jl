export Barotropic, ShallowWater, PrimitiveEquation, PrimitiveDry, PrimitiveWet, AbstractModel

abstract type AbstractSimulation{Model} end
abstract type AbstractModel end
abstract type Barotropic <: AbstractModel end
abstract type ShallowWater <: AbstractModel end
abstract type PrimitiveEquation <: AbstractModel end
abstract type PrimitiveDry <: PrimitiveEquation end
abstract type PrimitiveWet <: PrimitiveEquation end

abstract type AbstractModelComponent end

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
# TODO what happens if we have several concrete types under each abstract type?
model_type(::Type{<:Barotropic}) = BarotropicModel
model_type(::Type{<:ShallowWater}) = ShallowWaterModel
model_type(::Type{<:PrimitiveDry}) = PrimitiveDryModel
model_type(::Type{<:PrimitiveWet}) = PrimitiveWetModel
model_type(model::AbstractModel) = model_type(typeof(model))

initialize!(model::AbstractModel, ps::Union{ComponentVector,SpeedyParams}; kwargs...) = initialize!(reconstruct(model, ps); kwargs...)

function Base.show(io::IO, M::AbstractModel)
    println(io, "$(model_type(M)) <: $(model_class(M))")
    properties = propertynames(M)
    n = length(properties)
    for (i, key) in enumerate(properties)
        val = getfield(M, key)
        s = i == n ? "└" : "├"  # choose ending └ for last property
        p = i == n ? print : println
        a = "$s $key: $(typeof(val))"
        a = textwidth(a) > 100 ? string(a[1:97], "...") : a  # truncate long strings
        p(io, a)
    end
end