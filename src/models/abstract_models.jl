export Barotropic, ShallowWater, PrimitiveEquation, PrimitiveDry, PrimitiveWet, ModelSetup

abstract type AbstractSimulation{Model} end
abstract type ModelSetup end
abstract type Barotropic <: ModelSetup end
abstract type ShallowWater <: ModelSetup end
abstract type PrimitiveEquation <: ModelSetup end
abstract type PrimitiveDry <: PrimitiveEquation end
abstract type PrimitiveWet <: PrimitiveEquation end

abstract type AbstractModelComponent end

"""$(TYPEDSIGNATURES)
Returns true if the model `M` has a prognostic variable `var_name`, false otherwise.
The default fallback is that all variables are included. """
has(M::Type{<:ModelSetup}, var_name::Symbol) = var_name in (:vor, :div, :temp, :humid, :pres)
has(M::ModelSetup, var_name) = has(typeof(M), var_name)

# strip away the parameters of the model type
model_class(::Type{<:Barotropic}) = Barotropic
model_class(::Type{<:ShallowWater}) = ShallowWater
model_class(::Type{<:PrimitiveDry}) = PrimitiveDry
model_class(::Type{<:PrimitiveWet}) = PrimitiveWet
model_class(model::ModelSetup) = model_class(typeof(model))

function Base.show(io::IO,M::ModelSetup)
    println(io,"$(typeof(M))")
    for key in propertynames(M)[1:end-1]
        val = getfield(M,key)
        println(io,"├ $key: $(typeof(val))")
    end
    print(io,"└ feedback: $(typeof(M.feedback))")
end