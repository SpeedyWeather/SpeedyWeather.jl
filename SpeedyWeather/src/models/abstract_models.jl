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
finalize!(::Nothing, ::AbstractModel) = nothing

# fallback for model components: nothing to initialize
initialize!(::AbstractModelComponent, ::AbstractModel) = nothing
finalize!(::AbstractModelComponent, ::AbstractModel) = nothing

function Base.show(io::IO, P::AbstractModelComponent; values = true)
    type_str = split("$(typeof(P))", "{", limit = 2)
    type_itself = type_str[1]
    type_params = length(type_str) == 2 ? ("{" * type_str[2]) : ""
    type_params_short = length(type_params) > 30 ? first(type_params, 30) * "...}" : type_params
    println(io, styled"{warning:$type_itself}{note:$type_params_short}" * " <: $(supertype(typeof(P)))")
    keys = propertynames(P)
    print_fields(io, P, keys; values)
    return nothing
end

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

# pretty printing
function Base.show(io::IO, M::AbstractModel)
    properties = propertynames(M)
    n = length(properties)
    Msize = prettymemory(Base.summarysize(M))
    s = styled"{warning:$(model_type(M))}" * "{...} <: $(model_class(M)) " * styled"{note:($Msize)}"
    n == 0 ? print(io, s) : println(io, s)
    for (i, key) in enumerate(properties)
        val = getfield(M, key)
        s = i == n ? "└" : "├"  # choose ending └ for last property
        p = i == n ? print : println
        t = split("$(typeof(val))", "{", limit = 2)
        t1 = t[1]
        t2 = length(t) == 2 ? ("{" * t[2]) : ""
        a = "$s " * styled"{info:$key}" * "::$t1" * styled"{note:$t2}"
        a_short = textwidth(a) > 75 ? first(a, 75) * "..." : a
        p(io, a_short)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Extract the number of soil layers from the model. The fallback is 0 soil layers, i.e. no land model."""
@inline get_soil_layers(model::AbstractModel) = (hasproperty(model, :land) && !isnothing(model.land)) ? get_nlayers(model.land) : 0
@inline get_nlayers(model::AbstractModel) = model.spectral_grid.nlayers
