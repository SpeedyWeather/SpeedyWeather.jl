abstract type AbstractParameterization <: AbstractModelComponent end

# print all fields with type <: Number
function Base.show(io::IO,P::AbstractModelComponent)
    println(io,"$(typeof(P)) <: $(supertype(typeof(P)))")
    keys = propertynames(P)
    print_fields(io,P,keys)
end