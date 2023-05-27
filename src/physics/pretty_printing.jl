# print all fields with type <: Number
function Base.show(io::IO,P::AbstractParameterization)
    print(io,"$(typeof(P)):")
    for key in propertynames(P)
        val = getfield(P,key)
        val isa Number && print(io,"\n $key::$(typeof(val)) = $val")
    end
end