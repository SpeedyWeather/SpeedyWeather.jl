module SpeedyWeatherFiniteDifferencesExt 

using SpeedyWeather 
import FiniteDifferences
import FiniteDifferences: to_vec 

# FiniteDifferences needs to be able to convert data structures to Vectors and back 
# This doesn't work out of the box with our data types, so we'll define those 
# conversions here.
function FiniteDifferences.to_vec(x::Grid) where Grid <: AbstractGridArray
    x_vec, from_vec = FiniteDifferences.to_vec(Array(x))

    function GridArray_from_vec(x_vec)
        return Grid(reshape(from_vec(x_vec), size(x)), x.nlat_half)
    end 

    return x_vec, GridArray_from_vec
end 

function FiniteDifferences.to_vec(x::LTA) where LTA <: LowerTriangularArray
    x_vec, from_vec = FiniteDifferences.to_vec(x.data)

    function LowerTriangularArray_from_vec(x_vec)
        return LowerTriangularArray(reshape(from_vec(x_vec), size(x)), x.m, x.n)
    end 

    return x_vec, LowerTriangularArray_from_vec
end

# needs an extra modification because an empty vector yields Any[] with to_vec for Particle (which isn't the case for number types)
function FiniteDifferences.to_vec(x::Vector{Particle{NF}}) where NF 
    if isempty(x) 
        return NF[], identity
    else # the else statement is the unmodified to_vec(::DenseVector)
        x_vecs_and_backs = map(to_vec, x)
        x_vecs, backs = first.(x_vecs_and_backs), last.(x_vecs_and_backs)
        function Vector_from_vec(x_vec)
            sz = cumsum(map(length, x_vecs))
            x_Vec = [backs[n](x_vec[sz[n] - length(x_vecs[n]) + 1:sz[n]]) for n in eachindex(x)]
            return oftype(x, x_Vec)
        end
        # handle empty x
        x_vec = isempty(x_vecs) ? eltype(eltype(x_vecs))[] : reduce(vcat, x_vecs)
        return x_vec, Vector_from_vec
    end 
end

end 