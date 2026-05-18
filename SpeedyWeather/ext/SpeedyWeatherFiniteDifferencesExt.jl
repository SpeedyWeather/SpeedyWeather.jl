module SpeedyWeatherFiniteDifferencesExt

using SpeedyWeather
import FiniteDifferences
import FiniteDifferences: to_vec

# FiniteDifferences needs to be able to convert data structures to Vectors and back
# This doesn't work out of the box with our data types, so we'll define those
# conversions here.

# Vector{Particle} needs an extra modification because an empty vector yields Any[] with to_vec for Particle (which isn't the case for number types)
function FiniteDifferences.to_vec(x::Vector{Particle{NF}}) where {NF}
    if isempty(x)
        return NF[], identity
    else # the else statement is the unmodified to_vec(::DenseVector)
        x_vecs_and_backs = map(to_vec, x)
        x_vecs, backs = first.(x_vecs_and_backs), last.(x_vecs_and_backs)
        function Vector_from_vec(x_vec)
            sz = cumsum(map(length, x_vecs))
            x_Vec = [backs[n](x_vec[(sz[n] - length(x_vecs[n]) + 1):sz[n]]) for n in eachindex(x)]
            return oftype(x, x_Vec)
        end
        # handle empty x
        x_vec = isempty(x_vecs) ? eltype(eltype(x_vecs))[] : reduce(vcat, x_vecs)
        return x_vec, Vector_from_vec
    end
end

# TODO: We used to have an adaption here that replaced NaN's as FiniteDifferences can't deal with them, maybe we have to reintroduce it later

#=
# in the ocean and land variables we have NaNs, FiniteDifferences can't deal with those, so we replace them
function replace_NaN(x_type::T, vec) where {T <: NamedTuple}
    nan_indices = isnan.(vec)
    vec[nan_indices] .= 0
    return vec
end

# fallback, we really only want to replace the NaNs in ocean and land variables
replace_NaN(type, vec) = vec
=#

# By default FiniteDifferences doesn't include this, even though Integers can't be varied.
# there's an old GitHub issue and PR about this
function FiniteDifferences.to_vec(x::Integer)
    Integer_from_vec(v) = x
    return Bool[], Integer_from_vec
end

end
