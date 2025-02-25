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

# Vector{Particle} needs an extra modification because an empty vector yields Any[] with to_vec for Particle (which isn't the case for number types)
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

# A version of the generic fallback from FiniteDifferences that excludes some of the fields 
# that we don't want to be varied for our big data structures 
# also replaces NaNs that are expected in land and ocean variables
function FiniteDifferences.to_vec(x::T) where {T <: Union{PrognosticVariables, PrognosticVariablesOcean, PrognosticVariablesLand, DiagnosticVariables, Tendencies, GridVariables, DynamicsVariables, PhysicsVariables, ParticleVariables}}

    excluded_fields_pre, included_fields, excluded_fields_post = determine_included_fields(T)

    val_vecs_and_backs = map(name -> to_vec(getfield(x, name)), included_fields)
    vals = first.(val_vecs_and_backs)
    backs = last.(val_vecs_and_backs)

    vals_excluded_pre = map(name -> getfield(x, name), excluded_fields_pre)
    vals_excluded_post = map(name -> getfield(x, name), excluded_fields_post)

    v, vals_from_vec = to_vec(vals)
    v = replace_NaN(x, v)

    function structtype_from_vec(v::Vector{<:Real})
        val_vecs = vals_from_vec(v)
        values = map((b, v) -> b(v), backs, val_vecs)
        
        T(vals_excluded_pre..., values..., vals_excluded_post...)
    end
    return v, structtype_from_vec
end

function determine_included_fields(T::Type)
    names = fieldnames(T)

    included_field_types = Union{SpeedyWeather.AbstractDiagnosticVariables, 
    SpeedyWeather.AbstractPrognosticVariables, SpeedyWeather.ColumnVariables,
    NTuple, Dict{Symbol, <:Tuple}, Dict{Symbol, <:AbstractArray}, AbstractArray}

    excluded_fields_pre = []
    included_fields = []
    excluded_fields_post = []

    for name in names 
        if fieldtype(T, name) <: included_field_types
            push!(included_fields, name)
        else 
            if isempty(included_fields)
                push!(excluded_fields_pre, name)
            else 
                push!(excluded_fields_post, name)
            end 
        end 
    end 

    return excluded_fields_pre, included_fields, excluded_fields_post
end 

# in the ocean and land variables we have NaNs, FiniteDifferences can't deal with those, so we replace them
function replace_NaN(x_type::T, vec) where {T <: Union{PrognosticVariablesOcean, PrognosticVariablesLand, PhysicsVariables}}
    nan_indices = isnan.(vec)
    vec[nan_indices] .= 0 
    return vec
end 

# fallback, we really only want to replace the NaNs in ocean and land variables 
replace_NaN(type, vec) = vec

# By default FiniteDifferences doesn't include this, even though Integers can't be varied. 
# there's an old GitHub issue and PR about this 
function FiniteDifferences.to_vec(x::Integer)
    Integer_from_vec(v) = x
    return Bool[], Integer_from_vec
end

end 