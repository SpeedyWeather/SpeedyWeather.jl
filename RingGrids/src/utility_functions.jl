"""
$(TYPEDSIGNATURES)
Check whether elements of a vector `v` are strictly increasing."""
function isincreasing(v::AbstractVector)
    return all(diff(v) .> 0)
end

"""
$(TYPEDSIGNATURES)

Check whether elements of a vector `v` are strictly decreasing."""
function isdecreasing(v::AbstractVector)
    return all(diff(v) .< 0)
end

"""$(TYPEDSIGNATURES)
For every element vᵢ in v does a<=vi<=b hold?"""
function extrema_in(v::AbstractVector, a::Real, b::Real)
    vmin, vmax = extrema(v)
    return (vmin >= a) && (vmax <= b)
end

# MATRIX rotations
"""$(TYPEDSIGNATURES)
Rotate indices `i, j` of a square matrix of size s x s anti-clockwise by 90˚.""" 
@inline function rotate_matrix_indices_90(i::Integer, j::Integer, s::Integer)
    @boundscheck 0 < i <= s || throw(BoundsError)
    @boundscheck 0 < j <= s || throw(BoundsError)
    i_new = s+1-j   # new i from rotation
    j_new = i       # corresponding new j
    return i_new, j_new
end

"""$(TYPEDSIGNATURES)
Rotate indices `i, j` of a square matrix of size s x s by 180˚.""" 
@inline function rotate_matrix_indices_180(i::Integer, j::Integer, s::Integer)
    @boundscheck 0 < i <= s || throw(BoundsError)
    @boundscheck 0 < j <= s || throw(BoundsError)
    i_new = s+1-i   # new i from rotation
    j_new = s+1-j   # corresponding new j
    return i_new, j_new
end

"""$(TYPEDSIGNATURES)
Rotate indices `i, j` of a square matrix of size s x s anti-clockwise by 270˚.""" 
@inline function rotate_matrix_indices_270(i::Integer, j::Integer, s::Integer)
    @boundscheck 0 < i <= s || throw(BoundsError)
    @boundscheck 0 < j <= s || throw(BoundsError)
    i_new = j       # new i from rotation
    j_new = s+1-i   # corresponding new j
    return i_new, j_new
end

@inline function rotate_matrix_indices(i::Integer, j::Integer, s::Integer, r::Integer)
    r = mod(r, 4)    # map 4 to 0 rotation, 5 to 1 rotation etc.
    r == 0 && return i, j
    r == 1 && return rotate_matrix_indices_90(i, j, s)
    r == 2 && return rotate_matrix_indices_180(i, j, s)
    r == 3 && return rotate_matrix_indices_270(i, j, s)
end