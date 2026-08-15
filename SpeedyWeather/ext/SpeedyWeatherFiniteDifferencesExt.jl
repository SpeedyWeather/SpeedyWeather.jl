module SpeedyWeatherFiniteDifferencesExt

using SpeedyWeather
import FiniteDifferences
import FiniteDifferences: to_vec

# FiniteDifferences needs to be able to convert data structures to Vectors and back
# This doesn't work out of the box with our data types, so we'll define those
# conversions here.

# --- Variables (fused) ---------------------------------------------------------------------------
# The fused `Variables` store some prognostic/grid leaves as VIEWS into a shared parent buffer that
# lives under `vars.fused.*`. FiniteDifferences' generic struct `to_vec` mishandles this: it (a)
# vectorizes both the view and its parent (double-counting the same data) and (b) cannot reconstruct
# the `SubArray`-typed view field (`convert(::SubArray, ::Array)` error inside `NamedTuple_from_vec`).
# We mirror `Base.copy!(::Variables, ::Variables)`: vectorize only the non-view floating/complex
# array leaves (the fuse parents are non-view and hold the shared data), and reconstruct in place on
# a `deepcopy` template. `deepcopy` preserves the view→parent aliasing, so filling the parents
# updates the aliasing views automatically; the views themselves are skipped.
function FiniteDifferences.to_vec(vars::SpeedyWeather.Variables)
    template = deepcopy(vars)
    leaves = _fd_diff_leaves(vars)
    vecs_backs = [FiniteDifferences.to_vec(leaf) for leaf in leaves]
    lengths = Int[length(first(vb)) for vb in vecs_backs]
    x_vec = isempty(vecs_backs) ? Float64[] : reduce(vcat, first.(vecs_backs))

    function Variables_from_vec(v)
        out = deepcopy(template)
        outleaves = _fd_diff_leaves(out)
        offset = 0
        for (i, leaf) in enumerate(outleaves)
            n = lengths[i]
            copyto!(leaf, vecs_backs[i][2](v[(offset + 1):(offset + n)]))
            offset += n
        end
        return out
    end
    return x_vec, Variables_from_vec
end

# Collect the differentiable (floating/complex) non-view array backing of a Variables tree, in the
# same deterministic order for the original and any deepcopy. Skips view leaves (their data lives in
# a fuse parent, which is itself collected) — see `SpeedyWeather.is_view_entry` and `copy!`.
function _fd_diff_leaves(vars::SpeedyWeather.Variables)
    acc = AbstractArray[]
    for group in SpeedyWeather.ALL_VARIABLE_GROUPS
        _fd_push_leaves!(acc, getfield(vars, group))
    end
    return acc
end
_fd_push_leaves!(acc, x::NamedTuple) = foreach(k -> _fd_push_leaves!(acc, getfield(x, k)), keys(x))
_fd_push_leaves!(acc, x::SpeedyWeather.FusedParent) = _fd_push_leaves!(acc, x.parent)   # the fuse buffer
_fd_push_leaves!(acc, x::LowerTriangularArray) = SpeedyWeather.is_view_entry(x) || push!(acc, x.data)
_fd_push_leaves!(acc, x::SpeedyWeather.RingGrids.AbstractField) = SpeedyWeather.is_view_entry(x) || push!(acc, x.data)
_fd_push_leaves!(acc, ::SubArray) = nothing
_fd_push_leaves!(acc, x::AbstractArray) = eltype(x) <: Union{AbstractFloat, Complex} && push!(acc, x)
_fd_push_leaves!(acc, x) = nothing

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
