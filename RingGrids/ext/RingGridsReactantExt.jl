module RingGridsReactantExt

using RingGrids
using Reactant
using KernelAbstractions
using SpeedyWeatherInternals.Architectures: ReactantDevice, nonparametric_type
using SpeedyWeatherInternals.Utils: _jit
using SpeedyWeatherInternals.KernelLaunching: launch!, LinearWorkOrder

const AnyReactantArray = Union{Reactant.AnyConcretePJRTArray, Reactant.AnyConcreteIFRTArray, Reactant.AnyTracedRArray}
const ReactantGridGeometry = RingGrids.GridGeometry{<:AbstractGrid{<:ReactantDevice}}

# BROADCASTING
# Fields with Reactant arrays (concrete or traced) get their own broadcast style, analogous to
# FieldGPUStyle: out-of-place broadcasts allocate a Field (reusing the grid of the first field),
# in-place broadcasts unwrap all fields and broadcast into the data via Reactant, compiled
# outside of and traced inside of Reactant's compilation.
struct FieldReactantStyle{N, Grid} <: Base.Broadcast.AbstractArrayStyle{N} end

# ::Val{M} for broadcasting with M-dimensional arrays, choose the larger dimensionality
FieldReactantStyle{N, Grid}(::Val{M}) where {N, Grid, M} = FieldReactantStyle{M, Grid}()

# `Type{<:...}` to match the concrete Field type including all its type parameters
Base.Broadcast.BroadcastStyle(::Type{<:AbstractField{T, N, ArrayType, Grid}}) where {T, N, ArrayType <: AnyReactantArray, Grid} =
    FieldReactantStyle{N, nonparametric_type(Grid)}()

# combine with Reactant's own concrete and traced array styles (choosing the field style)
const ReactantArrayStyle = Union{
    Base.Broadcast.ArrayStyle{ConcretePJRTArray}, Base.Broadcast.ArrayStyle{ConcreteIFRTArray},
    Reactant.TracedRArrayOverrides.AbstractReactantArrayStyle,
}
Base.Broadcast.BroadcastStyle(s::FieldReactantStyle, ::ReactantArrayStyle) = s
Base.Broadcast.BroadcastStyle(::ReactantArrayStyle, s::FieldReactantStyle) = s
Base.Broadcast.BroadcastStyle(::FieldReactantStyle{M, Grid}, ::FieldReactantStyle{N, Grid}) where {M, N, Grid} =
    FieldReactantStyle{max(M, N), Grid}()

function Base.similar(bc::Base.Broadcast.Broadcasted{FieldReactantStyle{N, Grid}}, ::Type{T}) where {N, Grid, T}
    field = RingGrids.find_field(bc)
    # parent of broadcasted arrays is used because we don't want e.g. a view as a result
    data = similar(parent(field.data), Reactant.unwrapped_eltype(T), axes(bc))
    return nonparametric_type(typeof(field))(data, field.grid, field.dims)
end

# replace fields by their data
unwrap_fields(field::AbstractField) = field.data
unwrap_fields(bc::Base.Broadcast.Broadcasted) = Base.Broadcast.broadcasted(bc.f, map(unwrap_fields, bc.args)...)
unwrap_fields(x) = x

broadcast_into!(dest, f, args...) = (dest .= f.(args...); nothing)

function Base.copyto!(dest::AbstractArray, bc::Base.Broadcast.Broadcasted{<:FieldReactantStyle})
    dest_data = dest isa AbstractField ? dest.data : dest
    bc_data = Base.Broadcast.flatten(unwrap_fields(bc))
    _jit(ReactantDevice(), broadcast_into!, dest_data, bc_data.f, bc_data.args...)
    return dest
end

# Bounds checks on index arrays: values of traced arrays are unknown during tracing so the check
# is skipped, concrete arrays are checked on the host (`extrema` doesn't work on Reactant arrays,
# see https://github.com/EnzymeAD/Reactant.jl/issues/2387)
RingGrids.extrema_in(v::AbstractVector{<:Reactant.TracedRNumber}, a::Real, b::Real) = true
RingGrids.extrema_in(v::Union{ConcretePJRTArray{T, 1}, ConcreteIFRTArray{T, 1}}, a::Real, b::Real) where {T} =
    RingGrids.extrema_in(Array(v), a, b)

# Reshaped Reactant arrays can't be passed to kernels (`reshape(::CuTracedArray, dims)` is broken
# in Reactant, its `_derived_array` is undefined), so only "flatten" matrices which don't need a reshape.
# Vectors and views are interpolated with the single-layer kernel, higher-dimensional arrays per layer
# TODO: remove once reshape(::CuTracedArray, dims) works in Reactant
const ReactantArray = Union{Reactant.TracedRArray, ConcretePJRTArray, ConcreteIFRTArray}
RingGrids.is_flattenable(x::ReactantArray) = ndims(x) == 2
RingGrids.is_flattenable(x::SubArray{<:Any, <:Any, <:ReactantArray}) = false

# Pole averages of concrete (not traced) arrays: `mean` falls back to scalar indexing (`first`), which
# is disallowed on GPU, so average on the host instead. Only used eagerly outside of `@compile`
# (e.g. interpolating boundary conditions in `initialize!`), traced arrays use Reactant's `mean`
RingGrids.average_on_poles(A::Union{ConcretePJRTArray{T, 1}, ConcreteIFRTArray{T, 1}}, rings) where {T} =
    RingGrids.average_on_poles(Array(A), rings)

# Same for `extrema` of fields with concrete data (also views, e.g. into fused parents), used eagerly
# by the progress feedback between compiled time steps
Base.extrema(F::RingGrids.AbstractField{T, N, <:Union{Reactant.AnyConcretePJRTArray, Reactant.AnyConcreteIFRTArray}}) where {T, N} =
    extrema(Array(F.data))

# Force `unsafe=true` in `update_locator!` for Reactant: the safety checks in `find_rings!`
# call `extrema(θs)` and assertions on `latd` that don't trace through Reactant.
function RingGrids.update_locator!(
        locator::RingGrids.AbstractLocator,
        geometry::ReactantGridGeometry,
        λs::AbstractVector,
        θs::AbstractVector;
        unsafe::Bool = true,   # ignored, always unsafe with Reactant
    )
    (; latd) = geometry
    (; js, Δys) = locator
    RingGrids.find_rings_unsafe!(js, Δys, θs, latd, geometry.grid.architecture)
    return RingGrids.find_grid_indices!(locator, geometry, λs, geometry.grid.architecture)
end

# Branchless kernel for finding ring indices: counts how many latd entries are >= θ
# to obtain the ring index j such that latd[j] >= θ > latd[j+1]. This avoids the
# data-dependent control flow (while / break / elseif) of the binary search in the main-package
# `find_rings_kernel!`, which is faster on CPU/GPU but fails to compile within Reactant's `@jit`
# (Reactant 0.2.292: "Deactivating wrong context"), hence this Reactant-only variant.
@kernel inbounds = true function find_rings_kernel_branchless!(
        js,        # Out: ring indices j
        Δys,       # Out: distance fractions to ring further south
        θs,        # latitudes to interpolate onto
        latd,      # latitudes of the rings on the original grid (strictly decreasing)
        n::Int     # length(latd), passed explicitly to avoid dynamic dispatch
    )
    k = @index(Global, Linear)
    θ = θs[k]

    j = 0
    for i in 1:n
        j += ifelse(latd[i] >= θ, 1, 0)
    end

    # clamp to valid range [1, n-1] to avoid out-of-bounds on latd[j+1]
    j = max(j, 1)
    j = min(j, n - 1)

    js[k] = j - 1                                       # 0 = north pole, nlat = south pole
    Δys[k] = (latd[j] - θ) / (latd[j] - latd[j + 1])
end

# Dispatch find_rings_unsafe! for ReactantDevice to use the branchless kernel
function RingGrids.find_rings_unsafe!(
        js::AbstractArray,
        Δys::AbstractArray,
        θs::AbstractArray,
        latd::AbstractArray,
        architecture::ReactantDevice,
    )
    # TODO: include boundschecks again when issue with Reactant resolved
    return launch!(
        architecture,
        LinearWorkOrder,
        size(js),
        find_rings_kernel_branchless!,
        js,
        Δys,
        θs,
        latd,
        length(latd),
    )
end

end
