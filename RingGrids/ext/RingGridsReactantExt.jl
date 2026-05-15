module RingGridsReactantExt

using RingGrids
using Reactant
using KernelAbstractions
using SpeedyWeatherInternals.Architectures: ReactantDevice
using SpeedyWeatherInternals.KernelLaunching: launch!, LinearWorkOrder

const AnyReactantArray = Union{Reactant.AnyConcretePJRTArray, Reactant.AnyConcreteIFRTArray, Reactant.AnyTracedRArray}

# Route Field broadcasting through the underlying Reactant array's style so that
# broadcasts (e.g. `field .*= scale`) get JIT-compiled instead of iterating scalar.
# Must dispatch on `Type{F} where F <: AbstractField{...}` (not `Type{AbstractField{...}}`)
# to match the concrete `Field{...}` type Julia sees at a broadcast site.
Base.Broadcast.BroadcastStyle(::Type{F}) where {F <: AbstractField{T, N, ArrayType, Grid}} where {T, N, ArrayType <: AnyReactantArray, Grid} = Base.Broadcast.BroadcastStyle(ArrayType)

# Unwrap Fields with Reactant-backed data to their underlying `.data` so that Reactant's
# broadcast JIT compiler sees plain ConcretePJRTArrays it knows how to trace (instead of
# falling back to scalar `field[i]` via Field's AbstractArray interface).
Base.Broadcast.broadcastable(f::AbstractField{T, N, ArrayType}) where {T, N, ArrayType <: AnyReactantArray} = f.data

# Forward `copyto!(field, bc)` to `copyto!(field.data, bc)` so Reactant's broadcast copyto!
# can write directly into the ConcretePJRTArray instead of materializing into a Vector and
# scalar-iterating into the Field (which Reactant's AbstractArray fallback would do).
for AT in (:ConcretePJRTArray, :ConcreteIFRTArray)
    @eval function Base.copyto!(
            dest::AbstractField{T, N, ArrayType},
            bc::Base.Broadcast.Broadcasted{Base.Broadcast.ArrayStyle{$AT}},
        ) where {T, N, ArrayType <: AnyReactantArray}
        Base.copyto!(dest.data, bc)
        return dest
    end
end

# Branchless kernel for finding ring indices, compatible with Reactant tracing.
# It's slower than the regular kernel on both CPU and GPU, that's why it's only used here
# Counts how many latd entries are >= θ to find the ring index j
# such that latd[j] >= θ > latd[j+1]. This avoids data-dependent control flow
# (while, break, if/elseif) which is incompatible with Reactant tracing.
@kernel inbounds = true function find_rings_kernel_branchless!(
        js,                # Out: ring indices j
        Δys,               # Out: distance fractions to ring further south
        θs,        # latitudes to interpolate onto
        latd,      # latitudes of the rings on the original grid
        n::Int     # length(latd), passed explicitly to avoid dynamic dispatch
    )
    k = @index(Global, Linear)
    θ = θs[k]

    # latd is strictly decreasing, count how many entries satisfy latd[i] >= θ
    # this gives j such that latd[j] >= θ > latd[j+1]
    j = 0
    for i in 1:n
        j += ifelse(latd[i] >= θ, 1, 0)
    end

    # clamp to valid range [1, n-1] to avoid out-of-bounds on latd[j+1]
    j = max(j, 1)
    j = min(j, n - 1)

    # Store results
    js[k] = j - 1  # 0 = north pole, nlat = south pole
    Δys[k] = (latd[j] - θ) / (latd[j] - latd[j + 1])
end

# With Reactant, always use unsafe=true to skip safety checks that use
# data-dependent control flow (assertions, extrema) incompatible with tracing
const ReactantGridGeometry = RingGrids.GridGeometry{<:AbstractGrid{<:ReactantDevice}}

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

#TODO: Most of these below are temporrary and will be consolidated with the regular code again, once everythign else works
# Cast the element type of a TracedRArray via stablehlo.convert (device-local, no CPU hop)
_cast_eltype(::Type{T}, x) where {T} = Reactant.Ops.convert(Reactant.TracedRArray{T, ndims(x)}, x)

# GPU-friendly longitude index lookup. `λ` assumed in [0, 360), `λ₀` in [0, 360).
# Uses `unsafe_trunc` (not `floor(Int, x)`) and arithmetic wrap-around to avoid ops
# that emit Julia type-object references in PTX (`jl_float32_type`).
@inline function _find_lon_indices_reactant(λ, λ₀, nlon)
    NF = typeof(λ)
    Δλ = NF(360) / NF(nlon)
    # shift so (λ - λ₀)/Δλ is non-negative, then trunc == floor
    ix = (λ - λ₀ + NF(360)) / Δλ
    i = unsafe_trunc(Int, ix)
    Δ = ix - NF(i)
    # periodic wrap via subtraction: i ∈ [0, 2*nlon) after the +360 shift above
    i = ifelse(i >= nlon, i - nlon, i)
    i_a = i + 1
    i_b = i + 2
    i_b = ifelse(i_b > nlon, i_b - nlon, i_b)
    return i_a, i_b, Δ
end

# Reactant-specific kernel: `ring_starts[j]` is the first grid-point index of ring j,
# so `rings[j][i] = ring_starts[j] + i - 1`. This avoids the large static-tuple
# dynamic indexing that breaks PTX codegen.
@kernel inbounds = true function find_grid_indices_kernel_reactant!(
        js,
        ij_as, ij_bs,
        ij_cs, ij_ds,
        Δabs, Δcds,
        λs,
        lon_offsets,
        nlons,
        nlat,
        ring_starts,
    )
    k = @index(Global, Linear)

    j = js[k]
    λ = λs[k]

    # NORTHERN POINTS a, b
    if j == 0
        ij_as[k] = 0
        ij_bs[k] = 0
    else
        i_a, i_b, Δ = _find_lon_indices_reactant(λ, lon_offsets[j], nlons[j])
        ij_as[k] = ring_starts[j] + i_a - 1
        ij_bs[k] = ring_starts[j] + i_b - 1
        Δabs[k] = Δ
    end

    # SOUTHERN POINTS c, d
    if j == nlat
        ij_cs[k] = -1
        ij_ds[k] = -1
    else
        i_c, i_d, Δ = _find_lon_indices_reactant(λ, lon_offsets[j + 1], nlons[j + 1])
        ij_cs[k] = ring_starts[j + 1] + i_c - 1
        ij_ds[k] = ring_starts[j + 1] + i_d - 1
        Δcds[k] = Δ
    end
end

# For ReactantDevice, override find_grid_indices! to convert λs to the element type of
# lon_offsets via Reactant.Ops.convert rather than a broadcast (broadcasting triggers
# Reactant JIT of `convert.(T, ConcretePJRTArray)` which fails because EnsureReturnType{T}
# rejects a TracedRNumber{T} return value from T(::TracedRNumber{S}))
function RingGrids.find_grid_indices!(
        locator::RingGrids.AbstractLocator,
        geometry::ReactantGridGeometry,
        λs::AbstractArray,
        architecture::ReactantDevice,
    )
    (; js, ij_as, ij_bs, ij_cs, ij_ds) = locator
    (; Δabs, Δcds) = locator
    (; nlons, lon_offsets, nlat) = geometry
    (; rings) = geometry

    T = eltype(lon_offsets)
    λs_T = eltype(λs) == T ? λs : Reactant.@jit(_cast_eltype(T, λs))

    # Flatten the tuple of UnitRanges to a plain Int array of per-ring starts.
    # Dynamic indexing into NTuple{nlat, UnitRange} breaks PTX codegen on GPU.
    # This runs once at interpolator setup, so overhead is negligible.
    ring_starts = Reactant.to_rarray(Int[first(r) for r in rings])

    launch!(
        architecture,
        LinearWorkOrder,
        size(js),
        find_grid_indices_kernel_reactant!,
        js,
        ij_as, ij_bs,
        ij_cs, ij_ds,
        Δabs, Δcds,
        λs_T,
        lon_offsets,
        nlons,
        nlat,
        ring_starts,
    )
    return nothing
end

# Override `average_on_poles` for Reactant: `Statistics.mean` iterates scalars on
# `ConcretePJRTArray` which is disallowed. Download the small polar ring slices to
# CPU and average there. This runs once during interpolator setup, so the CPU hop
# is negligible.
function RingGrids.average_on_poles(A::Reactant.AnyConcretePJRTArray, rings)
    north = @allowscalar A[rings[1]]
    south = @allowscalar A[rings[end]]
    return sum(north) / length(north), sum(south) / length(south)
end

# Dispatch find_rings_unsafe! for ReactantDevice to use the branchless kernel
function RingGrids.find_rings_unsafe!(
        js::AbstractArray,
        Δys::AbstractArray,
        θs::AbstractArray,
        latd::AbstractArray,
        architecture::ReactantDevice
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
        length(latd)
    )
end

end
