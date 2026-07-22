module RingGridsReactantExt

using RingGrids
using Reactant
using KernelAbstractions
using SpeedyWeatherInternals.Architectures: ReactantDevice
using SpeedyWeatherInternals.KernelLaunching: launch!, LinearWorkOrder

const AnyReactantArray = Union{Reactant.AnyConcretePJRTArray, Reactant.AnyConcreteIFRTArray, Reactant.AnyTracedRArray}
const ReactantGridGeometry = RingGrids.GridGeometry{<:AbstractGrid{<:ReactantDevice}}

# Route Field broadcasting through the underlying Reactant array's style so that
# broadcasts (e.g. `field .*= scale`) get JIT-compiled instead of iterating scalar.
# Must dispatch on `Type{F} where F <: AbstractField{...}` (not `Type{AbstractField{...}}`)
# to match the concrete `Field{...}` type Julia sees at a broadcast site.
Base.Broadcast.BroadcastStyle(::Type{<:AbstractField{T, N, ArrayType, Grid}}) where {T, N, ArrayType <: AnyReactantArray, Grid} = Base.Broadcast.BroadcastStyle(ArrayType)

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
# to obtain the ring index j such that latd[j] >= θ > latd[j+1]. This avoids
# data-dependent control flow (while / break / elseif) which is incompatible with
# Reactant tracing. The main-package `find_rings_kernel!` uses a binary search which
# is faster on CPU/GPU but cannot trace, hence this Reactant-only variant.
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
