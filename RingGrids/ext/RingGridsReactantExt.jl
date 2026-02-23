module RingGridsReactantExt

using RingGrids
using Reactant
using KernelAbstractions
using SpeedyWeatherInternals.Architectures: ReactantDevice
using SpeedyWeatherInternals.Utils: launch!, LinearWorkOrder

const AnyReactantArray = Union{Reactant.AnyConcretePJRTArray, Reactant.AnyConcreteIFRTArray, Reactant.AnyTracedRArray}

Base.Broadcast.BroadcastStyle(::Type{AbstractField{T, N, ArrayType, Grid}}) where {T, N, ArrayType <: AnyReactantArray, Grid} = Base.Broadcast.BroadcastStyle(ArrayType)

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
