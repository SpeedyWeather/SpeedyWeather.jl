module SpeedyWeatherEnzymeExt

using SpeedyWeather
using Enzyme
using Enzyme.EnzymeCore                 # Annotation, Const, ...
using Enzyme.EnzymeCore: Annotation
import Enzyme.EnzymeCore.EnzymeRules    # `EnzymeRules.forward` is extended below (qualified)
using SpeedyWeather.ProgressMeter

# currently not needed anymore:
#=
function __init__()
    # On Julia > 1.10, Enzyme's type analysis fails on large aggregates 
    # differentiated kernels (EnzymeNoTypeError, see https://github.com/EnzymeAD/Enzyme.jl/issues/3275).
    # We can fix this by setting a high maxtypeoffset (4096), but only on Julia 1.11+
    if VERSION >= v"1.11"
        Enzyme.API.maxtypeoffset!(4096)
    end
    return nothing
end
=# 

###
# implement make_zero where the default one fails

# this lock is part of the ProgressMeter that's part of the Feedback of all models
@inline function Enzyme.make_zero(
        ::Type{ProgressMeter.ProgressCore},
        seen::IdDict,
        prev::ProgressMeter.ProgressCore,
        ::Val{copy_if_inactive} = Val(false),
    )::ProgressMeter.ProgressCore where {copy_if_inactive}
    return prev
end

# View-preserving `make_zero` for the fused `Variables`.
#
# Some prognostic/grid leaves are `SubArray`-backed views into a shared parent that lives under
# `vars.fused.*`. Enzyme's default `make_zero` materialises those views into plain `Array`s, so the
# shadow's type no longer matches the view-backed primal. That mismatch breaks AD.
#
# Build the shadow as a `deepcopy` (which preserves the view→fused-parent aliasing and all types),
# then zero the differentiable data through the fuse parents / non-view leaves — the aliasing views
# are zeroed automatically and keep their `SubArray` type, matching the primal.
function Enzyme.make_zero(prev::SpeedyWeather.Variables)
    z = deepcopy(prev)
    for group in SpeedyWeather.ALL_VARIABLE_GROUPS
        _make_zero_view!(getfield(z, group))
    end
    return z
end
_make_zero_view!(x::NamedTuple) = foreach(k -> _make_zero_view!(getfield(x, k)), keys(x))
_make_zero_view!(x::SpeedyWeather.FusedParent) = _make_zero_view!(x.parent)   # the fuse buffer
_make_zero_view!(x::LowerTriangularArray) = SpeedyWeather.is_view_entry(x) || fill!(x.data, 0)
_make_zero_view!(x::SpeedyWeather.RingGrids.AbstractField) = SpeedyWeather.is_view_entry(x) || fill!(x.data, 0)
_make_zero_view!(::SubArray) = nothing
_make_zero_view!(x::AbstractArray) = eltype(x) <: Union{AbstractFloat, Complex} && fill!(x, 0)
# the transform scratch is a plain struct, so it needs its own methods to be reached at all
_make_zero_view!(x::SpeedyWeather.SpeedyTransforms.ScratchMemory) =
    (_make_zero_view!(x.north); _make_zero_view!(x.south); _make_zero_view!(x.column); nothing)
_make_zero_view!(x::SpeedyWeather.SpeedyTransforms.ColumnScratchMemory) =
    (_make_zero_view!(x.north); _make_zero_view!(x.south); nothing)
_make_zero_view!(x::Base.RefValue) = (x[] isa Number && (x[] = zero(x[])); nothing)
_make_zero_view!(x) = nothing

### 
# TODO
# Enzyme FORWARD-mode rule for iterating the tracer registry
#
# The time stepping iterates `model.tracers::TRACER_DICT` in 9 places (time_integration.jl,
# time_stepping/transform.jl, horizontal_diffusion.jl, vertical_advection.jl, tendencies.jl).
# Enzyme's forward mode fails with an internal error compiling `iterate(::Dict)` there, which
# blocks `autodiff(Forward, time_step!, ...)` for every model.
#
# This is a temporary workaround that works for differnetiating the model in forward mode without tracers
function EnzymeRules.forward(
        config::EnzymeRules.FwdConfig, func::Const{typeof(Base.iterate)}, ::Type{RT},
        dict::Annotation{<:SpeedyWeather.TRACER_DICT}, args::Annotation...,
    ) where {RT <: Annotation}
    primal = func.val(dict.val, map(a -> a.val, args)...)
    # the return holds only Symbol/Tracer, so no shadow is ever requested
    return EnzymeRules.needs_primal(config) ? primal : nothing
end

end
