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

###
# Enzyme FORWARD-mode rule for `reconstruct` on a model
#
# `reconstruct(model, p)` cannot be differentiated in forward mode: it is `@generated` and expands
# to a `setproperties` over the entire model type, so the layout Enzyme has to flatten
# the very complex model structre. Raising `Enzyme.API.maxtypeoffset!` also clears it (8192 is needed for
# `PrimitiveWetModel`, 2048 is not), but that is a global that must be set before any
# differentiation and only moves a threshold that drifts as the model type grows.
#
# The rule removes the dependence on type size, and is exact because `reconstruct` does NO
# arithmetic: it is a pure structural scatter, where parameter-valued fields of the output come
# from `values` and every other field comes from `obj`. It is therefore linear, and its tangent is
# the same scatter over the tangents:
#
#     d/dε reconstruct(obj + ε·dobj, values + ε·dvalues) = reconstruct(dobj, dvalues)
#
# RESTRICTED TO `AbstractModel` ON PURPOSE. Enzyme aborts the process
# ("Assertion failed: needsReRooting", FixupJuliaCallingConvention) one some comments with the rule (bug)

# `b`-th tangent of an annotated argument. An inactive (Const) argument contributes a ZERO tangent,
# which for a structural scatter means a zeroed copy — reusing the primal instead would leak primal
# values into the shadow's non-parameter fields and corrupt downstream tangents.
@inline _tangent(x::Union{Duplicated, DuplicatedNoNeed}, ::Int, _) = x.dval
@inline _tangent(x::Union{BatchDuplicated, BatchDuplicatedNoNeed}, b::Int, _) = x.dval[b]
@inline _tangent(::Const, ::Int, zeroed) = zeroed

function EnzymeRules.forward(
        config::EnzymeRules.FwdConfig, func::Const{typeof(SpeedyWeather.reconstruct)}, ::Type{RT},
        obj::Annotation{<:SpeedyWeather.AbstractModel}, values::Annotation,
    ) where {RT <: Annotation}

    # Each return path is written out separately, and `primal` is only bound on the path that uses
    # it: binding it up front as `needs_primal(config) ? ... : nothing` would give it type
    # `Union{Nothing, T}` and feed that Union into `Duplicated`.
    if !EnzymeRules.needs_shadow(config)
        return EnzymeRules.needs_primal(config) ? func.val(obj.val, values.val) : nothing
    end

    # zeroed stand-ins for inactive arguments, built once rather than per batch member
    # (NOTE: `make_zero` on a whole model is not free; it is only paid when `obj` is Const)
    zobj = obj isa Const ? Enzyme.make_zero(obj.val) : nothing
    zval = values isa Const ? Enzyme.make_zero(values.val) : nothing

    shadow = ntuple(
        b -> func.val(_tangent(obj, b, zobj), _tangent(values, b, zval)),
        Val(EnzymeRules.width(config)),      # Val so the width stays a compile-time constant
    )

    if !EnzymeRules.needs_primal(config)
        return EnzymeRules.width(config) == 1 ? only(shadow) : shadow
    end
    primal = func.val(obj.val, values.val)
    return EnzymeRules.width(config) == 1 ? Duplicated(primal, only(shadow)) :
        BatchDuplicated(primal, shadow)
end

### REVERSE-mode rule for `reconstruct` on a model
#
# Needed for NESTED AD (reverse over forward): once a function carries a custom rule, Enzyme looks
# one up for whichever mode is running, so the forward rule above has to be matched by a reverse
# one. Plain single-mode reverse over `reconstruct` works without any rule and is unchanged by this.
#
# The pullback follows from `reconstruct` being a structural scatter. For each field of the output:
#   * a parameter-valued field comes from `values` -> its cotangent belongs to `values`
#   * every other field comes from `obj`           -> its cotangent belongs to `obj`
#
# The non-parameter fields are copied BY REFERENCE (`setproperties` shares the arrays), so building
# the output shadow as `reconstruct(obj.dval, <zeroed values>)` makes those fields alias `obj`'s own
# shadow: anything the caller accumulates into them lands in `obj.dval` directly, with no explicit
# adjoint accumulation and no extra allocation. That leaves only the scalar parameter slots, which
# are stored by value and so must be gathered explicitly in `reverse` — `parameters(shadow)` reads
# exactly those slots back out, in the same layout as `values`.
function EnzymeRules.augmented_primal(
        config::EnzymeRules.RevConfig, func::Const{typeof(SpeedyWeather.reconstruct)}, ::Type{RT},
        obj::Annotation{<:SpeedyWeather.AbstractModel}, values::Annotation,
    ) where {RT <: Annotation}

    primal = func.val(obj.val, values.val)

    shadow = if EnzymeRules.needs_shadow(config)
        # alias obj's shadow for the non-parameter fields (see above); zero the parameter slots so
        # they start from zero and accumulate
        base = obj isa Const ? Enzyme.make_zero(obj.val) : obj.dval
        func.val(base, Enzyme.make_zero(values.val))
    else
        nothing
    end

    # the shadow is also the tape: `reverse` reads the accumulated cotangent back out of it
    return EnzymeRules.AugmentedReturn(
        EnzymeRules.needs_primal(config) ? primal : nothing, shadow, shadow,
    )
end

function EnzymeRules.reverse(
        ::EnzymeRules.RevConfig, ::Const{typeof(SpeedyWeather.reconstruct)}, ::Type{RT}, tape,
        obj::Annotation{<:SpeedyWeather.AbstractModel}, values::Annotation,
    ) where {RT <: Annotation}

    # gather the parameter slots of the output cotangent into the `values` cotangent. The obj
    # cotangent needs nothing here: its non-parameter fields alias the shadow built above.
    if !(values isa Const) && tape !== nothing
        values.dval .+= vec(SpeedyWeather.parameters(tape))
    end

    return (nothing, nothing)
end

end
