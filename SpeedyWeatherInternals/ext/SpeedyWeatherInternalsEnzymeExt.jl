module SpeedyWeatherInternalsEnzymeExt

using SpeedyWeatherInternals
using SpeedyWeatherInternals.ParameterEditing: reconstruct
using Enzyme
using Enzyme.EnzymeCore
using Enzyme.EnzymeCore: Annotation
import Enzyme.EnzymeCore.EnzymeRules   # `EnzymeRules.forward` is extended below (qualified)

### Enzyme FORWARD-mode rule for `reconstruct`
#
# `reconstruct(obj, values)` cannot be differentiated in forward mode when `obj` is a full model:
# it is `@generated` and expands to a `setproperties` over the entire model type, so the layout
# Enzyme has to flatten (SpectralTransform with its scratch and FFTW plans, NetCDFOutput with its
# interpolator, Geometry, Feedback, ...) overruns Enzyme's type-analysis size budget and fails with
# `EnzymeNoTypeError`. Raising `Enzyme.API.maxtypeoffset!` clears it, but that is a global that has
# to be set before any differentiation and only moves the threshold as the model type grows.
#
# Reverse mode is unaffected: it already works at the default budget and gets no rule here.

# `b`-th tangent of an annotated argument. An inactive (Const) argument contributes a ZERO tangent,
# which for a structural scatter means a zeroed copy — using the primal here instead would leak
# primal values into the shadow's non-parameter fields and silently corrupt downstream tangents.
@inline _tangent(x::Union{Duplicated, DuplicatedNoNeed}, ::Int, _) = x.dval
@inline _tangent(x::Union{BatchDuplicated, BatchDuplicatedNoNeed}, b::Int, _) = x.dval[b]
@inline _tangent(::Const, ::Int, zeroed) = zeroed

function EnzymeRules.forward(
        config::EnzymeRules.FwdConfig, func::Const{typeof(reconstruct)}, ::Type{RT},
        obj::Annotation, values::Annotation,
    ) where {RT <: Annotation}

    primal = EnzymeRules.needs_primal(config) ? func.val(obj.val, values.val) : nothing
    EnzymeRules.needs_shadow(config) || return primal

    # zeroed stand-ins for inactive arguments, built once rather than per batch member
    # (NOTE: `make_zero` on a whole model is not free; it is only paid when `obj` is Const)
    zobj = obj isa Const ? Enzyme.make_zero(obj.val) : nothing
    zval = values isa Const ? Enzyme.make_zero(values.val) : nothing

    B = EnzymeRules.width(config)
    shadow = ntuple(b -> func.val(_tangent(obj, b, zobj), _tangent(values, b, zval)), B)

    EnzymeRules.needs_primal(config) || return B == 1 ? shadow[1] : shadow
    return B == 1 ? Duplicated(primal, shadow[1]) : BatchDuplicated(primal, shadow)
end

end
