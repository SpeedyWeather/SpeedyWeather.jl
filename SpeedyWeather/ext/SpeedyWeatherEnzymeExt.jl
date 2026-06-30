module SpeedyWeatherEnzymeExt

using SpeedyWeather
using Enzyme
using Enzyme.EnzymeCore
import .EnzymeRules: augmented_primal, reverse
using .EnzymeRules
using SpeedyWeather.ProgressMeter

import SpeedyWeather: get_step

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

###
# Custom reverse rule for `get_step(var, step)`.
#
# `get_step(var, step)` returns a *view* into `var` selecting the leapfrog/time step `step`
# along the last dimension. With a runtime `step` this view carries a runtime offset, and on
# Julia ≥ 1.11 Enzyme's type analysis can fail to type the view construction once the view
# participates in a differentiated kernel (`EnzymeNoTypeError`). Making `get_step` a custom-rule
# boundary sidesteps that: we hand Enzyme the primal view and, for a `Duplicated` input, the
# *same* view into the shadow array. Because that shadow view aliases `var`'s shadow, cotangents
# written through it already land in `var.dval`, so the reverse pass has nothing to accumulate.
#
# Keeping the primal `get_step` branchless (no step-to-literal ladder) is also what we want for
# Reactant, which traces the plain view and does not use these Julia-level rules.

function augmented_primal(
        config::EnzymeRules.RevConfigWidth{1},
        func::Const{typeof(get_step)},
        ::Type{RT},
        var::Annotation{<:AbstractArray},
        step::Const,
    ) where {RT}
    primal = func.val(var.val, step.val)
    shadow = var isa Duplicated ? func.val(var.dval, step.val) : nothing
    primal_return = EnzymeRules.needs_primal(config) ? primal : nothing
    return EnzymeRules.AugmentedReturn(primal_return, shadow, nothing)
end

function reverse(
        config::EnzymeRules.RevConfigWidth{1},
        func::Const{typeof(get_step)},
        ::Type{RT},
        tape,
        var::Annotation{<:AbstractArray},
        step::Const,
    ) where {RT}
    # the returned view aliases `var`; its shadow aliases `var.dval`, so derivatives have
    # already been accumulated there. Nothing to add for either argument (var, step).
    return (nothing, nothing)
end

end
