module SpeedyWeatherEnzymeExt

using SpeedyWeather
using Enzyme
using SpeedyWeather.ProgressMeter

function __init__()
    # On Julia > 1.10, Enzyme's type analysis cannot prove the element type of a prognostic step
    # view (`get_step`) once it is reached inside a differentiated kernel, throwing an
    # `EnzymeNoTypeError`. (Increasing the type-analysis size/depth via maxtypeoffset!/maxtypedepth!
    # does not resolve it — it is an unprovable-type issue, not a search-space one.) SpeedyWeather's
    # differentiated arrays are homogeneous (all Float32 / ComplexF32), so letting Enzyme take its
    # best guess is safe here and keeps the extensively-used `get_step`/`get_*_step` views and the
    # kernels untouched. Julia 1.10 resolves these types natively and needs no change.
    if VERSION >= v"1.11"
        Enzyme.API.looseTypeAnalysis!(true)
    end
    return nothing
end

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

end
