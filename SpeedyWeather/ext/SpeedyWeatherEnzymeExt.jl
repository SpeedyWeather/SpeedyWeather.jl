module SpeedyWeatherEnzymeExt

using SpeedyWeather
using Enzyme
using SpeedyWeather.ProgressMeter

function __init__()
    # On Julia > 1.10, Enzyme's type analysis fails on large aggregates 
    # differentiated kernels (EnzymeNoTypeError, see https://github.com/EnzymeAD/Enzyme.jl/issues/3275).
    # We can fix this by setting a high maxtypeoffset (4096), but only on Julia 1.11+
    if VERSION >= v"1.11"
        Enzyme.API.maxtypeoffset!(4096)
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
