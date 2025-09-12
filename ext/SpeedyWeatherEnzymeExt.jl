module SpeedyWeatherEnzymeExt

using SpeedyWeather
using Enzyme
using SpeedyWeather.ProgressMeter

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