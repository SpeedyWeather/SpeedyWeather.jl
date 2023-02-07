"""Coefficients of the generalised logistic function to describe the vertical coordinate.
Default coefficients A,K,C,Q,B,M,ν are fitted to the old L31 configuration at ECMWF.
See geometry.jl and function vertical_coordinate for more informaiton.

Following the notation of https://en.wikipedia.org/wiki/Generalised_logistic_function (Dec 15 2021).

Change default parameters for more/fewer levels in the stratosphere vs troposphere vs boundary layer."""
Base.@kwdef struct GenLogisticCoefs <: Coefficients
    A::Float64 = -0.283     # obtained from a fit in /input_date/vertical_coordinate/vertical_resolution.ipynb
    K::Float64 = 0.871
    C::Float64 = 0.414
    Q::Float64 = 6.695
    B::Float64 = 10.336
    M::Float64 = 0.602
    ν::Float64 = 5.812
end

"""Generalised logistic function based on the coefficients in `coefs`."""
function generalised_logistic(x,coefs::GenLogisticCoefs)
    @unpack A,K,C,Q,B,M,ν = coefs
    return @. A + (K-A)/(C+Q*exp(-B*(x-M)))^inv(ν)
end