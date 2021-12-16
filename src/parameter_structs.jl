"""Coefficients of the generalised logistic function to describe the vertical coordinate.
Default coefficients A,K,C,Q,B,M,ν are fitted to the old L31 configuration at ECMWF.
See geometry.jl and function vertical_coordinate for more informaiton.

Following the notation of https://en.wikipedia.org/wiki/Generalised_logistic_function (Dec 15 2021).

Change default parameters for more/fewer levels in the stratosphere vs troposphere vs boundary layer."""
@with_kw struct GenLogisticCoefs
    A::Real=-0.283     # obtained from a fit in /input_date/vertical_coordinate/vertical_resolution.ipynb
    K::Real= 0.871
    C::Real= 0.414
    Q::Real= 6.695
    B::Real=10.336
    M::Real= 0.602
    ν::Real= 5.812
end