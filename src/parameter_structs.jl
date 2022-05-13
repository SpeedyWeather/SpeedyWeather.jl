abstract type Parameters end

"""Coefficients of the generalised logistic function to describe the vertical coordinate.
Default coefficients A,K,C,Q,B,M,ν are fitted to the old L31 configuration at ECMWF.
See geometry.jl and function vertical_coordinate for more informaiton.

Following the notation of https://en.wikipedia.org/wiki/Generalised_logistic_function (Dec 15 2021).

Change default parameters for more/fewer levels in the stratosphere vs troposphere vs boundary layer."""
@with_kw struct GenLogisticCoefs <: Parameters
    A::Real = -0.283     # obtained from a fit in /input_date/vertical_coordinate/vertical_resolution.ipynb
    K::Real = 0.871
    C::Real = 0.414
    Q::Real = 6.695
    B::Real = 10.336
    M::Real = 0.602
    ν::Real = 5.812
end

"""
Parameters for computing saturation specific humidity.

These parameters, and the formula where they are used, do not appear in the original
documentation.
"""
@with_kw struct HumidityCoefs <: Parameters
    e₀::Real = 6.108e-3
    C₁::Real = 17.269
    C₂::Real = 21.875
    T₀::Real = 273.16
    T₁::Real = 35.86
    T₂::Real = 7.66
end
