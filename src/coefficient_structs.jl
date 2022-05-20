abstract type Coefficients end

"""Coefficients of the generalised logistic function to describe the vertical coordinate.
Default coefficients A,K,C,Q,B,M,ν are fitted to the old L31 configuration at ECMWF.
See geometry.jl and function vertical_coordinate for more informaiton.

Following the notation of https://en.wikipedia.org/wiki/Generalised_logistic_function (Dec 15 2021).

Change default parameters for more/fewer levels in the stratosphere vs troposphere vs boundary layer."""
@with_kw struct GenLogisticCoefs{NF<:Real} <: Coefficients
    A::NF = -0.283     # obtained from a fit in /input_date/vertical_coordinate/vertical_resolution.ipynb
    K::NF = 0.871
    C::NF = 0.414
    Q::NF = 6.695
    B::NF = 10.336
    M::NF = 0.602
    ν::NF = 5.812
end

"""
Parameters for computing saturation specific humidity.

These parameters, and the formula where they are used, do not appear in the original
Speedy documentation.
"""
@with_kw struct HumidityCoefs{NF<:Real} <: Coefficients
    e₀::NF = 6.108e-3
    C₁::NF = 17.269
    C₂::NF = 21.875
    T₀::NF = 273.16
    T₁::NF = 35.86
    T₂::NF = 7.66
end
