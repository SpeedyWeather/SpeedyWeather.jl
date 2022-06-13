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
Parameters for computing saturation vapour pressure using the August-Roche-Magnus formula,

    eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T - Tᵢ)),

    where T is in Kelvin and i = 1,2 for saturation with respect to water and ice,
    respectively.
"""
@with_kw struct MagnusCoefs{NF<:Real} <: Coefficients
    e₀::NF = 6.108   # Saturation vapour pressure at 0°C
    T₀::NF = 273.16  # 0°C in Kelvin
    T₁::NF = 35.86
    T₂::NF = 7.66
    C₁::NF = 17.269
    C₂::NF = 21.875
end
