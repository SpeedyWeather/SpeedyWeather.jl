abstract type Coefficients end

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

"""
Parameters for computing saturation vapour pressure using the August-Roche-Magnus formula,

    eᵢ(T) = e₀ * exp(Cᵢ * (T - T₀) / (T - Tᵢ)),

    where T is in Kelvin and i = 1,2 for saturation with respect to water and ice,
    respectively.
"""
Base.@kwdef struct MagnusCoefs{NF<:Real} <: Coefficients
    e₀::NF = 6.108   # Saturation vapour pressure at 0°C
    T₀::NF = 273.16  # 0°C in Kelvin
    T₁::NF = 35.86
    T₂::NF = 7.66
    C₁::NF = 17.269
    C₂::NF = 21.875
end

"""
Parameters for radiation parameterizations.
"""
Base.@kwdef struct RadiationCoefs{NF<:Real} <: Coefficients
    epslw::NF = 0.05    # Fraction of blackbody spectrum absorbed/emitted by PBL only
    emisfc::NF = 0.98   # Longwave surface emissivity

    p0::Real = 1e5      # Reference pressure (Pa)

    # Shortwave radiation: sol_oz
    solc::NF = 342.0    # Solar constant (area averaged) [W/m^2]
    epssw::NF = 0.020   # Fraction of incoming solar radiation absorbed by ozone

    # Shortwave radiation: cloud
    rhcl1::NF = 0.30    # Relative humidity threshold corresponding to cloud cover = 0
    rhcl2::NF = 1.00    # Relative humidity correponding to cloud cover = 1
    rrcl::NF = 1 / (rhcl2 - rhcl1)
    qcl::NF = 0.20      # Specific humidity threshold for cloud cover
    pmaxcl::NF = 10.0   # Maximum value of precipitation (mm/day) contributing to cloud cover
    wpcl::NF = 0.2      # Cloud cover weight for the square-root of precipitation (for p = 1 mm/day)
    gse_s1::NF = 0.40   # Gradient of dry static energy corresponding to stratiform cloud cover = 1
    gse_s0::NF = 0.25   # Gradient of dry static energy corresponding to stratiform cloud cover = 0
    clsmax::NF = 0.60   # Maximum stratiform cloud cover
    clsminl::NF = 0.15  # Minimum stratiform cloud cover over land (for RH = 1)

    # Shortwave radiation: radsw
    albcl::NF = 0.43    # Cloud albedo (for cloud cover = 1)
    albcls::NF = 0.50   # Stratiform cloud albedo (for st. cloud cover = 1)
    abscl1::NF = 0.015  # Absorptivity of clouds (visible band, maximum value)
    abscl2::NF = 0.15   # Absorptivity of clouds (visible band, for dq_base = 1 g/kg)

    absdry::NF = 0.033  # Absorptivity of dry air (visible band)
    absaer::NF = 0.033  # Absorptivity of aerosols (visible band)
    abswv1::NF = 0.022  # Absorptivity of water vapour (visible band, for dq = 1 g/kg)
    abswv2::NF = 15.0   # Absorptivity of water vapour (near IR band, for dq = 1 g/kg)

    ablwin::NF = 0.3    # Absorptivity of air in "window" band
    ablco2::NF = 6.0    # Absorptivity of air in CO2 band
    ablwv1::NF = 0.7    # Absorptivity of water vapour in H2O band 1 (weak), (for dq = 1 g/kg)
    ablwv2::NF = 50.0   # Absorptivity of water vapour in H2O band 2 (strong), (for dq = 1 g/kg)
    ablcl2::NF = 0.6    # Absorptivity of "thin" upper clouds in window and H2O bands
    ablcl1::NF = 12.0   # Absorptivity of "thick" clouds in window band (below cloud top)
end    
