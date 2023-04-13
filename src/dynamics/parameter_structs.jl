"""Coefficients of the generalised logistic function to describe the vertical coordinate.
Default coefficients A,K,C,Q,B,M,ν are fitted to the old L31 configuration at ECMWF.
See geometry.jl and function `vertical_coordinate` for more informaiton.

Following the notation of [https://en.wikipedia.org/wiki/Generalised_logistic_function](https://en.wikipedia.org/wiki/Generalised_logistic_function) (Dec 15 2021).

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

struct StartFromFile <: InitialConditions end
struct StartFromRest <: InitialConditions end
struct StartWithVorticity <: InitialConditions end

"""
    Z = ZonalJet(;kwargs...) <: InitialConditions

Create a struct that contains all parameters for the Galewsky et al, 2004 zonal jet
intitial conditions for the shallow water model. Default values as in Galewsky."""
Base.@kwdef struct ZonalJet <: InitialConditions
    # jet
    latitude = 45               # degrees north [˚N]
    width = (1/4-1/7)*180       # ≈ 19.29˚ as in Galewsky et al, 2004 
    umax = 80                   # [m/s]
    
    # perturbation
    perturb_lat = latitude          # [˚N], position in jet by default
    perturb_lon = 270               # [˚E]
    perturb_xwidth = 1/3*360/2π     # ≈ 19.1˚E zonal extent [˚E]
    perturb_ywidth = 1/15*360/2π    # ≈ 3.8˚N meridional extent [˚N]
    perturb_height = 120            # amplitude [m]
end

"""
    Z = ZonalWind(;kwargs...) <: InitialConditions

Create a struct that contains all parameters for the Jablonowski and Williamson, 2006
intitial conditions for the primitive equation model. Default values as in Jablonowski."""
Base.@kwdef struct ZonalWind <: InitialConditions
    
    # vertical
    η₀ = 0.252                  # conversion from σ to Jablonowski's ηᵥ-coordinates
    u₀ = 35                     # max amplitude of zonal wind [m/s]
    
    # perturbation
    perturb_lat = 40            # Gaussian profile perturbation centred at [˚N]
    perturb_lon = 20            # and [˚E]
    perturb_uₚ = 1              # strength of perturbation [m/s]
    perturb_radius = 1/10       # radius of Gaussian perturbation in units of Earth's radius [1]

    # temperature
    ΔT = 4.8e5                  # temperature difference used for stratospheric lapse rate [K]
end

