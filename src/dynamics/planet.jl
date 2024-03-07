abstract type AbstractPlanet <: AbstractModelComponent end

const DEFAULT_ROTATION = 7.29e-5    # default angular frequency of Earth's rotation [1/s]
const DEFAULT_GRAVITY = 9.81        # default gravitational acceleration on Earth [m/s²]

export Earth

"""
$(TYPEDSIGNATURES)
Create a struct `Earth<:AbstractPlanet`, with the following physical/orbital
characteristics. Note that `radius` is not part of it as this should be chosen
in `SpectralGrid`. Keyword arguments are
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Earth{NF<:AbstractFloat} <: AbstractPlanet

    "angular frequency of Earth's rotation [rad/s]"
    rotation::NF = DEFAULT_ROTATION
    
    "gravitational acceleration [m/s^2]"
    gravity::NF = DEFAULT_GRAVITY              
    
    "switch on/off daily cycle"
    daily_cycle::Bool = true
    
    "Seconds in a daily rotation"
    length_of_day::Second = Hour(24)             

    "switch on/off seasonal cycle"
    seasonal_cycle::Bool = true

    "Seconds in an orbit around the sun"
    length_of_year::Second = Day(365.25)
    
    "time of spring equinox (year irrelevant)"
    equinox::DateTime = DateTime(2000,3,20) 

    "angle [˚] rotation axis tilt wrt to orbit"
    axial_tilt::NF = 23.4

    "Total solar irradiance at the distance of 1 AU [W/m²]"
    solar_constant::NF = 1365
end

Earth(;kwargs...) = Earth{DEFAULT_NF}(;kwargs...)
Earth(SG::SpectralGrid;kwargs...) = Earth{SG.NF}(;kwargs...)
Earth(::Type{NF};kwargs...) where NF = Earth{NF}(;kwargs...)