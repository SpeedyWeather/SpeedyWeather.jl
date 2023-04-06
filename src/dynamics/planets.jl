Base.@kwdef struct Earth <: Planet
    radius::Float64 = 6.371e6               # radius of Earth [m]
    rotation::Float64 = 7.29e-5             # angular frequency of Earth's rotation [rad/s]
    gravity::Float64 = 9.81                 # gravitational acceleration [m/s^2]
    
    daily_cycle::Bool = true                # daily cycle?
    length_of_day::Float64 = 24             # time [hrs] of a day

    seasonal_cycle::Bool = true             # Seasonal cycle?
    length_of_year::Float64 = 365.25        # time [days] of a year
    equinox::DateTime = DateTime(2000,3,20) # Spring equinox (year irrelevant)
    axial_tilt::Float64 = 23.4              # angle [˚] rotation axis tilt wrt to orbit
end