# struct Constants{T<:AbstractFloat}
#     R_earth::T      # Radius of Earth
#     Ω::T            # Angular frequency of Earth's rotation
#     g::T            # Gravitational acceleration
#     akap::T         # Ratio of gas constant to specific heat of dry air at constant pressure
#     R::T            # Gas constant
#     γ::T            # Reference temperature lapse rate (-dT/dz in deg/km)
#     hscale::T       # Reference scale height for pressure (in km)
#     hshum::T        # Reference scale height for specific humidity (in km)
#     refrh1::T       # Reference relative humidity of near-surface air
# end

@with_kw struct Constants
    R_earth::Real=6.371e6       # Radius of Earth [m]
    Ω::Real=7.292e-5            # Angular frequency of Earth's rotation [1/s]
    g::Real=9.81                # Gravitational acceleration [m/s^2]
    akap::Real=2/7              # Ratio of gas constant to specific heat of dry air at constant pressure
    R::Real=287.058             # Specific gas constant for dry air [J/kg/K]
    γ::Real=6.0                 # Reference temperature lapse rate -dT/dz [deg/km]
    hscale::Real=7.5            # Reference scale height for pressure (in km)
    hshum::Real=2.5             # Reference scale height for specific humidity (in km)
    refrh1::Real=0.7            # Reference relative humidity of near-surface air
end
