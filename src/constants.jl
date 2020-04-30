struct Constants{T<:AbstractFloat}
    R_earth::T      # Radius of Earth
    Ω::T            # Angular frequency of Earth's rotation
    g::T            # Gravitational acceleration
    akap::T         # Ratio of gas constant to specific heat of dry air at constant pressure
    R::T            # Gas constant
    γ::T            # Reference temperature lapse rate (-dT/dz in deg/km)
    hscale::T       # Reference scale height for pressure (in km)
    hshum::T        # Reference scale height for specific humidity (in km)
    refrh1::T       # Reference relative humidity of near-surface air
end

@with_kw struct Constants{T<:AbstractFloat}
    R_earth::T      # Radius of Earth
    Ω::T            # Angular frequency of Earth's rotation
    g::T            # Gravitational acceleration
    akap::T         # Ratio of gas constant to specific heat of dry air at constant pressure
    R::T            # Gas constant
    γ::T            # Reference temperature lapse rate (-dT/dz in deg/km)
    hscale::T       # Reference scale height for pressure (in km)
    hshum::T        # Reference scale height for specific humidity (in km)
    refrh1::T       # Reference relative humidity of near-surface air
end
