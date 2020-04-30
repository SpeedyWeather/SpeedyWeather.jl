struct Params
    n_diag::Int             # Frequency at which to output diagnostics (in number of time steps)
    n_steps_day::Int        # Number of time steps in a day
    Δt::Int                   # Time step in seconds
    # Coefficient for semi-implicit computations
    # 0 -> forward step for gravity wave terms,
    # 1 -> backward implicit, 0.5 -> centered implicit
    α::T
end

@with_kw struct Params
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
