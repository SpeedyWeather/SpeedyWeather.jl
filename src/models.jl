struct Model
    T::DataType
    params::Parameters
    constants::Constants
    geometry::Geometry
    current_datetime::DateTime
    end_datetime::DateTime
    spectral_trans::SpectralTrans
    boundaries::Boundaries
    prognostics::Prognostics
    horizontal_diffusion::HorizontalDiffusion
    implicit::Implicit
end

function Model(;
    # Model grid point resolution
    nlon, nlat,
    # Number of model levels
    nlev,
    # Model spectral resolution
    trunc,
    # Number of time steps per day
    n_steps_day,
    # Real number type
    real_type = Float64,
    # Radius of Earth
    Rₑ = 6.371e+6,
    # Angular frequency of Earth's rotation
    rotation_earth = 7.292e-05,
    # Gravitational acceleration
    g = 9.81,
    # Start and end datetimes
    start_datetime = DateTime(1982,1,1),
    end_datetime = DateTime(1982,1,2),
    # Name of restart file
    restart_file = nothing,
    # Number of tracers (including humidity)
    n_tr = 1,
    # Frequency of diagnostic checks
    n_diag = 36*5
    )

    # Default gas parameters
    akap = 2.0/7.0
    R    = akap*1000.4

    # Default dynamical constant parameters
    γ      = 6.0       # Reference temperature lapse rate (-dT/dz in deg/km)
    hscale = 7.5       # Reference scale height for pressure (in km)
    hshum  = 2.5       # Reference scale height for specific humidity (in km)
    refrh1 = 0.7       # Reference relative humidity of near-surface air

    current_datetime = start_datetime

    params = Parameters(n_diag, n_steps_day, 86400.0/real_type(n_steps_day), real_type(0.5),)
    constants = Constants(real_type, Rₑ, rotation_earth, g, akap, R, γ, hscale, hshum, refrh1)
    geometry = Geometry(real_type, constants, nlon, nlat, nlev, trunc)
    spectral_trans = SpectralTrans(real_type, geometry, constants.Rₑ)
    boundaries = Boundaries(real_type, geometry, spectral_trans, g)

    if restart_file == nothing
        prognostics = initialize_from_rest(real_type, geometry, constants, boundaries.ϕ₀ₛ,
                                           spectral_trans, n_tr)
       # Write initial data
       output(geometry, params, spectral_trans, current_datetime, 1, prognostics)
    else
        throw("Restart functionality not yet implemented")
    end

    horizontal_diffusion = HorizontalDiffusion(real_type, geometry, constants)
    implicit = Implicit(real_type, geometry, constants, params, horizontal_diffusion)

    Model(real_type, params, constants, geometry, current_datetime, end_datetime, spectral_trans, boundaries,
          prognostics, horizontal_diffusion, implicit)
end
