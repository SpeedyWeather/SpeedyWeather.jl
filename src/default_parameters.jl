"""
Parameter struct that holds all parameters that define the default model setup.
With keywords such that default values can be changed at creation.
"""
@with_kw struct Parameters

    # NUMBER FORMATS
    NF::DataType=Float32        # number format

    # RESOLUTION
    trunc::Int=31                       # spectral truncation
    nlon::Int=roundup_fft(3*trunc+1)    # number of longitudes
    nlat::Int=nlon÷2                    # number of latitudes
    nlev::Int=8                         # number of vertical levels

    # PHYSICAL CONSTANTS
    model::Symbol=:barotropic   # :barotropic, :shallowwater, or :primitive
    radius_earth::Real=6.371e6          # radius of Earth [m]
    rotation_earth::Real=7.292e-5       # angular frequency of Earth's rotation [1/s]
    gravity::Real=9.81          # gravitational acceleration [m/s^2]
    akap::Real=2/7              # ratio of gas constant to specific heat of dry air
                                # at constant pressure = 1 - 1/γ where γ is the
                                # heat capacity ratio of a perfect diatomic gas (7/5)
    cp::Real=1004               # specific heat at constant pressure [J/K/kg]
    R_gas::Real=akap*cp         # specific gas constant for dry air [J/kg/K]
    alhc::Real=2501             # latent heat of condensation [J/g] for consistency with
                                # specific humidity [g/Kg]
    alhs::Real=2801             # latent heat of sublimation [?]
    sbc::Real=5.67e-8           # stefan-Boltzmann constant [W/m^2/K^4]

    # STANDARD ATMOSPHERE
    lapse_rate::Real=6          # reference temperature lapse rate -dT/dz [K/km]
    temp_ref::Real=288          # reference absolute temperature at surface z=0 [K]
    temp_top::Real=216          # reference absolute temperature in stratosphere [K]
    scale_height::Real=7.5      # reference scale height for pressure [km]
    pres_ref::Real=1013         # reference surface pressure [hPa]
    scale_height_humid::Real=2.5# reference scale height for specific humidity [km]
    relhumid_ref::Real=0.7      # reference relative humidity of near-surface air [1]
    water_pres_ref::Real=17     # reference saturation water vapour pressure [Pa]

    # VERTICAL COORDINATES
    # of the nlev vertical levels, defined by a generalised logistic function,
    # interpolating ECMWF's L31 configuration
    GLcoefs::GenLogisticCoefs=GenLogisticCoefs()
    n_stratosphere_levels::Int=2    # number of vertical levels used for the stratosphere

    # DIFFUSION AND DRAG
    diffusion_power::Real=4                 # Power n of Laplacian in horizontal diffusion ∇²ⁿ
    diffusion_time::Real=2.4                # Diffusion time scale [hrs] for temperature and vorticity
    diffusion_time_div::Real=diffusion_time # Diffusion time scale [hrs] for divergence
    diffusion_time_strat::Real=12           # Diffusion time scale [hrs] for extra ∇² in the stratosphere
    damping_time_strat::Real=24*30          # Damping time [hrs] for drag on zonal-mean wind in the stratosphere

    # PARAMETRIZATIONS
    seasonal_cycle::Bool=true                  # Seasonal cycle?
    n_shortwave::Int=3                         # Compute shortwave radiation every n steps
    sppt_on::Bool=false                        # Turn on SPPT?
    magnus_coefs::MagnusCoefs = MagnusCoefs()  # For computing saturation vapour pressure

    # Large-scale condensation (from table B10)
    RH_thresh_boundary::Real = 0.95  # Relative humidity threshold for boundary layer
    RH_thresh_range::Real = 0.1      # Vertical range of relative humidity threshold
    RH_thresh_max::Real = 0.9        # Maximum relative humidity threshold
    humid_relax_time::Real = 4.0     # Relaxation time for humidity (hours)

    # TIME STEPPING
    Δt_at_T85::Real=20          # time step in minutes for T85, scale linearly for specified trunc
    n_days::Real=10             # number of days to integrate for

    # NUMERICS
    robert_filter::Real=0.05    # Robert (1966) time filter coefficeint for suppress comput. mode
    williams_filter::Real=0.53  # William's time filter (Amezcua 2011) coefficient for 3rd order acc
    α::Real=0.5                 # coefficient for semi-implicit computations
                                # 0 -> forward step for gravity wave terms,
                                # 1 -> backward implicit
                                # 0.5 -> centered implicit
    recompute_legendre::Bool=false

    # BOUNDARY FILES
    boundary_path::String=""    # package location is default
    orography_path::String=boundary_path
    orography_file::String="orography_F512.nc"

    # INITIAL CONDITIONS
    initial_conditions::Symbol=:barotropic_vorticity    # :test, :rest, :barotropic_vorticity or :restart

    # OUTPUT
    verbose::Bool=true          # print dialog for feedback
    output::Bool=false          # Store data in netCDF?
    output_dt::Real=6           # output time step [hours]
    output_startdate::DateTime=DateTime(2000,1,1)
    out_path::String=pwd()      # path to output folder
    output_vars::Vector{String}=["u","v","T","humid","logp0"]
    compression_level::Int=3    # 1=low but fast, 9=high but slow
    keepbits::Int=10            # mantissa bits to keep for every variable

    # TODO assert not allowed parameter values
    @assert α in [0,0.5,1] "Only semi-implicit α = 0, 0.5 or 1 allowed."
end
