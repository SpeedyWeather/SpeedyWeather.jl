"""
Parameter struct that holds all parameters that define the default model setup.
With keywords such that default values can be changed at creation.
"""
@with_kw struct Parameters

    # NUMBER FORMATS
    NF::DataType=Float64    # number format

    # RESOLUTION
    trunc::Int=31                      # spectral truncation
    nlon::Int=roundup_fft(3*trunc+1)   # number of longitudes 
    nlat::Int=nlon÷2                   # number of latitudes
    
    nlev::Int=8             # number of vertical levels
    ntracers::Int=1         # number of tracers (specific humidity is one)

    # PHYSICAL CONSTANTS
    R_earth::Real=6.371e6   # Radius of Earth [m]
    Ω::Real=7.292e-5        # Angular frequency of Earth's rotation [1/s]
    gravity::Real=9.81      # Gravitational acceleration [m/s^2]
    akap::Real=2/7          # Ratio of gas constant to specific heat of dry air
                            # at constant pressure = 1 - 1/γ where γ is the
                            # heat capacity ratio of a perfect diatomic gas (7/5)
    cp::Real=1004.0         # Specific heat at constant pressure [J/K/kg]
    R_gas::Real=akap*cp     # Specific gas constant for dry air [J/kg/K]
    alhc::Real=2501.0       # Latent heat of condensation [J/g] for consistency with
                            # specific humidity [g/Kg]
    alhs::Real=2801.0       # Latent heat of sublimation [?]
    sbc::Real=5.67e-8       # Stefan-Boltzmann constant [W/m^2/K^4]

    # STANDARD ATMOSPHERE
    lapse_rate::Real=6.0    # Reference temperature lapse rate -dT/dz [K/km]
    temp_ref::Real=288      # Reference absolute temperature at surface z=0 [K]
    temp_top::Real=216      # Reference absolute temperature in stratosphere [K]
    scale_height::Real=7.5  # Reference scale height for pressure [km]
    pres_ref::Real=1013     # Reference surface pressure [hPa]
    scale_height_humid::Real=2.5    # Reference scale height for specific humidity [km]
    relhumid_ref::Real=0.7          # Reference relative humidity of near-surface air [1]
    water_pres_ref::Real=17         # Reference saturation water vapour pressure [Pa]

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
    seasonal_cycle::Bool=true   # Seasonal cycle?
    n_shortwave::Int=3          # Compute shortwave radiation every n steps
    sppt_on::Bool=false         # Turn on SPPT?

    # TIME STEPPING
    nstepsday::Int=36           # number of time steps in one day
    Δt::Real=86400/nstepsday    # time step in seconds
    ndays::Real=10              # number of days to integrate for

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
    initial_conditions::Symbol=:rest    # :rest or :restart

    # OUTPUT
    verbose::Bool=true          # print dialog for feedback
    output::Bool=false          # Store data in netCDF?
    output_dt::Real=6           # output time step [hours]
    output_startdate::DateTime=DateTime(2000,1,1)
    outpath::String=pwd()       # path to output folder
    output_vars::Vector{String}=["u","v","T","humid","logp0"]
    compression_level::Int=3    # 1=low but fast, 9=high but slow
    keepbits::Int=10            # mantissa bits to keep for every variable 

    # TODO assert not allowed parameter values
    @assert α in [0,0.5,1] "Only semi-implicit α = 0, 0.5 or 1 allowed."
end