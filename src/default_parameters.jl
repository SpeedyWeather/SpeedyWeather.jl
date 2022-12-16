"""
    P = Parameters(kwargs...)

A struct to hold all model parameters that may be changed by the user.
The struct uses keywords such that default values can be changed at creation.
The default values of the keywords define the default model setup.
"""
@with_kw struct Parameters

    # NUMBER FORMATS
    NF::DataType                        # number format (default defined in run_speedy)

    # MODEL (BarotropicModel, ShallowWaterModel, or PrimitiveEquationModel)
    model::Type{<:ModelSetup} = Barotropic
    dry_core::Bool = true               # true = switches off humidity calculations

    # RESOLUTION AND GRID
    trunc::Int = 31                                     # spectral truncation
    Grid::Type{<:AbstractGrid} = FullGaussianGrid       # grid used
    
    # EARTH'S PROPERTIES
    radius_earth::Real = 6.371e6            # radius of Earth [m]
    rotation_earth::Real = 7.29e-5          # angular frequency of Earth's rotation [rad/s]
    gravity::Real = 9.81                    # gravitational acceleration [m/s^2]
    seasonal_cycle::Bool = true             # Seasonal cycle?
    equinox::DateTime = DateTime(2000,3,20) # Spring equinox (year irrelevant)
    tropic_cancer::Real = 23.5              # latitude [˚N] of the tropic of cancer
    
    # ATMOSPHERE
    mol_mass_dry_air = 28.964917        # molar mass of dry air [g/mol]
    mol_mass_vapour = 18.01528          # molar mass of water vapour [g/mol]
    cₚ::Real = 1004.0                   # specific heat at constant pressure [J/K/kg]
    R_gas::Real = 8.31446261815324      # universal gas constant [J/K/mol]
    R_dry::Real = 1000*R_gas/mol_mass_dry_air   # specific gas constant for dry air [J/kg/K]
    R_vapour::Real = 1000*R_gas/mol_mass_vapour # specific gas constant for water vapour [J/kg/K]
    alhc::Real = 2501                   # latent heat of condensation [J/g] for consistency with
                                        # specific humidity [g/Kg]
    alhs::Real = 2801                   # latent heat of sublimation [?]
    sbc::Real = 5.67e-8                 # stefan-Boltzmann constant [W/m^2/K^4]

    # STANDARD ATMOSPHERE (reference values)
    lapse_rate::Real = 6                # moist adiabatic temperature lapse rate -dT/dz [K/km]
    temp_ref::Real = 288                # absolute temperature at surface z=0 [K]
    temp_top::Real = 216                # absolute temperature in stratosphere [K]
    scale_height::Real = 7.5            # scale height for pressure [km]
    pres_ref::Real = 1013               # surface pressure [hPa]
    scale_height_humid::Real = 2.5      # scale height for specific humidity [km]
    relhumid_ref::Real = 0.7            # relative humidity of near-surface air [1]
    water_pres_ref::Real = 17           # saturation water vapour pressure [Pa]
    layer_thickness::Real = 8.5         # layer thickness for the shallow water model [km]

    # VERTICAL COORDINATES
    # of the nlev vertical levels, defined by a generalised logistic function,
    # interpolating ECMWF's L31 configuration
    GLcoefs::GenLogisticCoefs = GenLogisticCoefs()
    n_stratosphere_levels::Integer = 2                  # number of vertical levels used for the stratosphere
    σ_levels_half::Vector{Real} = []                    # vector of σ half levels, only if set manually, otherwise an empty vector
    nlev::Integer = nlev_default(model, σ_levels_half)  # number of vertical levels 

    # DIFFUSION AND DRAG
    diffusion_power::Real=4                 # Power n of Laplacian in horizontal diffusion ∇²ⁿ
    diffusion_time::Real=2.4                # Diffusion time scale [hrs] for temperature and vorticity
    diffusion_time_div::Real=diffusion_time # Diffusion time scale [hrs] for divergence
    diffusion_time_strat::Real=12           # Diffusion time scale [hrs] for extra ∇² in the stratosphere
    damping_time_strat::Real=24*30          # Damping time [hrs] for drag on zonal-mean wind in the stratosphere

    # FORCING
    interface_relaxation::Bool = false      # turn on interface relaxation for shallow water?
    interface_relax_time::Real = 96         # time scale [hrs] of interface relaxation
    interface_relax_amplitude::Real = 300   # Amplitude [m] of interface relaxation

    # PARAMETRIZATIONS
    n_shortwave::Integer = 3                # Compute shortwave radiation every n steps
    sppt_on::Bool = false                   # Turn on SPPT?
    magnus_coefs::MagnusCoefs = MagnusCoefs{NF}()  # For computing saturation vapour pressure

    # Large-Scale Condensation (from table B10)
    k_lsc::Int = 2                    # Index of atmospheric level at which large-scale condensation begins
    RH_thresh_pbl_lsc::Real = 0.95    # Relative humidity threshold for boundary layer
    RH_thresh_range_lsc::Real = 0.1   # Vertical range of relative humidity threshold
    RH_thresh_max_lsc::Real = 0.9     # Maximum relative humidity threshold
    humid_relax_time_lsc::Real = 4.0  # Relaxation time for humidity (hours)

    # Convection
    pres_thresh_cnv::Real = 0.8             # Minimum (normalised) surface pressure for the occurrence of convection
    RH_thresh_pbl_cnv::Real = 0.9           # Relative humidity threshold for convection in PBL
    RH_thresh_trop_cnv::Real = 0.7          # Relative humidity threshold for convection in the troposphere
    humid_relax_time_cnv::Real = 6.0        # Relaxation time for PBL humidity (hours)
    max_entrainment::Real = 0.5             # Maximum entrainment as a fraction of cloud-base mass flux
    ratio_secondary_mass_flux::Real = 0.8   # Ratio between secondary and primary mass flux at cloud-base

    # Radiation
    radiation_coefs::RadiationCoefs = RadiationCoefs{NF}()

    # TIME STEPPING
    startdate::DateTime = DateTime(2000,1,1)# time at which the integration starts
    n_days::Real = 10                       # number of days to integrate for
    Δt_at_T31::Real = 60                    # time step in minutes for T31, scale linearly to trunc

    # NUMERICS
    robert_filter::Real = 0.05          # Robert (1966) time filter coefficeint to suppress comput. mode
    williams_filter::Real = 0.53        # William's time filter (Amezcua 2011) coefficient for 3rd order acc
    implicit_α::Real = 0.5              # coefficient for semi-implicit computations to filter gravity waves
    
    # LEGENDRE TRANSFORM AND POLYNOMIALS
    recompute_legendre::Bool = false    # recomputation is slower but requires less memory
    legendre_NF::DataType = Float64     # which format to use to calculate the Legendre polynomials
    legendre_shortcut::Symbol = :linear # :linear, :quadratic, :cubic, :lincub_coslat, :linquad_coslat²

    # BOUNDARY FILES
    boundary_path::String = ""          # package location is default
    orography_path::String = boundary_path
    orography_file::String = "orography_F512.nc"

    # INITIAL CONDITIONS
    seed::Int = abs(rand(Int))          # random seed for the global random number generator
    initial_conditions::Symbol = :barotropic_vorticity    # :rest, :barotropic_vorticity or :restart

    # OUTPUT
    verbose::Bool = true            # print dialog for feedback
    output::Bool = false            # Store data in netCDF?
    output_dt::Real = 6             # output time step [hours]
    output_path::String = pwd()     # path to output folder

    # name of the output folder, defaults to 4-digit number counting up from run-0001
    run_id::Union{String,Integer} = get_run_id(output, output_path)    
    output_filename::String="output.nc"     # name of the output netcdf file
    
    # variables to output: :u, :v, :vor, :div, :temp, :humid
    output_vars::Vector{Symbol} = output_vars_default(model)
    compression_level::Integer = 3          # 1=low but fast, 9=high but slow
    keepbits::Integer = 7                   # mantissa bits to keep for every variable
    version::VersionNumber = pkgversion(SpeedyWeather)

    # OUTPUT GRID
    output_NF::DataType = NF        # number format used for output
    missing_value::Real = NaN       # missing value to be used in netcdf output
    output_grid::Symbol = :full     # :full, pick the corresponding full grid for reduced grids
                                    # or :matrix, sort gridpoints into a matrix
    output_quadrant_rotation::NTuple{4,Integer} = (0,1,2,3)
    output_matrix_quadrant::NTuple{4,Tuple{Integer,Integer}} = ((2,2),(1,2),(1,1),(2,1))

    # RESTART
    write_restart::Bool = output            # also write restart file if output==true?
    restart_path::String = output_path      # path for restart file
    restart_id::Union{String,Integer} = 1   # run_id of restart file in run-????/restart.jld2
end

"""
    nlev = nlev_default(model::Symbol, σ_levels_half::AbstractVector)

Number of vertical levels chosen either automatically based on `model`,
or from the length of `σ_levels_half` if not a 0-length vector
(default if not specified parameter).
"""
function nlev_default(model::Type{<:ModelSetup}, σ_levels_half::AbstractVector)
    if length(σ_levels_half) == 0   # choose nlev automatically 
        model <: Barotropic && return 1
        model <: ShallowWater && return 1
        model <: PrimitiveEquation && return 8
    else                            # use manually set levels 
        return length(σ_levels_half) - 1
    end
end

# default variables to output by model
output_vars_default(::Type{Barotropic}) = [:vor,:u]
output_vars_default(::Type{ShallowWater}) = [:vor,:u]
output_vars_default(::Type{PrimitiveEquation}) = [:vor,:u,:temp]
