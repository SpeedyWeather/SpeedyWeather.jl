const DEFAULT_NF = Float64          # number format
const DEFAULT_MODEL = Barotropic    # abstract model type

"""
    P = Parameters{M<:ModelSetup}(kwargs...) <: AbstractParameters{M}

A struct to hold all model parameters that may be changed by the user.
The struct uses keywords such that default values can be changed at creation.
The default values of the keywords define the default model setup.
"""
Base.@kwdef struct Parameters{Model<:ModelSetup} <: AbstractParameters{Model}

    # NUMBER FORMATS
    NF::DataType = DEFAULT_NF               # number format

    # RESOLUTION AND GRID
    trunc::Int = 31                                     # spectral truncation
    Grid::Type{<:AbstractGrid} = FullGaussianGrid       # grid used
    
    # EARTH'S PROPERTIES
    radius_earth::Float64 = 6.371e6         # radius of Earth [m]
    rotation_earth::Float64 = 7.29e-5       # angular frequency of Earth's rotation [rad/s]
    gravity::Float64 = 9.81                 # gravitational acceleration [m/s^2]
    seasonal_cycle::Bool = true             # Seasonal cycle?
    equinox::DateTime = DateTime(2000,3,20) # Spring equinox (year irrelevant)
    tropic_cancer::Float64 = 23.5           # latitude [˚N] of the tropic of cancer
    
    # ATMOSPHERE
    mol_mass_dry_air = 28.9649              # molar mass of dry air [g/mol]
    mol_mass_vapour = 18.0153               # molar mass of water vapour [g/mol]
    cₚ::Float64 = 1004                      # specific heat at constant pressure [J/K/kg]
    R_gas::Float64 = 8.3145                 # universal gas constant [J/K/mol]
    R_dry::Float64 = 1000*R_gas/mol_mass_dry_air   # specific gas constant for dry air [J/kg/K]
    R_vapour::Float64 = 1000*R_gas/mol_mass_vapour # specific gas constant for water vapour [J/kg/K]
    alhc::Float64 = 2501                    # latent heat of condensation [J/g] for consistency with
                                            # specific humidity [g/Kg]
    alhs::Float64 = 2801                    # latent heat of sublimation [?]
    sbc::Float64 = 5.67e-8                  # stefan-Boltzmann constant [W/m^2/K^4]

    # STANDARD ATMOSPHERE (reference values)
    lapse_rate::Float64 = 5                 # moist adiabatic temperature lapse rate -dT/dz [K/km]
    temp_ref::Float64 = 288                 # absolute temperature at surface z=0 [K]
    temp_top::Float64 = 216                 # absolute temperature in stratosphere [K]
    ΔT_stratosphere::Float64 = 4.8e5        # for stratospheric lapse rate [K] after Jablonowski
    scale_height::Float64 = 7.5             # scale height for pressure [km]
    pres_ref::Float64 = 1000                # surface pressure [hPa]
    scale_height_humid::Float64 = 2.5       # scale height for specific humidity [km]
    relhumid_ref::Float64 = 0.7             # relative humidity of near-surface air [1]
    water_pres_ref::Float64 = 17            # saturation water vapour pressure [Pa]
    layer_thickness::Float64 = 8.5          # layer thickness for the shallow water model [km]

    # VERTICAL COORDINATES
    # of the nlev vertical levels, defined by a generalised logistic function,
    # interpolating ECMWF's L31 configuration
    GLcoefs::Coefficients = GenLogisticCoefs()
    n_stratosphere_levels::Int = 2                  # number of vertical levels used for the stratosphere
    σ_tropopause::Float64 = 0.2                     # σ coordinate where the tropopause starts
    σ_levels_half::Vector{Float64} = []             # only used if set manually, otherwise empty
    nlev::Int = nlev_default(Model, σ_levels_half)  # number of vertical levels 

    # DIFFUSION AND DRAG
    diffusion_power::Float64=4                  # Power n of Laplacian in horizontal diffusion ∇²ⁿ
    diffusion_time::Float64=2.4                 # Diffusion time scale [hrs] for temperature and vorticity
    diffusion_time_div::Float64=diffusion_time  # Diffusion time scale [hrs] for divergence
    diffusion_time_strat::Float64=12            # Diffusion time scale [hrs] for extra ∇² in the stratosphere
    damping_time_strat::Float64=24*30           # Damping time [hrs] for drag on zonal-mean wind in the stratosphere

    # FORCING
    interface_relaxation::Bool = false          # turn on interface relaxation for shallow water?
    interface_relax_time::Float64 = 96          # time scale [hrs] of interface relaxation
    interface_relax_amplitude::Float64 = 300    # Amplitude [m] of interface relaxation

    # PARAMETRIZATIONS
    physics::Bool = true                        # en/disables the physics parameterizations

    n_shortwave::Int = 3                        # Compute shortwave radiation every n steps
    sppt_on::Bool = false                       # Turn on SPPT?
    magnus_coefs::Coefficients = MagnusCoefs{NF}()  # For computing saturation vapour pressure

    # Large-Scale Condensation (from table B10)
    k_lsc::Int = 2                    # Index of atmospheric level at which large-scale condensation begins
    RH_thresh_pbl_lsc::Float64 = 0.95    # Relative humidity threshold for boundary layer
    RH_thresh_range_lsc::Float64 = 0.1   # Vertical range of relative humidity threshold
    RH_thresh_max_lsc::Float64 = 0.9     # Maximum relative humidity threshold
    humid_relax_time_lsc::Float64 = 4.0  # Relaxation time for humidity (hours)

    # Convection
    pres_thresh_cnv::Float64 = 0.8              # Minimum (normalised) surface pressure for the occurrence of convection
    RH_thresh_pbl_cnv::Float64 = 0.9            # Relative humidity threshold for convection in PBL
    RH_thresh_trop_cnv::Float64 = 0.7           # Relative humidity threshold for convection in the troposphere
    humid_relax_time_cnv::Float64 = 6.0         # Relaxation time for PBL humidity (hours)
    max_entrainment::Float64 = 0.5              # Maximum entrainment as a fraction of cloud-base mass flux
    ratio_secondary_mass_flux::Float64 = 0.8    # Ratio between secondary and primary mass flux at cloud-base

    # Longwave radiation
    nband::Int = 4                          # Number of bands used to compute fband

    # Radiation
    radiation_coefs::Coefficients = RadiationCoefs{NF}()

    # BOUNDARY LAYER
    boundary_layer::BoundaryLayer = LinearDrag{NF}()
    temperature_relaxation::TemperatureRelaxation = HeldSuarez{NF}()

    # TIME STEPPING
    startdate::DateTime = DateTime(2000,1,1)    # time at which the integration starts
    n_days::Float64 = 10                        # number of days to integrate for
    Δt_at_T31::Float64 = 30                     # time step in minutes for T31, scale linearly to trunc

    # NUMERICS
    robert_filter::Float64 = 0.05       # Robert (1966) time filter coefficeint to suppress comput. mode
    williams_filter::Float64 = 0.53     # William's time filter (Amezcua 2011) coefficient for 3rd order acc
    implicit_α::Float64 = 1             # coefficient for semi-implicit computations to filter gravity waves
    
    # LEGENDRE TRANSFORM AND POLYNOMIALS
    recompute_legendre::Bool = false    # recomputation is slower but requires less memory
    legendre_NF::DataType = Float64     # which format to use to calculate the Legendre polynomials
    legendre_shortcut::Symbol = :linear # :linear, :quadratic, :cubic, :lincub_coslat, :linquad_coslat²

    # BOUNDARY FILES
    boundary_path::String = ""          # package location is default
    orography::AbstractOrography = EarthOrography()
    orography_scale::Float64 = 1        # scale orography by a factor
    orography_path::String = boundary_path
    orography_file::String = "orography_F512.nc"

    # INITIAL CONDITIONS
    seed::Int = 123456789           # random seed for the global random number generator
    initial_conditions::InitialConditions = initial_conditions_default(Model)
    pressure_on_orography::Bool = false # calculate the initial surface pressure from orography

    # OUTPUT
    verbose::Bool = true            # print dialog for feedback
    debug::Bool = true              # print debug info, NaR detection
    output::Bool = false            # Store data in netCDF?
    output_dt::Float64 = 6          # output time step [hours]
    output_path::String = pwd()     # path to output folder

    # name of the output folder, defaults to 4-digit number counting up from run-0001
    run_id::Union{String,Int} = get_run_id(output, output_path)    
    output_filename::String = "output.nc"   # name of the output netcdf file
    
    # variables to output: :u, :v, :vor, :div, :temp, :humid
    output_vars::Vector{Symbol} = output_vars_default(Model)
    compression_level::Int = 3          # 1=low but fast, 9=high but slow
    keepbits::Int = 7                   # mantissa bits to keep for every variable
    version::VersionNumber = pkgversion(SpeedyWeather)

    # OUTPUT GRID
    output_NF::DataType = Float32       # number format used for output
    output_nlat_half::Int = 0           # 0 = reuse nlat_half from dynamical core
    output_Grid::Type{<:AbstractFullGrid} = RingGrids.full_grid(Grid)
    output_Interpolator::Type{<:AbstractInterpolator} = DEFAULT_INTERPOLATOR
    output_matrix::Bool = false         # if true sort gridpoints into a matrix
    output_quadrant_rotation::NTuple{4,Int} = (0,1,2,3)
    output_matrix_quadrant::NTuple{4,Tuple{Int,Int}} = ((2,2),(1,2),(1,1),(2,1))
    
    missing_value::Float64 = NaN        # missing value to be used in netcdf output
    
    # RESTART
    write_restart::Bool = output        # also write restart file if output==true?
    restart_path::String = output_path  # path for restart file
    restart_id::Union{String,Int} = 1   # run_id of restart file in run-????/restart.jld2
end

Parameters(;kwargs...) = Parameters{default_concrete_model(DEFAULT_MODEL)}(;kwargs...)

"""
    nlev = nlev_default(Model::Type{<:ModelSetup}, σ_levels_half::AbstractVector)

Number of vertical levels chosen either automatically based on `Model`,
or from the length of `σ_levels_half` if not a 0-length vector
(default if not specified parameter).
"""
function nlev_default(Model::Type{<:ModelSetup}, σ_levels_half::AbstractVector)
    if length(σ_levels_half) == 0   # choose nlev automatically 
        Model <: Barotropic && return 1
        Model <: ShallowWater && return 1
        Model <: PrimitiveEquation && return 8
    else                            # use manually set levels 
        return length(σ_levels_half) - 1
    end
end

# default variables to output by model
output_vars_default(::Type{<:Barotropic}) = [:vor,:u]
output_vars_default(::Type{<:ShallowWater}) = [:vor,:u]
output_vars_default(::Type{<:PrimitiveDryCore}) = [:vor,:u,:temp,:pres]
output_vars_default(::Type{<:PrimitiveWetCore}) = [:vor,:u,:temp,:humid,:pres]

# default initial conditions by model
initial_conditions_default(::Type{<:Barotropic}) = StartWithVorticity()
initial_conditions_default(::Type{<:ShallowWater}) = ZonalJet()
initial_conditions_default(::Type{<:PrimitiveEquation}) = ZonalWind()