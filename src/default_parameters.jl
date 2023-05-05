using DocStringExtensions

const DEFAULT_NF = Float64          # number format
const DEFAULT_MODEL = Barotropic    # abstract model type

"""
    P = Parameters{M<:ModelSetup}(kwargs...) <: AbstractParameters{M}

A struct to hold all model parameters that may be changed by the user.
The struct uses keywords such that default values can be changed at creation.
The default values of the keywords define the default model setup.

$(TYPEDFIELDS)
"""
Base.@kwdef struct Parameters{Model<:ModelSetup} <: AbstractParameters{Model}

    "number format"
    NF::DataType = DEFAULT_NF


    # RESOLUTION AND GRID

    "spectral truncation"
    trunc::Int = 31

    "grid in use"
    Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "dealiasing factor, 1=linear, 2=quadratic, 3=cubic grid"
    dealiasing::Float64 = 2


    # PLANET'S PROPERTIES

    "planet"
    planet::Planet = Earth()


    # ATMOSPHERE

    "molar mass of dry air [g/mol]"
    mol_mass_dry_air = 28.9649

    "molar mass of water vapour [g/mol]"
    mol_mass_vapour = 18.0153

    "specific heat at constant pressure [J/K/kg]"
    cₚ::Float64 = 1004

    "universal gas constant [J/K/mol]"
    R_gas::Float64 = 8.3145

    "specific gas constant for dry air [J/kg/K]"
    R_dry::Float64 = 1000*R_gas/mol_mass_dry_air

    "specific gas constant for water vapour [J/kg/K]"
    R_vapour::Float64 = 1000*R_gas/mol_mass_vapour

    "latent heat of condensation [J/g] for consistency with specific humidity [g/Kg]"
    alhc::Float64 = 2501

    "latent heat of sublimation [?]"
    alhs::Float64 = 2801

    "stefan-Boltzmann constant [W/m²/K⁴]"
    sbc::Float64 = 5.67e-8


    # STANDARD ATMOSPHERE (reference values)

    "moist adiabatic temperature lapse rate ``-dT/dz`` [K/km]"
    lapse_rate::Float64 = 5

    "absolute temperature at surface ``z=0`` [K]"
    temp_ref::Float64 = 288

    "absolute temperature in stratosphere [K]"
    temp_top::Float64 = 216

    "for stratospheric lapse rate [K] after Jablonowski"
    ΔT_stratosphere::Float64 = 4.8e5

    "scale height for pressure [km]"
    scale_height::Float64 = 7.5

    "surface pressure [hPa]"
    pres_ref::Float64 = 1000

    "scale height for specific humidity [km]"
    scale_height_humid::Float64 = 2.5

    "relative humidity of near-surface air [1]"
    relhumid_ref::Float64 = 0.7

    "saturation water vapour pressure [Pa]"
    water_pres_ref::Float64 = 17

    "layer thickness for the shallow water model [km]"
    layer_thickness::Float64 = 8.5


    # VERTICAL COORDINATES

    "vertical coordinates of the nlev vertical levels, defined by a generalised logistic function, interpolating ECMWF's L31 configuration"
    GLcoefs::Coefficients = GenLogisticCoefs()

    "number of vertical levels used for the stratosphere"
    n_stratosphere_levels::Int = 2

    "σ coordinate where the tropopause starts"
    σ_tropopause::Float64 = 0.2

    "only used if set manually, otherwise empty"
    σ_levels_half::Vector{Float64} = []

    "number of vertical levels "
    nlev::Int = nlev_default(Model, σ_levels_half)


    # DIFFUSION AND DRAG

    "horizontal (hyper)-diffusion"
    diffusion::DiffusionParameters = HyperDiffusion()

    "vertical diffusion"
    vertical_diffusion::VerticalDiffusion = NoVerticalDiffusion()


    # FORCING

    "turn on interface relaxation for shallow water?"
    interface_relaxation::Bool = false

    "time scale [hrs] of interface relaxation"
    interface_relax_time::Float64 = 96

    "Amplitude [m] of interface relaxation"
    interface_relax_amplitude::Float64 = 300


    # PARAMETRIZATIONS

    "en/disables the physics parameterizations"
    physics::Bool = true


    "Compute shortwave radiation every `n` steps"
    n_shortwave::Int = 3

    "Turn on SPPT?"
    sppt_on::Bool = false

    "For computing saturation vapour pressure"
    magnus_coefs::Coefficients = MagnusCoefs{NF}()

    # Large-Scale Condensation (from table B10)
    "Index of atmospheric level at which large-scale condensation begins"
    k_lsc::Int = 2

    "Relative humidity threshold for boundary layer"
    RH_thresh_pbl_lsc::Float64 = 0.95

    "Vertical range of relative humidity threshold"
    RH_thresh_range_lsc::Float64 = 0.1

    "Maximum relative humidity threshold"
    RH_thresh_max_lsc::Float64 = 0.9

    "Relaxation time for humidity (hours)"
    humid_relax_time_lsc::Float64 = 4.0

    # Convection
    "Minimum (normalised) surface pressure for the occurrence of convection"
    pres_thresh_cnv::Float64 = 0.8

    "Relative humidity threshold for convection in PBL"
    RH_thresh_pbl_cnv::Float64 = 0.9

    "Relative humidity threshold for convection in the troposphere"
    RH_thresh_trop_cnv::Float64 = 0.7

    "Relaxation time for PBL humidity (hours)"
    humid_relax_time_cnv::Float64 = 6.0

    "Maximum entrainment as a fraction of cloud-base mass flux"
    max_entrainment::Float64 = 0.5

    "Ratio between secondary and primary mass flux at cloud-base"
    ratio_secondary_mass_flux::Float64 = 0.8

    # Longwave radiation
    "Number of bands used to compute `fband`"
    nband::Int = 4

    # Radiation
    "radiation coefficients"
    radiation_coefs::Coefficients = RadiationCoefs{NF}()

    # BOUNDARY LAYER
    "boundary layer drag"
    boundary_layer::BoundaryLayer{Float64} = LinearDrag()

    "temperature relaxation"
    temperature_relaxation::TemperatureRelaxation{Float64} = HeldSuarez()


    # TIME STEPPING

    "time at which the integration starts"
    startdate::DateTime = DateTime(2000,1,1)

    "number of days to integrate for"
    n_days::Float64 = 10

    "time step in minutes for T31, scale linearly to `trunc`"
    Δt_at_T31::Float64 = 30


    # NUMERICS

    "Robert (1966) time filter coefficeint to suppress comput. mode"
    robert_filter::Float64 = 0.05

    "William's time filter (Amezcua 2011) coefficient for 3rd order acc"
    williams_filter::Float64 = 0.53

    "coefficient for semi-implicit computations to filter gravity waves"
    implicit_α::Float64 = 1

    "recalculate implicit operators on temperature profile every n time steps"
    recalculate_implicit::Int = 100

    # LEGENDRE TRANSFORM AND POLYNOMIALS

    "recomputation is slower but requires less memory"
    recompute_legendre::Bool = false

    "which format to use to calculate the Legendre polynomials"
    legendre_NF::DataType = Float64

    "`:linear`, `:quadratic`, `:cubic`, `:lincub_coslat`, `:linquad_coslat²`"
    legendre_shortcut::Symbol = :linear


    # BOUNDARY FILES

    "package location is default"
    boundary_path::String = ""

    "orography"
    orography::AbstractOrography = EarthOrography(smoothing_strength=1e-2,smoothing_power=0.1)

    "scale orography by a factor"
    orography_scale::Float64 = 1

    "path of orography"
    orography_path::String = boundary_path

    "filename of orography"
    orography_file::String = "orography_F512.nc"


    # INITIAL CONDITIONS

    "random seed for the global random number generator"
    seed::Int = 123456789

    "initial conditions"
    initial_conditions::InitialConditions = initial_conditions_default(Model)

    "calculate the initial surface pressure from orography"
    pressure_on_orography::Bool = false


    # OUTPUT

    "print dialog for feedback"
    verbose::Bool = true

    "print debug info, NaR detection"
    debug::Bool = true

    "Store data in netCDF?"
    output::Bool = false

    "output time step [hours]"
    output_dt::Float64 = 6

    "path to output folder"
    output_path::String = pwd()

    "name of the output folder, defaults to 4-digit number counting up from `run-0001`"
    run_id::Union{String,Int} = get_run_id(output, output_path)

    "name of the output netcdf file"
    output_filename::String = "output.nc"
    
    "variables to output: `:u`, `:v`, `:vor`, `:div`, `:temp`, `:humid`"
    output_vars::Vector{Symbol} = output_vars_default(Model)

    "compression level; 1=low but fast, 9=high but slow"
    compression_level::Int = 3

    "mantissa bits to keep for every variable"
    keepbits::Int = 7

    "SpeedyWeather.jl version number"
    version::VersionNumber = pkgversion(SpeedyWeather)


    # OUTPUT GRID

    "number format used for output"
    output_NF::DataType = Float32

    "0 = reuse nlat_half from dynamical core"
    output_nlat_half::Int = 0

    "output grid"
    output_Grid::Type{<:AbstractFullGrid} = RingGrids.full_grid(Grid)

    "output interpolator"
    output_Interpolator::Type{<:AbstractInterpolator} = DEFAULT_INTERPOLATOR

    "if true sort gridpoints into a matrix"
    output_matrix::Bool = false

    "rotation of output quadrant"
    output_quadrant_rotation::NTuple{4,Int} = (0,1,2,3)

    "matrix of output quadrant"
    output_matrix_quadrant::NTuple{4,Tuple{Int,Int}} = ((2,2),(1,2),(1,1),(2,1))
    
    "missing value to be used in netcdf output"
    missing_value::Float64 = NaN


    # RESTART

    "also write restart file if output==true?"
    write_restart::Bool = output

    "path for restart file"
    restart_path::String = output_path

    "`run_id` of restart file in `run-????/restart.jld2`"
    restart_id::Union{String,Int} = 1
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