abstract type AbstractColumnVariables end
export ColumnVariables

"""
Mutable struct that contains all prognostic (copies thereof) and diagnostic variables in a single column
needed to evaluate the physical parametrizations. For now the struct is mutable as we will reuse the struct
to iterate over horizontal grid points. Most column vectors have `nlayers` entries, from [1] at the top to
[end] at the lowermost model level at the planetary boundary layer.
$(TYPEDFIELDS)"""
@kwdef mutable struct ColumnVariables{
    NF,
    VectorType,
    MatrixType,
    } <: AbstractColumnVariables

    # DIMENSIONS
    const nlayers::Int = 0                  # number of vertical levels
    const nbands_shortwave::Int = 0         # number of spectral bands for shortwave radiation
    const nbands_longwave::Int = 0          # number of spectral bands for longwave radiation

    # COORDINATES
    ij::Int = 0                             # grid-point index
    jring::Int = 0                          # latitude ring the column is on
    lond::NF = 0                            # longitude
    latd::NF = 0                            # latitude, needed for shortwave radiation
    land_fraction::NF = 0                   # fraction of the column that is over land
    orography::NF = 0                       # orography height [m]

    # PROGNOSTIC VARIABLES
    const u::VectorType = zeros(NF, nlayers)            # zonal velocity [m/s]
    const v::VectorType = zeros(NF, nlayers)            # meridional velocity [m/s]
    const temp::VectorType = zeros(NF, nlayers)         # absolute temperature [K]
    const humid::VectorType = zeros(NF, nlayers)        # specific humidity [kg/kg]

    # (log) pressure per layer, surface is prognostic, last element here, but precompute other layers too
    const ln_pres::VectorType = zeros(NF, nlayers+1)    # logarithm of pressure [log(Pa)]
    const pres::VectorType = zeros(NF, nlayers+1)       # pressure [Pa]

    # TENDENCIES to accumulate the parametrizations into
    const u_tend::VectorType = zeros(NF, nlayers)       # zonal velocity [m/s²]
    const v_tend::VectorType = zeros(NF, nlayers)       # meridional velocity [m/s²]
    const temp_tend::VectorType = zeros(NF, nlayers)    # absolute temperature [K/s]
    const humid_tend::VectorType = zeros(NF, nlayers)   # specific humidity [kg/kg/s]

    # FLUXES, arrays to be used for various parameterizations, on half levels incl top and bottom
    const flux_u_upward::VectorType = zeros(NF, nlayers+1)
    const flux_u_downward::VectorType = zeros(NF, nlayers+1)

    const flux_v_upward::VectorType = zeros(NF, nlayers+1)
    const flux_v_downward::VectorType = zeros(NF, nlayers+1)

    const flux_temp_upward::VectorType = zeros(NF, nlayers+1)
    const flux_temp_downward::VectorType = zeros(NF, nlayers+1)
    
    const flux_humid_upward::VectorType = zeros(NF, nlayers+1)
    const flux_humid_downward::VectorType = zeros(NF, nlayers+1)

    # random value (scalar) from random pattern controlled by model.random_process
    random_value::NF = 0

    # boundary layer
    boundary_layer_depth::Int = 0
    boundary_layer_drag::NF = 0
    surface_geopotential::NF = 0
    surface_u::NF = 0
    surface_v::NF = 0
    surface_temp::NF = 0
    surface_humid::NF = 0
    surface_wind_speed::NF = 0

    # land
    skin_temperature_sea::NF = 0
    skin_temperature_land::NF = 0
    soil_moisture_availability::NF = 0

    # surface fluxes
    evaporative_flux::NF = 0            # land-sea mask fraction-weighted flux 
    evaporative_flux_ocean::NF = 0      # flux from ocean only
    evaporative_flux_land::NF = 0       # and from land

    sensible_heat_flux::NF = 0          # land-sea mask fraction-weighted flux 
    sensible_heat_flux_ocean::NF = 0    # flux from ocean only
    sensible_heat_flux_land::NF = 0     # and from land

    # THERMODYNAMICS
    surface_air_density::NF = 0
    const sat_humid::VectorType = zeros(NF, nlayers)                # Saturation specific humidity [kg/kg]
    const dry_static_energy::VectorType = zeros(NF, nlayers)        # Dry static energy [J/kg]
    const temp_virt::VectorType = zeros(NF, nlayers)                # virtual temperature [K]
    const geopot::VectorType = zeros(NF, nlayers)                   # gepotential height [m]

    # CONVECTION AND PRECIPITATION
    cloud_top::Int = nlayers+1              # layer index k of top-most layer with clouds
    precip_convection::NF = 0               # Precipitation due to convection [m]
    precip_large_scale::NF = 0              # precipitation due to large-scale condensation [m]
    precip_rate_convection::NF = 0          # Precipitation rate due to convection [m/s]
    precip_rate_large_scale::NF = 0         # precipitation rate due to large-scale condensation [m/s]

    # RADIATION
    cos_zenith::NF = 0                      # cosine of solar zenith angle [1]
    albedo_ocean::NF = 0                    # surface albedo over ocean [1]
    albedo_land::NF = 0                     # surface albedo over land [1]    

    # surface fluxes, positive down
    surface_shortwave_down::NF = 0          # surface shortwave radiation down (into land/sea)
    surface_shortwave_down_ocean::NF = 0    # ocean only
    surface_shortwave_down_land::NF = 0     # land only
    
    surface_shortwave_up::NF = 0            # surface shortwave radiation up (reflected)
    surface_shortwave_up_ocean::NF = 0
    surface_shortwave_up_land::NF = 0
    
    surface_longwave_down::NF = 0           # surface longwave radiation down (into land/sea)
    surface_longwave_down_ocean::NF = 0
    surface_longwave_down_land::NF = 0
    
    # surface longwave radiation up (into atmosphere)
    surface_longwave_up::NF = 0             # land-sea mask weighted flux   
    surface_longwave_up_ocean::NF = 0       # ocean only
    surface_longwave_up_land::NF = 0        # land only

    # top-of-atmosphere fluxes, positive up (outgoing)
    outgoing_longwave_radiation::NF = 0     # OLR [W/m^2]
    outgoing_shortwave_radiation::NF = 0    # same for shortwave reflection [W/m^2]

    # optical depth of the atmosphere, on half levels, for shortwave and longwave radiation
    const optical_depth_shortwave::MatrixType = zeros(NF, nlayers, nbands_shortwave)
    const optical_depth_longwave::MatrixType = zeros(NF, nlayers, nbands_longwave)

    # WORK ARRAYS
    const a::VectorType = zeros(NF, nlayers)
    const b::VectorType = zeros(NF, nlayers)
    const c::VectorType = zeros(NF, nlayers)
    const d::VectorType = zeros(NF, nlayers)
end

# generator based on spectral grid
ColumnVariables(SG::SpectralGrid; kwargs...) = ColumnVariables{SG.NF, SG.VectorType, SG.MatrixType}(; nlayers=SG.nlayers, kwargs...)

# generator assuming Julia Arrays
ColumnVariables{NF}(; kwargs...) where NF = ColumnVariables{NF, Vector{NF}, Matrix{NF}}(; kwargs...)