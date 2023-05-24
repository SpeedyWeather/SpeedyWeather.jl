"""
    column = ColumnVariables{NF<:AbstractFloat}

Mutable struct that contains all prognostic (copies thereof) and diagnostic variables in a single column
needed to evaluate the physical parametrizations. For now the struct is mutable as we will reuse the struct
to iterate over horizontal grid points. Every column vector has `nlev` entries, from [1] at the top to
[end] at the lowermost model level at the planetary boundary layer."""
Base.@kwdef mutable struct ColumnVariables{NF<:AbstractFloat} <: AbstractColumnVariables{NF}

    # DIMENSIONS
    const nlev::Int = 0                     # number of vertical levels
    const nband::Int = 0                    # number of radiation bands, needed for radiation
    const n_stratosphere_levels::Int = 0    # number of stratospheric levels, needed for radiation
    
    # COORDINATES
    jring::Int = 0                          # latitude ring the column is on
    lond::NF = 0                            # longitude
    latd::NF = 0                            # latitude, needed for shortwave radiation

    # PROGNOSTIC VARIABLES
    const u::Vector{NF} = zeros(NF,nlev)            # zonal velocity
    const v::Vector{NF} = zeros(NF,nlev)            # meridional velocity
    const temp::Vector{NF} = zeros(NF,nlev)         # temperature
    const humid::Vector{NF} = zeros(NF,nlev)        # specific humidity

    # (log) pressure per layer, surface is prognostic, last element here, but precompute other layers too
    const ln_pres::Vector{NF} = zeros(NF,nlev+1)    # logarithm of pressure [log(hPa)]
    const pres::Vector{NF} = zeros(NF,nlev+1)       # pressure [hPa]

    # TENDENCIES to accumulate the parametrizations into
    const u_tend::Vector{NF} = zeros(NF,nlev)                   # zonal velocity [m]
    const v_tend::Vector{NF} = zeros(NF,nlev)                   # meridional velocity [m]
    const temp_tend::Vector{NF} = zeros(NF,nlev)                # absolute temperature [K]
    const humid_tend::Vector{NF} = zeros(NF,nlev)               # specific humidity

    # DIAGNOSTIC VARIABLES
    const geopot::Vector{NF} = zeros(NF,nlev)                   # gepotential height [m]

    # FLUXES, arrays to be used for various parameterizations, on half levels incl top and bottom
    const flux_u_upward::Vector{NF} = zeros(NF,nlev+1)
    const flux_u_downward::Vector{NF} = zeros(NF,nlev+1)

    const flux_v_upward::Vector{NF} = zeros(NF,nlev+1)
    const flux_v_downward::Vector{NF} = zeros(NF,nlev+1)

    const flux_temp_upward::Vector{NF} = zeros(NF,nlev+1)
    const flux_temp_downward::Vector{NF} = zeros(NF,nlev+1)
    
    const flux_humid_upward::Vector{NF} = zeros(NF,nlev+1)
    const flux_humid_downward::Vector{NF} = zeros(NF,nlev+1)

    # THERMODYNAMICS
    const sat_humid::Vector{NF} = zeros(NF,nlev)                # Saturation specific humidity
    const sat_vap_pres::Vector{NF} = zeros(NF,nlev)             # Saturation vapour pressure
    const dry_static_energy::Vector{NF} = zeros(NF,nlev)        # Dry static energy
    const moist_static_energy::Vector{NF} = zeros(NF,nlev)      # Moist static energy
    
    # an interpolated to half levels
    humid_half::Vector{NF} = zeros(NF,nlev)                    # Specific humidity interpolated to half-levels
    sat_humid_half::Vector{NF} = zeros(NF,nlev)                # Saturation specific humidity interpolated to half-levels
    sat_moist_static_energy::Vector{NF} = zeros(NF,nlev)       # Saturation moist static energy
    dry_static_energy_half::Vector{NF} = zeros(NF,nlev)        # Dry static energy interpolated to half-levels
    sat_moist_static_energy_half::Vector{NF} = zeros(NF,nlev)  # Saturation moist static energy interpolated to half-levels

    # Convection
    conditional_instability::Bool = false                    # Whether a conditional instability exists in this column (condition 1)
    activate_convection::Bool = false                        # Whether convection should be activated in this column (condition 2)
    cloud_top::Int = nlev+1                                  # Top-of-convection layer
    excess_humidity::NF = 0                                  # Excess humidity due to convection
    cloud_base_mass_flux::NF = 0                             # Mass flux at the top of the PBL
    precip_convection::NF = 0                                # Precipitation due to convection
    net_flux_humid::Vector{NF} = zeros(NF,nlev)              # Net fluxes of moisture in this column
    net_flux_dry_static_energy::Vector{NF} = zeros(NF,nlev)  # Net fluxes of dry static energy in this column
    entrainment_profile::Vector{NF} = zeros(NF,nlev)         # Entrainment coefficients

    # Large-scale condensation
    precip_large_scale::NF = 0  # Precipitation due to large-scale condensation

    # Longwave radiation
    ## New vars in radlw_down!
    wvi::Matrix{NF} = fill(NF(NaN), nlev, 2)  # Weights for vertical interpolation
    tau2::Matrix{NF} = fill(NF(NaN), nlev, nband) # Transmissivity of atmospheric layers
    dfabs::Vector{NF} = fill(NF(NaN), nlev)   # Flux of sw rad. absorbed by each atm. layer
    fsfcd::NF = NaN                       # Downward-only flux of sw rad. at the surface
    st4a::Matrix{NF} = fill(NF(NaN), nlev, 2) # Blackbody emission from full and half atmospheric levels
    flux::Vector{NF} = fill(NF(NaN), nband)       # Radiative flux in different spectral bands

    ## New vars in compute_bbe!
    fsfcu::NF = NaN # surface blackbody emission (upward)
    ts::NF = NaN    # surface temperature

    ## New vars in radlw_up!
    fsfc::NF = NaN # Net (downw.) flux of sw rad. at the surface
    ftop::NF = NaN # Net (downw.) flux of sw rad. at the atm. top
    stratc::Vector{NF} = fill(NF(NaN), n_stratosphere_levels) # Stratospheric correction term 
    # Shortwave radiation: solar
    tyear::NF = NF(NaN) # time as fraction of year (0-1, 0 = 1jan.h00)
    csol::NF = NF(NaN)  # FIXME
    topsr::NF = NF(NaN) # FIXME
    # Shortwave radiation: solar_oz
    fsol::NF = NF(NaN)   # Flux of incoming solar radiation
    ozupp::NF = NF(NaN)  # Flux absorbed by ozone (upper stratos.)
    ozone::NF = NF(NaN)  # Flux absorbed by ozone (lower stratos.)
    zenit::NF = NF(NaN)  # Optical depth ratio (function of solar zenith angle)
    stratz::NF = NF(NaN) # Stratospheric correction for polar night
    # Shortwave radiation: radsw
    albsfc::NF = NF(NaN) # Combined surface albedo (land + sea)
    ssrd::NF = NF(NaN)   # Surface shortwave radiation (downward-only)
    ssr::NF = NF(NaN)    # Surface shortwave radiation (net downward)
    tsr::NF = NF(NaN)    # Top-of-atm. shortwave radiation (downward)
    tend_t_rsw::Vector{NF} = fill(NF(NaN), nlev) # Tempterature tendency
    norm_pres::NF = NF(NaN) # Normalized pressure (p/1000 hPa)
    # Shortwave radiation: cloud
    icltop::Int = typemax(Int) # Cloud top level (all clouds)
    cloudc::NF = NF(NaN)       # Total cloud cover (fraction)
    clstr::NF = NF(NaN)        # Stratiform cloud cover (fraction)
    qcloud::NF = NF(NaN)       # Equivalent specific humidity of clouds
    fmask::NF = NF(NaN)        # Fraction of land
    # Shortwave radiation: shortwave_radiation
    rel_hum::Vector{NF} = fill(NF(NaN), nlev) # Relative humidity
    grad_dry_static_energy::NF = NF(NaN)      # gradient of dry static energy
end

# use Float64 if not provided
ColumnVariables(;kwargs...) = ColumnVariables{Float64}(;kwargs...)