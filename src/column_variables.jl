"""
    column = ColumnVariables{NF<:AbstractFloat}

Mutable struct that contains all prognostic (copies thereof) and diagnostic variables in a single column
needed to evaluate the physical parametrizations. For now the struct is mutable as we will reuse the struct
to iterate over horizontal grid points. Every column vector has `nlev` entries, from [1] at the top to
[end] at the lowermost model level at the planetary boundary layer."""
@with_kw mutable struct ColumnVariables{NF <: AbstractFloat}

    # COORDINATES
    nlev::Int = 0                        # number of vertical levels
    lon::NF = NaN                        # longitude
    lat::NF = NaN                        # latitude, needed for shortwave radiation
    nband::Int = 0                       # number of radiation bands, needed for radiation
    n_stratosphere_levels::Int = 0       # number of stratospheric levels, needed for radiation

    # PROGNOSTIC VARIABLES
    u::Vector{NF} = zeros(NF, nlev)      # zonal velocity
    v::Vector{NF} = zeros(NF, nlev)      # meridional velocity
    temp::Vector{NF} = zeros(NF, nlev)   # temperature
    humid::Vector{NF} = zeros(NF, nlev)  # specific humidity

    # surface layer prognostic variables
    log_pres::NF = 0                    # logarithm of surface pressure [log(hPa)]
    pres::NF = 0                        # surface pressure [hPa]

    # TENDENCIES to accumulate the parametrizations into
    u_tend::Vector{NF} = zeros(NF, nlev)
    v_tend::Vector{NF} = zeros(NF, nlev)
    temp_tend::Vector{NF} = zeros(NF, nlev)
    humid_tend::Vector{NF} = zeros(NF, nlev)

    # DIAGNOSTIC VARIABLES
    geopot::Vector{NF} = zeros(NF, nlev)

    ## PARAMETERIZATIONS
    # Thermodynamics
    humid_half::Vector{NF} = zeros(NF, nlev)                    # Specific humidity interpolated to half-levels
    sat_humid::Vector{NF} = zeros(NF, nlev)                     # Saturation specific humidity
    sat_humid_half::Vector{NF} = zeros(NF, nlev)                # Saturation specific humidity interpolated to half-levels
    sat_vap_pres::Vector{NF} = zeros(NF, nlev)                  # Saturation vapour pressure
    dry_static_energy::Vector{NF} = zeros(NF, nlev)             # Dry static energy
    dry_static_energy_half::Vector{NF} = zeros(NF, nlev)        # Dry static energy interpolated to half-levels
    moist_static_energy::Vector{NF} = zeros(NF, nlev)           # Moist static energy
    sat_moist_static_energy::Vector{NF} = zeros(NF, nlev)       # Saturation moist static energy
    sat_moist_static_energy_half::Vector{NF} = zeros(NF, nlev)  # Saturation moist static energy interpolated to half-levels

    # Convection
    conditional_instability::Bool = false                    # Whether a conditional instability exists in this column (condition 1)
    activate_convection::Bool = false                        # Whether convection should be activated in this column (condition 2)
    cloud_top::Int = nlev + 1                                  # Top-of-convection layer
    excess_humidity::NF = 0                                  # Excess humidity due to convection
    cloud_base_mass_flux::NF = 0                             # Mass flux at the top of the PBL
    precip_convection::NF = 0                                # Precipitation due to convection
    net_flux_humid::Vector{NF} = zeros(NF, nlev)              # Net fluxes of moisture in this column
    net_flux_dry_static_energy::Vector{NF} = zeros(NF, nlev)  # Net fluxes of dry static energy in this column
    entrainment_profile::Vector{NF} = zeros(NF, nlev)         # Entrainment coefficients

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
ColumnVariables(; kwargs...) = ColumnVariables{Float64}(; kwargs...)

"""
    get_column!(C,D,ij,G)

Update `C::ColumnVariables` by copying the prognostic variables from `D::DiagnosticVariables`
at gridpoint index `ij`. Provide `G::Geometry` for coordinate information."""
function get_column!(C::ColumnVariables,
                     D::DiagnosticVariables,
                     ij::Int,            # grid point index
                     G::Geometry)
    @boundscheck C.nlev == D.nlev || throw(BoundsError)

    C.lat = G.lats[ij]      # pull latitude, longitude for gridpoint ij from Geometry
    C.lon = G.lons[ij]
    coslat⁻¹ = 1 / cos(C.lat)

    # surface pressure (logarithm used in dynamics, convert back here)
    C.log_pres = D.surface.pres_grid[ij]
    C.pres = exp(C.log_pres)

    @inbounds for (k, layer) in enumerate(D.layers)
        C.u[k] = layer.grid_variables.U_grid[ij] * coslat⁻¹
        C.v[k] = layer.grid_variables.V_grid[ij] * coslat⁻¹
        C.temp[k] = layer.grid_variables.temp_grid[ij]
        C.humid[k] = layer.grid_variables.humid_grid[ij]
        C.geopot[k] = layer.grid_variables.geopot_grid[ij]
    end
end

"""
    write_column_tendencies!(D,C,ij,G)

Write the parametrization tendencies from `C::ColumnVariables` into the horizontal fields
of tendencies stored in `D::DiagnosticVariables` at gridpoint index `ij`."""
function write_column_tendencies!(D::DiagnosticVariables,
                                  C::ColumnVariables,
                                  ij::Int)
    @boundscheck C.nlev == D.nlev || throw(BoundsError)

    @inbounds for (k, layer) in enumerate(D.layers)
        layer.tendencies.u_tend[ij] = C.u_tend[k]
        layer.tendencies.v_tend[ij] = C.v_tend[k]
        layer.tendencies.temp_tend_grid[ij] = C.temp_tend[k]
        layer.tendencies.humid_tend_grid[ij] = C.humid_tend[k]
    end

    D.surface.precip_large_scale[ij] = C.precip_large_scale
    D.surface.precip_convection[ij] = C.precip_convection

    return nothing
end

"""
    reset_column!(column::ColumnVariables)

Set the accumulators (tendencies but also vertical sums and similar) back to zero
for `column` to be reused at other grid points."""
function reset_column!(column::ColumnVariables{NF}) where {NF}
    fill!(column.u_tend, 0)      # set tendencies to 0 for += accumulation
    fill!(column.v_tend, 0)
    fill!(column.temp_tend, 0)
    fill!(column.humid_tend, 0)

    # Convection
    column.cloud_top = column.nlev + 1
    column.conditional_instability = false
    column.activate_convection = false
    column.precip_convection = zero(NF)
    fill!(column.net_flux_humid, 0)
    fill!(column.net_flux_dry_static_energy, 0)

    # Large-scale condensation
    column.precip_large_scale = zero(NF)
    return nothing
end

# iterator for convenience
eachlayer(column::ColumnVariables) = Base.OneTo(column.nlev)
