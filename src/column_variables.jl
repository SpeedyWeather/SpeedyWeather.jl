"""
    column = ColumnVariables{NF<:AbstractFloat}

Mutable struct that contains all prognostic (copies thereof) and diagnostic variables in a single column
needed to evaluate the physical parametrizations. For now the struct is mutable as we will reuse the struct
to iterate over horizontal grid points. Every column vector has `nlev` entries, from [1] at the top to
[end] at the lowermost model level at the planetary boundary layer."""
@with_kw mutable struct ColumnVariables{NF<:AbstractFloat}

    # COORDINATES
    lat::NF = 0                         # latitude [˚N], needed for radiation?
    lon::NF = 0                         # longitude [˚E], needed for radiation?
    nlev::Int = 0                       # number of vertical levels

    # PROGNOSTIC VARIABLES
    u::Vector{NF} = zeros(NF,nlev)      # zonal velocity
    v::Vector{NF} = zeros(NF,nlev)      # meridional velocity
    temp::Vector{NF} = zeros(NF,nlev)   # temperature
    humid::Vector{NF} = zeros(NF,nlev)  # specific humidity

    # surface layer prognostic variables
    log_pres::NF = 0                    # logarithm of surface pressure [log(hPa)]
    pres::NF = 0                        # surface pressure [hPa]

    # TENDENCIES to accumulate the parametrizations into
    u_tend::Vector{NF} = zeros(NF,nlev)
    v_tend::Vector{NF} = zeros(NF,nlev)
    temp_tend::Vector{NF} = zeros(NF,nlev)
    humid_tend::Vector{NF} = zeros(NF,nlev)

    # DIAGNOSTIC VARIABLES
    geopot::Vector{NF} = zeros(NF,nlev)

    ## PARAMETERIZATIONS
    # Thermodynamics
    humid_half::Vector{NF} = zeros(NF,nlev)                    # Specific humidity interpolated to half-levels
    sat_humid::Vector{NF} = zeros(NF,nlev)                     # Saturation specific humidity
    sat_humid_half::Vector{NF} = zeros(NF,nlev)                # Saturation specific humidity interpolated to half-levels
    sat_vap_pres::Vector{NF} = zeros(NF,nlev)                  # Saturation vapour pressure
    dry_static_energy::Vector{NF} = zeros(NF,nlev)             # Dry static energy
    dry_static_energy_half::Vector{NF} = zeros(NF,nlev)        # Dry static energy interpolated to half-levels
    moist_static_energy::Vector{NF} = zeros(NF,nlev)           # Moist static energy
    sat_moist_static_energy::Vector{NF} = zeros(NF,nlev)       # Saturation moist static energy
    sat_moist_static_energy_half::Vector{NF} = zeros(NF,nlev)  # Saturation moist static energy interpolated to half-levels

    # Convection
    conditional_instability::Bool = false                    # Whether a conditional instability exists in this column
    activate_convection::Bool = false                        # Whether convection is activated in this column
    cloud_top::Int = nlev+1                                  # Top-of-convection layer
    excess_humidity::NF = 0                                  # Excess humidity due to convection
    cloud_base_mass_flux::NF = 0                             # Mass flux at the top of the PBL
    precip_cnv::NF = 0                                       # Precipitation due to convection
    net_flux_humid::Vector{NF} = zeros(NF,nlev)              # Fluxes of moisture in this column
    net_flux_dry_static_energy::Vector{NF} = zeros(NF,nlev)  # Fluxes of dry static energy in this column
    entrainment_profile::Vector{NF} = zeros(NF,nlev)         # Entrainment coefficients

    # Large-scale condensation
    precip_large_scale::NF = 0  # Precipitation due to large-scale condensation
end

# use Float64 if not provided
ColumnVariables(;kwargs...) = ColumnVariables{Float64}(;kwargs...)

"""
    get_column!(C,D,ij,G)

Update `C::ColumnVariables` by copying the prognostic variables from `D::DiagnosticVariables`
at gridpoint index `ij`. Provide `G::Geometry` for coordinate information."""
function get_column!(   C::ColumnVariables,
                        D::DiagnosticVariables,
                        ij::Int,            # grid point index
                        G::Geometry,
                        )

    @boundscheck C.nlev == D.nlev || throw(BoundsError)

    C.lat = G.lats[ij]      # pull latitude, longitude for gridpoint ij from Geometry
    C.lon = G.lons[ij]
    coslat⁻¹ = 1/cos(C.lat)

    # surface pressure (logarithm used in dynamics, convert back here)
    C.log_pres = D.surface.pres_grid[ij]
    C.pres = exp(C.log_pres)

    @inbounds for (k,layer) = enumerate(D.layers)
        C.u[k] = layer.grid_variables.U_grid[ij]*coslat⁻¹
        C.v[k] = layer.grid_variables.V_grid[ij]*coslat⁻¹
        C.temp[k] = layer.grid_variables.temp_grid[ij]
        C.humid[k] = layer.grid_variables.humid_grid[ij]
        C.geopot[k] = layer.grid_variables.geopot_grid[ij]
    end
end

"""
    write_column_tendencies!(D,C,ij,G)

Write the parametrization tendencies from `C::ColumnVariables` into the horizontal fields
of tendencies stored in `D::DiagnosticVariables` at gridpoint index `ij`."""
function write_column_tendencies!(  D::DiagnosticVariables,
                                    C::ColumnVariables,
                                    ij::Int,            # grid point index
                                    )

    @boundscheck C.nlev == D.nlev || throw(BoundsError)

    @inbounds for (k,layer) =  enumerate(D.layers)
        layer.tendencies.u_tend[ij] = C.u_tend[k]
        layer.tendencies.v_tend[ij] = C.v_tend[k]
        layer.tendencies.temp_grid_tend[ij] = C.temp_tend[k]
        layer.tendencies.humid_grid_tend[ij] = C.humid_tend[k]
    end

    D.surface.precip_large_scale[ij] = C.precip_large_scale
    D.surface.precip_convection[ij] = C.precip_convection

    return nothing
end

"""
    reset_column!(column::ColumnVariables)

Set the accumulators (tendencies but also vertical sums and similar) back to zero
for `column` to be reused at other grid points."""
function reset_column!(column::ColumnVariables{NF}) where NF

    fill!(column.u_tend,0)      # set tendencies to 0 for += accumulation
    fill!(column.v_tend,0)
    fill!(column.temp_tend,0)
    fill!(column.humid_tend,0)

    # Convection
    column.cloud_top = column.nlev+1
    column.conditional_instability = false
    column.activate_convection = false
    column.excess_humidity = zero(NF)
    column.precip_convection = zero(NF)
    column.cloud_base_mass_flux = zero(NF)
    fill!(column.net_flux_humid, 0)
    fill!(column.net_flux_dry_static_energy, 0)

    # Large-scale condensation
    column.precip_large_scale = zero(NF)
    return nothing
end

# iterator for convenience
eachlayer(column::ColumnVariables) = Base.OneTo(column.nlev)
