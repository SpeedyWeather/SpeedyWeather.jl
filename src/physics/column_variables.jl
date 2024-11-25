# function barrier
function get_column!(   
    C::ColumnVariables,
    D::DiagnosticVariables,
    P::PrognosticVariables,
    ij::Integer,        # grid point index
    jring::Integer,     # ring index 1 around North Pole to J around South Pole
    model::PrimitiveEquation,
)
    get_column!(C, D, P, ij, jring, 
        model.geometry,
        model.planet,
        model.orography,
        model.land_sea_mask,
        model.albedo,
        model.implicit)
end

"""
$(TYPEDSIGNATURES)
Update `C::ColumnVariables` by copying the prognostic variables from `D::DiagnosticVariables`
at gridpoint index `ij`. Provide `G::Geometry` for coordinate information."""
function get_column!(   
    C::ColumnVariables,
    D::DiagnosticVariables,
    P::PrognosticVariables,
    ij::Integer,        # grid point index
    jring::Integer,     # ring index 1 around North Pole to J around South Pole
    geometry::Geometry,
    planet::AbstractPlanet,
    orography::AbstractOrography,
    land_sea_mask::AbstractLandSeaMask,
    albedo::AbstractAlbedo,
    implicit::AbstractImplicit,
)

    (; σ_levels_full, ln_σ_levels_full) = geometry
    (; temp_profile) = implicit     # reference temperature

    @boundscheck C.nlayers == D.nlayers || throw(BoundsError)

    C.latd = geometry.latds[ij]     # pull latitude, longitude [˚N, ˚E] for gridpoint ij from Geometry
    C.lond = geometry.londs[ij]
    C.ij = ij                       # grid-point index
    C.jring = jring                 # ring index j of column, used to index latitude vectors
    C.land_fraction = land_sea_mask.mask[ij]
    C.orography = orography.orography[ij]
    C.surface_geopotential = C.orography * planet.gravity

    # pressure [Pa] or [log(Pa)]
    lnpₛ = D.grid.pres_grid[ij]             # logarithm of surf pressure used in dynamics
    pₛ = exp(lnpₛ)                          # convert back here
    C.ln_pres .= ln_σ_levels_full .+ lnpₛ   # log pressure on every level ln(p) = ln(σ) + ln(pₛ)
    C.pres[1:end-1] .= σ_levels_full.*pₛ    # pressure on every level p = σ*pₛ
    C.pres[end] = pₛ                        # last element is surface pressure pₛ

    (; u_grid_prev, v_grid_prev, temp_grid_prev, humid_grid_prev) = D.grid

    @inbounds for k in eachgrid(u_grid_prev, v_grid_prev, temp_grid_prev, humid_grid_prev)
        # read out prognostic variables on grid at previous time step
        # for numerical stability
        C.u[k] = u_grid_prev[ij, k]
        C.v[k] = v_grid_prev[ij, k]
        C.temp[k] = temp_grid_prev[ij, k] + temp_profile[k]
        C.humid[k] = humid_grid_prev[ij, k] 
    end

    # extract value from random pattern for this column
    C.random_value = D.grid.random_pattern[ij]

    # TODO skin = surface approximation for now
    C.skin_temperature_sea = P.ocean.sea_surface_temperature[ij]
    C.skin_temperature_land = P.land.land_surface_temperature[ij]
    C.soil_moisture_availability = D.physics.soil_moisture_availability[ij]

    # RADIATION
    C.cos_zenith = D.physics.cos_zenith[ij]
    C.albedo = albedo.albedo[ij]
end

"""Recalculate ring index if not provided."""
function get_column!(
    C::ColumnVariables,
    D::DiagnosticVariables,
    P::PrognosticVariables,
    ij::Int,            # grid point index
    model::PrimitiveEquation
)
    rings = eachring(D.grid.vor_grid)
    jring = whichring(ij, rings)
    get_column!(C, D, P, ij, jring, model)
end

function get_column(    S::AbstractSimulation,
                        ij::Integer,
                        verbose::Bool = true)
    (; prognostic_variables, diagnostic_variables, model) = S

    column = deepcopy(S.diagnostic_variables.columns[1])
    reset_column!(column)

    get_column!(column,
                diagnostic_variables,
                prognostic_variables,
                ij,
                model)

    # execute all parameterizations for this column to return a consistent state
    parameterization_tendencies!(column, S.model)

    verbose && @info "Receiving column at $(column.latd)˚N, $(column.lond)˚E."
    return column
end

"""
$(TYPEDSIGNATURES)
Write the parametrization tendencies from `C::ColumnVariables` into the horizontal fields
of tendencies stored in `D::DiagnosticVariables` at gridpoint index `ij`."""
function write_column_tendencies!(
    diagn::DiagnosticVariables,
    column::ColumnVariables,
    planet::AbstractPlanet,
    ij::Integer,                                # grid point index
)
    (; nlayers) = column
    @boundscheck nlayers == diagn.nlayers || throw(BoundsError)

    @inbounds for k in eachgrid(diagn.grid.vor_grid)
        diagn.tendencies.u_tend_grid[ij, k] = column.u_tend[k]
        diagn.tendencies.v_tend_grid[ij, k] = column.v_tend[k]
        diagn.tendencies.temp_tend_grid[ij, k] = column.temp_tend[k]
        diagn.tendencies.humid_tend_grid[ij, k] = column.humid_tend[k]
    end

    # accumulate precipitation [m]
    diagn.physics.precip_large_scale[ij] += column.precip_large_scale
    diagn.physics.precip_convection[ij] += column.precip_convection

    # Cloud top in height [m] from geopotential height divided by gravity, 0 for no clouds
    diagn.physics.cloud_top[ij] = column.cloud_top == nlayers+1 ? 0 : column.geopot[column.cloud_top]
    diagn.physics.cloud_top[ij] /= planet.gravity
    
    # just use layer index 1 (top) to nlayers (surface) for analysis, but 0 for no clouds
    # diagn.physics.cloud_top[ij] = column.cloud_top == nlayers+1 ? 0 : column.cloud_top

    # net surface fluxes of humidity and temperature, defined as positive downward (i.e. into ocean/land)
    diagn.physics.surface_flux_heat[ij] = column.flux_temp_downward[nlayers+1] -
                                            column.flux_temp_upward[nlayers+1]
    diagn.physics.surface_flux_humid[ij] = column.flux_humid_downward[nlayers+1] -
                                            column.flux_humid_upward[nlayers+1]

    return nothing
end

"""
$(TYPEDSIGNATURES)
Set the accumulators (tendencies but also vertical sums and similar) back to zero
for `column` to be reused at other grid points."""
function reset_column!(column::ColumnVariables{NF}) where NF

    # set tendencies to 0 for += accumulation
    column.u_tend .= 0
    column.v_tend .= 0
    column.temp_tend .= 0
    column.humid_tend .= 0

    # set fluxes to 0 for += accumulation
    column.flux_u_upward .= 0
    column.flux_u_downward .= 0
    column.flux_v_upward .= 0
    column.flux_v_downward .= 0
    column.flux_humid_upward .= 0
    column.flux_humid_downward .= 0
    column.flux_temp_upward .= 0
    column.flux_temp_downward .= 0

    # Convection and precipitation
    column.cloud_top = column.nlayers+1
    column.precip_convection = 0
    column.precip_large_scale = 0

    return nothing
end

# iterator for convenience
eachlayer(column::ColumnVariables) = eachindex(column)
Base.eachindex(column::ColumnVariables) = Base.OneTo(column.nlayers)