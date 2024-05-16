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
    (; temp_profile) = implicit     # reference temperature on this layer

    @boundscheck C.nlev == D.nlev || throw(BoundsError)

    C.latd = geometry.latds[ij]     # pull latitude, longitude [˚N, ˚E] for gridpoint ij from Geometry
    C.lond = geometry.londs[ij]
    C.ij = ij                       # grid-point index
    C.jring = jring                 # ring index j of column, used to index latitude vectors
    C.land_fraction = land_sea_mask.mask[ij]
    C.orography = orography.orography[ij]
    C.surface_geopotential = C.orography * planet.gravity

    # pressure [Pa] or [log(Pa)]
    lnpₛ = D.surface.pres_grid[ij]          # logarithm of surf pressure used in dynamics
    pₛ = exp(lnpₛ)                          # convert back here
    C.ln_pres .= ln_σ_levels_full .+ lnpₛ   # log pressure on every level ln(p) = ln(σ) + ln(pₛ)
    C.pres[1:end-1] .= σ_levels_full.*pₛ    # pressure on every level p = σ*pₛ
    C.pres[end] = pₛ                        # last element is surface pressure pₛ

    @inbounds for (k, layer) = enumerate(D.layers)   # read out prognostic variables on grid
        C.u[k] = layer.grid_variables.u_grid[ij]
        C.v[k] = layer.grid_variables.v_grid[ij]
        C.temp[k] = layer.grid_variables.temp_grid[ij]
        C.temp_virt[k] = layer.grid_variables.temp_virt_grid[ij]    # actually diagnostic
        C.humid[k] = layer.grid_variables.humid_grid[ij] 

        # and at previous time step, add temp reference profile back in as temp_grid_prev is anomaly
        C.u_prev[k] = layer.grid_variables.u_grid_prev[ij]
        C.v_prev[k] = layer.grid_variables.v_grid_prev[ij]
        C.temp_prev[k] = layer.grid_variables.temp_grid_prev[ij] + temp_profile[k]
        C.humid_prev[k] = layer.grid_variables.humid_grid_prev[ij] 
    end

    # TODO skin = surface approximation for now
    C.skin_temperature_sea = P.ocean.sea_surface_temperature[ij]
    C.skin_temperature_land = P.land.land_surface_temperature[ij]
    C.soil_moisture_availability = D.surface.soil_moisture_availability[ij]

    # RADIATION
    C.cos_zenith = D.surface.cos_zenith[ij]
    C.albedo = albedo.albedo[ij]
end

"""Recalculate ring index if not provided."""
function get_column!(   C::ColumnVariables,
                        D::DiagnosticVariables,
                        P::PrognosticVariables,
                        ij::Int,            # grid point index
                        model::PrimitiveEquation)

    SG = model.spectral_grid
    rings = eachring(SG.Grid, SG.nlat_half)
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
function write_column_tendencies!(  D::DiagnosticVariables,
                                    column::ColumnVariables,
                                    planet::AbstractPlanet,
                                    ij::Int)            # grid point index

    (; nlev) = column
    @boundscheck nlev == D.nlev || throw(BoundsError)

    @inbounds for (k, layer) = enumerate(D.layers)
        layer.tendencies.u_tend_grid[ij] = column.u_tend[k]
        layer.tendencies.v_tend_grid[ij] = column.v_tend[k]
        layer.tendencies.temp_tend_grid[ij] = column.temp_tend[k]
        layer.tendencies.humid_tend_grid[ij] = column.humid_tend[k]
    end

    # accumulate (set back to zero when netcdf output)
    D.surface.precip_large_scale[ij] += column.precip_large_scale
    D.surface.precip_convection[ij] += column.precip_convection

    # Output cloud top in height [m] from geopotential height divided by gravity,
    # but NaN for no clouds
    D.surface.cloud_top[ij] = column.cloud_top == nlev+1 ? 0 : column.geopot[column.cloud_top]
    D.surface.cloud_top[ij] /= planet.gravity
    
    # just use layer index 1 (top) to nlev (surface) for analysis, but 0 for no clouds
    # D.surface.cloud_top[ij] = column.cloud_top == nlev+1 ? 0 : column.cloud_top
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
    column.cloud_top = column.nlev+1            # also diagnostic from condensation
    column.precip_convection = 0
    column.precip_large_scale = 0

    return nothing
end

# iterator for convenience
eachlayer(column::ColumnVariables) = Base.OneTo(column.nlev)
Base.eachindex(column::ColumnVariables) = Base.OneTo(column.nlev)