"""
$(TYPEDSIGNATURES)
Compute tendencies for u,v,temp,humid from physical parametrizations.
Extract for each vertical atmospheric column the prognostic variables
(stored in `diagn` as they are grid-point transformed), loop over all
grid-points, compute all parametrizations on a single-column basis,
then write the tendencies back into a horizontal field of tendencies.
"""
function parameterization_tendencies!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    time::DateTime,
    model::PrimitiveEquation,
)
    # TODO move into shortwave radiation code
    cos_zenith!(time,model)

    G = model.geometry
    rings = eachring(G.Grid,G.nlat_half)

    @floop for ij in eachgridpoint(diagn)       # loop over all horizontal grid points

        thread_id = Threads.threadid()          # not two threads should use the same ColumnVariables
        column = diagn.columns[thread_id]
        jring = whichring(ij,rings)             # ring index gridpoint ij is on

        # extract current column for contiguous memory access
        reset_column!(column)                   # set accumulators back to zero for next grid point
        get_column!(column,diagn,progn,ij,jring,model)  
        
        # execute all parameterizations
        parameterization_tendencies!(column,model)

        # write tendencies from parametrizations back into horizontal fields
        write_column_tendencies!(diagn,column,model.constants,ij)
    end
end

function parameterization_tendencies!(
    column::ColumnVariables,
    model::PrimitiveEquation
    )

    # Pre-compute thermodynamic quantities
    get_thermodynamics!(column,model)

    temperature_relaxation!(column,model)
    boundary_layer_drag!(column,model)

    # VERTICAL DIFFUSION
    # diffusion_coefficient!(column,model)
    # momentum_diffusion!(column,model)
    static_energy_diffusion!(column,model)
    humidity_diffusion!(column,model)

    # Calculate parametrizations
    convection!(column,model)
    large_scale_condensation!(column,model)
    # clouds!(column, model)
    shortwave_radiation!(column,model)
    longwave_radiation!(column,model)
    surface_fluxes!(column,model)

    # sum fluxes on half levels up and down for every layer
    fluxes_to_tendencies!(column,model.geometry,model.constants)
end

"""
$(TYPEDSIGNATURES)
Convert the fluxes on half levels to tendencies on full levels."""
function fluxes_to_tendencies!(
    column::ColumnVariables,
    geometry::Geometry,
    constants::DynamicsConstants,
)
    
    (;nlev, u_tend, flux_u_upward, flux_u_downward) = column
    (;      v_tend, flux_v_upward, flux_v_downward) = column
    (;humid_tend, flux_humid_upward, flux_humid_downward) = column
    (;temp_tend,  flux_temp_upward,  flux_temp_downward) = column

    Δσ = geometry.σ_levels_thick
    pₛ = column.pres[end]               # surface pressure
    (;radius) = constants               # used for scaling

    # for g/Δp and g/(Δp*cₚ), see Fortran SPEEDY documentation eq. (3,5)
    g_pₛ = constants.gravity/pₛ

    # fluxes are defined on half levels including top k=1/2 and surface k=nlev+1/2
    @inbounds for k in 1:nlev

        # Absorbed flux in a given layer, i.e. flux in minus flux out from above and below
        # Fortran SPEEDY documentation eq. (2)
        ΔF_u = (flux_u_upward[k+1] - flux_u_upward[k]) +
            (flux_u_downward[k] - flux_u_downward[k+1])
        
        ΔF_v = (flux_v_upward[k+1] - flux_v_upward[k]) +
            (flux_v_downward[k] - flux_v_downward[k+1])

        ΔF_humid = (flux_humid_upward[k+1] - flux_humid_upward[k]) +
            (flux_humid_downward[k] - flux_humid_downward[k+1])

        ΔF_temp = (flux_temp_upward[k+1] - flux_temp_upward[k]) +
            (flux_temp_downward[k] - flux_temp_downward[k+1])

        # convert absorbed flux to tendency, accumulate with
        # non-flux tendencies and scale with radius
        g_Δp = g_pₛ/Δσ[k]
        g_Δp_cₚ = g_Δp/constants.cₚ
        u_tend[k] = radius*(u_tend[k] + g_Δp*ΔF_u)
        v_tend[k] = radius*(v_tend[k] + g_Δp*ΔF_v)
        humid_tend[k] = radius*(humid_tend[k] + g_Δp*ΔF_humid)
        temp_tend[k] = radius*(temp_tend[k] + g_Δp_cₚ*ΔF_temp)
    end

    return nothing
end