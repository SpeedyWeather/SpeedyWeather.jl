"""$(TYPEDSIGNATURES)
Compute tendencies for u, v, temp, humid from physical parametrizations.
Extract for each vertical atmospheric column the prognostic variables
(stored in `diagn` as they are grid-point transformed), loop over all
grid-points, compute all parametrizations on a single-column basis,
then write the tendencies back into a horizontal field of tendencies.
"""
function parameterization_tendencies!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    model::PrimitiveEquation,
)
    # parameterizations with their own kernel
    (; time) = progn.clock
    cos_zenith!(diagn, time, model)

    # all other parameterizations are fused into a single kernel over horizontal grid point index ij
    (; architecture, npoints) = model.spectral_grid
    launch!(architecture, LinearWorkOrder, (npoints,), parameterization_tendencies_kernel!, diagn, progn, model)
    return nothing
end

@kernel function parameterization_tendencies_kernel!(diagn, progn, model)
    
    ij = @index(Global, Linear)     # every horizontal grid point ij

    perturb_inputs!(            ij, diagn, progn, model)

    albedo!(                    ij, diagn, progn, model)
    thermodynamics!(            ij, diagn, progn, model)
    temperature_relaxation!(    ij, diagn, progn, model)
    boundary_layer_drag!(       ij, diagn, progn, model)
    vertical_diffusion!(        ij, diagn, progn, model)
    convection!(                ij, diagn, progn, model)
    large_scale_condensation!(  ij, diagn, progn, model)
    optical_depth!(             ij, diagn, progn, model)
    shortwave_radiation!(       ij, diagn, progn, model)
    longwave_radiation!(        ij, diagn, progn, model)
    surface_fluxes!(            ij, diagn, progn, model)

    perturb_tendencies!(        ij, diagn, progn, model)
end

"""
$(TYPEDSIGNATURES)
Convert the fluxes on half levels to tendencies on full levels."""
function fluxes_to_tendencies!(
    column::ColumnVariables,
    geometry::Geometry,
    planet::AbstractPlanet,
    atmosphere::AbstractAtmosphere,
)
    
    (; u_tend, flux_u_upward, flux_u_downward) = column
    (; v_tend, flux_v_upward, flux_v_downward) = column
    (; humid_tend, flux_humid_upward, flux_humid_downward) = column
    (; temp_tend,  flux_temp_upward,  flux_temp_downward) = column

    Δσ = geometry.σ_levels_thick
    pₛ = column.pres[end]               # surface pressure
    (; radius) = planet                 # used for scaling

    # for g/Δp and g/(Δp*c_p), see Fortran SPEEDY documentation eq. (3, 5)
    g_pₛ = planet.gravity/pₛ
    cₚ = atmosphere.heat_capacity

    # fluxes are defined on half levels including top k=1/2 and surface k=nlayers+1/2
    @inbounds for k in eachindex(u_tend, v_tend, humid_tend, temp_tend)

        # Absorbed flux in a given layer, i.e. flux in minus flux out from above and below
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
        g_Δp_cₚ = g_Δp/cₚ
        u_tend[k] = radius*(u_tend[k] + g_Δp*ΔF_u)
        v_tend[k] = radius*(v_tend[k] + g_Δp*ΔF_v)
        humid_tend[k] = radius*(humid_tend[k] + g_Δp*ΔF_humid)
        temp_tend[k] = radius*(temp_tend[k] + g_Δp_cₚ*ΔF_temp)
    end

    return nothing
end