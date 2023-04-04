"""
    parameterization_tendencies!(   diagn::DiagnosticVariables,
                                    time::DateTime,
                                    M::PrimitiveEquation)

Compute tendencies for u,v,temp,humid from physical parametrizations.
Extract for each vertical atmospheric column the prognostic variables
(stored in `diagn` as they are grid-point transformed), loop over all
grid-points, compute all parametrizations on a single-column basis,
then write the tendencies back into a horizontal field of tendencies.
"""
function parameterization_tendencies!(  diagn::DiagnosticVariables,
                                        time::DateTime,
                                        model::PrimitiveEquation)

    G = model.geometry
    boundary_layer_scheme = model.parameters.boundary_layer
    temperature_relax_scheme = model.parameters.temperature_relaxation

    rings = eachring(G.Grid,G.nlat_half)

    @floop for ij in eachgridpoint(diagn)       # loop over all horizontal grid points

        thread_id = Threads.threadid()          # not two threads should use the same ColumnVariable
        column = diagn.columns[thread_id]
        jring = whichring(ij,rings)             # ring index gridpoint ij is on

        reset_column!(column)                   # set accumulators back to zero for next grid point
        get_column!(column,diagn,ij,jring,G)    # extract column for contiguous memory access

        # HELD-SUAREZ
        temperature_relaxation!(column,temperature_relax_scheme,model)
        boundary_layer!(column,boundary_layer_scheme,model)

        # Pre-compute thermodynamic quantities
        # get_thermodynamics!(column,model)

        # Calculate parametrizations (order of execution is important!)
        # convection!(column,model)
        # large_scale_condensation!(column,model)
        # clouds!(column, model)
        # shortwave_radiation!(column,model)
        # longwave_radiation!(column,model)
        # surface_fluxes!(column,model)
        # vertical_diffusion!(column,M)

        # write tendencies from parametrizations back into horizontal fields
        write_column_tendencies!(diagn,column,ij)
    end
end
