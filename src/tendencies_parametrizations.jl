"""
    parametrization_tendencies!(diagn::DiagnosticVariables,
                                M::PrimitiveEquationModel)

Compute tendencies for u,v,temp,humid from physical parametrizations.
Extract for each vertical atmospheric column the prognostic variables
(stored in `diagn` as they are grid-point transformed), loop over all
grid-points, compute all parametrizations on a single-column basis,
then write the tendencies back into a horizontal field of tendencies.
"""
function parametrization_tendencies!(diagn::DiagnosticVariables{NF},
                                     time::DateTime,
                                     model::PrimitiveEquationModel) where {NF}
    G = model.geometry
    column = ColumnVariables{NF}(nlev = diagn.nlev)

    for ij in eachgridpoint(diagn)      # loop over all horizontal grid points
        reset_column!(column)           # set accumulators back to zero for next grid point
        get_column!(column, diagn, ij, G)  # extract an atmospheric column for contiguous memory access

        # Pre-compute thermodynamic quantities
        get_thermodynamics!(column, model)

        # Calculate parametrizations (order of execution is important!)
        convection!(column, model)
        large_scale_condensation!(column, model)
        # clouds!(column, model)
        shortwave_radiation!(column, model)
        longwave_radiation!(column, model)
        # surface_fluxes!(column,model)
        # vertical_diffusion!(column,M)

        # write tendencies from parametrizations back into horizontal fields
        write_column_tendencies!(diagn, column, ij)
    end
end
