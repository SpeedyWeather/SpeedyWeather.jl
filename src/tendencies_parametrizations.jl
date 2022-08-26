"""
    parametrization_tendencies!(diagn::DiagnosticVariables,
                                M::PrimitiveEquationModel)

Compute tendencies for u,v,temp,humid from physical parametrizations.
Extract for each vertical atmospheric column the prognostic variables
(stored in `diagn` as they are grid-point transformed), loop over all
grid-points, compute all parametrizations on a single-column basis,
then write the tendencies back into a horizontal field of tendencies.
"""
function parametrization_tendencies!(
    diagn::DiagnosticVariables{NF},
    M::PrimitiveEquationModel,
) where {NF}
    G = M.geometry
    column = ColumnVariables{NF}(nlev = diagn.nlev)

    for ij in eachgridpoint(diagn)      # loop over all horizontal grid points

        reset_column!(column)           # set accumulators back to zero for next grid point
        get_column!(column, diagn, ij, G)  # extract an atmospheric column for contiguous memory access

        # Pre-compute thermodynamic quantities
        get_thermodynamics!(column, M)

        # Calculate parametrizations (order of execution is important!)
        convection!(column, M)
        large_scale_condensation!(column, M)
        # clouds!(column, M)
        # shortwave_radiation!(column, M)
        # longwave_radiation!(column, M)
        # surface_fluxes!(column,M)
        # vertical_diffusion!(column,M)

        # write tendencies from parametrizations back into horizontal fields
        write_column_tendencies!(diagn, column, ij)
    end
end
