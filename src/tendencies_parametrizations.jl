"""
    parametrization_tendencies!(diagn::DiagnosticVariables,
                                M::PrimitiveEquationModel)

Compute tendencies for u,v,temp,humid from physical parametrizations.
Extract for each vertical atmospheric column the prognostic variables
(stored in `diagn` as they are grid-point transformed), loop over all
grid-points, compute all parametrizations on a single-column basis,
then write the tendencies back into a horizontal field of tendencies.
"""
function parametrization_tendencies!(   diagn::DiagnosticVariables{NF},
                                        M::PrimitiveEquationModel,
                                        ) where NF

    G = M.geometry
    column = ColumnVariables{NF}(nlev=diagn.nlev)

    for ij in eachgridpoint(diagn)

        # extract a single atmospheric column for contiguous memory access
        get_column!(column,diagn,ij,G)

        fill!(column.u_tend,0)      # set tendencies to 0 to accumulate +=
        fill!(column.v_tend,0)      # parametrizations into
        fill!(column.temp_tend,0)
        fill!(column.humid_tend,0)

        # calculate parametrizations
        # convection!(column, M)
        large_scale_condensation!(column, M)
        # clouds!(column, M)
        # shortwave_radiation!(column, M)
        # longwave_radiation!(column, M)
        # surface_fluxes!(column,M)
        # vertical_diffusion!(column,M)

        # write tendencies from parametrizations back into horizontal fields
        write_column_tendencies!(diagn,column,ij)
    end
end
