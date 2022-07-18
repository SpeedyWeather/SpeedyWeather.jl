"""
Compute physical parametrization tendencies.
"""
function parametrization_tendencies!(
    Prog::PrognosticVariables{NF},
    Diag::DiagnosticVariables{NF},
    M,
) where {NF<:AbstractFloat}
    @unpack pres_grid = Diag.grid_variables

    # The prognostic variable pres has units of log(hPa), whereas the grid-point pressure
    # field which we require for the parametrizations has units of hPa, so we take the
    # exponential here.
    @. pres_grid = exp(pres_grid)

    get_large_scale_condensation_tendencies!(Diag, M)

    #Calculate
    #utend: u-wind tendency (gp)
    #vtend: v-wind tendency (gp)
    #ttend: temp. tendency (gp)
    #htend: spec. hum. tendency (gp)
end
