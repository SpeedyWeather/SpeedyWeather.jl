"""
Compute physical parametrization tendencies.
"""
function parametrization_tendencies!(
    Prog::PrognosticVariables{NF},
    Diag::DiagnosticVariables{NF},
    M,
) where {NF<:AbstractFloat}
    get_large_scale_condensation_tendencies!(Diag, M)

    #Calculate
    #utend: u-wind tendency (gp)
    #vtend: v-wind tendency (gp)
    #ttend: temp. tendency (gp)
    #htend: spec. hum. tendency (gp)
end
