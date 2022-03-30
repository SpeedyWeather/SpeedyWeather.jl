struct parametrization_tendencies{T<:AbstractFloat}

#Just a placeholder

end

"""
Compute physical parametrization tendencies
"""
function parametrization_tendencies!(Prog::PrognosticVariables{NF}, # Prognostic variables
                                     Diag::DiagnosticVariables{NF}, # Diagnostic variables
                                     M
                                     ) where {NF<:AbstractFloat}


    #Calculate 
    #utend: u-wind tendency (gp)
    #vtend: v-wind tendency (gp)
    #ttend: temp. tendency (gp)
    #htend: spec. hum. tendency (gp)


end
