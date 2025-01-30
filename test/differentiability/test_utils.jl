# additional utitility functions 

function Base.one(diag::DiagnosticVariables{NF}) where NF
    vec, re = to_vec(diag)
    vec .= NF(1)
    return re(vec)
end 