"""
Parameters for radiation parameterizations.
"""
@with_kw struct ParameterizationConstants{NF<:AbstractFloat}
    fband::Matrix{Real}
end    

function ParameterizationConstants(P::Parameters)
    @unpack NF, nband = P
    fband = fill(NF(NaN), 400, nband)   # Energy fraction emitted in each longwave band = f(T)

    # conversion to number format NF happens here
    K = ParameterizationConstants{NF}(fband)
    initialise_longwave_radiation!(P, K)
    return K
end
