"""
Parameters for radiation parameterizations.
"""
struct ParameterizationConstants{NF<:AbstractFloat}
    fband::Matrix{NF}
end    

function ParameterizationConstants(P::Parameters)
    @unpack NF, nband = P
    fband = fill(NF(NaN), 400, nband)   # Energy fraction emitted in each longwave band = f(T)

    # conversion to number format NF happens here
    K = ParameterizationConstants{NF}(fband)
    initialise_longwave_radiation!(P, K)
    return K
end
