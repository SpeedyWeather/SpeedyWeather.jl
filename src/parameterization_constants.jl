"""
Parameters for radiation parameterizations.
"""
struct ParameterizationConstants{NF <: AbstractFloat}
    fband::Matrix{NF}
end

function ParameterizationConstants(P::Parameters)
    @unpack NF, nband = P
    fband = nans(NF, 400, nband)   # Energy fraction emitted in each longwave band = f(T)

    K = ParameterizationConstants(fband)
    initialise_longwave_radiation!(K, P)
    return K
end
