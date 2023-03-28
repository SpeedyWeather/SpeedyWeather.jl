struct ParameterizationConstants{NF<:AbstractFloat} <: AbstractParameterizationConstants{NF}
    # BOUNDARY LAYER
    drag_coefs::Vector{NF}     # (inverse) drag time scale per layer

    # RADIATION
    # fband::Matrix{NF}
end

function Base.zeros(::Type{ParameterizationConstants},P::Parameters)
    (;NF,nlev) = P
    drag_coefs = zeros(NF,nlev)
    return ParameterizationConstants{NF}(drag_coefs)
end

function ParameterizationConstants(P::Parameters,G::AbstractGeometry)
    K = zeros(ParameterizationConstants,P)
    initialise_boundary_layer!(K,P.boundary_layer,P,G)
    # initialise_longwave_radiation!(K,P)
    return K
end
