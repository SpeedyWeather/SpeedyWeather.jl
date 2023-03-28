struct ParameterizationConstants{NF<:AbstractFloat} <: AbstractParameterizationConstants{NF}
    # BOUNDARY LAYER
    drag_coefs::Vector{NF}          # (inverse) drag time scale per layer

    # TEMPERATURE RELAXATION
    temp_equilibrium::Matrix{NF}    # equilibrium temperature function of latitude and pressure

    # RADIATION
    # fband::Matrix{NF}
end

function Base.zeros(::Type{ParameterizationConstants},G::AbstractGeometry{NF}) where NF
    (;nlev,nlat) = G
    drag_coefs = zeros(NF,nlev)
    temp_equilibrium = zeros(NF,nlev,nlat)
    return ParameterizationConstants{NF}(   drag_coefs,
                                            temp_equilibrium,
                                        )
end

function ParameterizationConstants(P::Parameters,G::AbstractGeometry)
    K = zeros(ParameterizationConstants,P,G)
    initialise_boundary_layer!(K,P.boundary_layer,P,G)
    initialise_temperature_relaxation!(K,P.temperature_relaxation,P,G)
    # initialise_longwave_radiation!(K,P)
    return K
end
