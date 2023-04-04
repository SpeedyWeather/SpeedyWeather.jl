struct ParameterizationConstants{NF<:AbstractFloat} <: AbstractParameterizationConstants{NF}
    # BOUNDARY LAYER
    drag_coefs::Vector{NF}          # (inverse) drag time scale per layer

    # TEMPERATURE RELAXATION (Held and Suarez, 1996)
    temp_relax_freq::Matrix{NF}     # (inverse) relaxation time scale per layer and latitude
    temp_equil_a::Vector{NF}        # two terms to calculate equilibrium temperature as function
    temp_equil_b::Vector{NF}        # of latitude and pressure

    # RADIATION
    fband::Matrix{NF}
end

function Base.zeros(::Type{ParameterizationConstants},
                    P::Parameters,
                    G::AbstractGeometry{NF}) where NF
    (;nlev,nlat) = G
    (;nband) = P

    drag_coefs = zeros(NF,nlev)     # coefficient for boundary layer drag
    
    temp_relax_freq = zeros(NF,nlev,nlat)   # (inverse) relaxation time scale per layer and latitude
    temp_equil_a = zeros(NF,nlat)   # terms to calculate equilibrium temperature as function
    temp_equil_b = zeros(NF,nlat)   # of latitude and pressure

    fband = zeros(NF,400,nband)

    return ParameterizationConstants{NF}(   drag_coefs,
                                            temp_relax_freq,
                                            temp_equil_a,
                                            temp_equil_b,
                                            fband,
                                        )
end

function ParameterizationConstants(P::Parameters,G::AbstractGeometry)
    K = zeros(ParameterizationConstants,P,G)
    initialize_boundary_layer!(K,P.boundary_layer,P,G)
    initialize_temperature_relaxation!(K,P.temperature_relaxation,P,G)
    initialize_longwave_radiation!(K,P)
    return K
end
