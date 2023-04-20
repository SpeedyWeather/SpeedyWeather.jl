struct ParameterizationConstants{NF<:AbstractFloat} <: AbstractParameterizationConstants{NF}
    # BOUNDARY LAYER
    drag_coefs::Vector{NF}          # (inverse) drag time scale per layer

    # TEMPERATURE RELAXATION (Held and Suarez, 1996)
    temp_relax_freq::Matrix{NF}     # (inverse) relaxation time scale per layer and latitude
    temp_equil_a::Vector{NF}        # two terms to calculate equilibrium temperature as function
    temp_equil_b::Vector{NF}        # of latitude and pressure

    # TEMPERATURE RELAXATION (Jablonowski and Williamson, 2006)
    temp_equil::Matrix{NF}

    # VERTICAL DIFFUSION
    vert_diff_∇²_above::Vector{NF}
    vert_diff_∇²_below::Vector{NF}  # 
    vert_diff_ν::Vector{NF}         # diffusion coefficient to be updated for every column
    vert_diff_Δσ::Vector{NF}

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

    temp_equil = zeros(NF,nlev,nlat)

    vert_diff_∇²_above = zeros(NF,nlev-1)   # defined on half levels, top, and bottom are =0 though
    vert_diff_∇²_below = zeros(NF,nlev-1)   # defined on half levels, top, and bottom are =0 though
    vert_diff_ν = zeros(NF,nlev-1)          # diffusion coefficient, to be reused
    vert_diff_Δσ = zeros(NF,nlev-1)         # vertical gradient operator wrt σ coordinates

    fband = zeros(NF,400,nband)

    return ParameterizationConstants{NF}(   drag_coefs,
                                            temp_relax_freq,
                                            temp_equil_a,
                                            temp_equil_b,
                                            temp_equil,
                                            vert_diff_∇²_above,
                                            vert_diff_∇²_below,
                                            vert_diff_ν,
                                            vert_diff_Δσ,
                                            fband,
                                        )
end

function ParameterizationConstants(P::Parameters,G::AbstractGeometry)
    K = zeros(ParameterizationConstants,P,G)
    initialize_boundary_layer!(K,P.boundary_layer,P,G)
    initialize_temperature_relaxation!(K,P.temperature_relaxation,P,G)
    initialize_vertical_diffusion!(K,P.vertical_diffusion,P,G)
    initialize_longwave_radiation!(K,P)
    return K
end
