"""
    Δp_geopot_half, Δp_geopot_full = initialise_geopotential(   σ_levels_full::Vector,
                                                                σ_levels_half::Vector,
                                                                R_gas::Real)

Precomputes """
function initialise_geopotential(   σ_levels_full::Vector,
                                    σ_levels_half::Vector,
                                    R_gas::Real)

    nlev = length(σ_levels_full)    # number of vertical levels
    @assert nlev+1 == length(σ_levels_half) "σ half levels must have length nlev+1"

    Δp_geopot_half = zeros(nlev)    # allocate arrays
    Δp_geopot_full = zeros(nlev)

    # 1. integration onto half levels
    for k in 1:nlev-1               # k is full level index, 1=top, nlev=bottom
        # used for: Φ_{k+1/2} = Φ_{k+1} + R*T_{k+1}*(ln(p_{k+1}) - ln(p_{k+1/2}))
        Δp_geopot_half[k+1] = R_gas*log(σ_levels_full[k+1]/σ_levels_half[k+1])
    end

    for k in 1:nlev
        # 2. integration onto full levels (same formula but k -> k-1/2)
        # used for: Φ_k = Φ_{k+1/2} + R*T_k*(ln(p_{k+1/2}) - ln(p_k))
        Δp_geopot_full[k] = R_gas*log(σ_levels_half[k+1]/σ_levels_full[k])
    end

    return Δp_geopot_half, Δp_geopot_full
end

"""
    geopotential!(diagn,progn,B,G)

Compute spectral geopotential `geopot` from spectral temperature `temp`
and spectral surface geopotential `geopot_surf` (orography*gravity).
"""
function geopotential!( diagn::DiagnosticVariables{NF},
                        progn::PrognosticVariables{NF},
                        lf::Int,                # leapfrog step
                        B::Boundaries{NF},      # contains surface geopotential
                        G::Geometry{NF}         # contains precomputed layer-thickness arrays
                        ) where NF              # number format NF

    @unpack geopot_surf = B                     # = orography*gravity
    @unpack Δp_geopot_half, Δp_geopot_full = G  # = R*Δlnp either on half or full levels
    @unpack lapserate_correction = G
    @unpack nlev = G                            # number of vertical levels

    # BOTTOM FULL LAYER
    temp = progn.layers[end].leapfrog[lf].temp
    geopot = diagn.layers[end].dynamics_variables.geopot
    
    for lm in eachharmonic(geopot,geopot_surf,temp)
        geopot[lm] = geopot_surf[lm] + Δp_geopot_full[end]*temp[lm]
    end

    # OTHER FULL LAYERS, integrate two half-layers from bottom to top
    for k in nlev-1:-1:1
        temp_k    = progn.layers[k].leapfrog[lf].temp
        temp_k1   = progn.layers[k+1].leapfrog[lf].temp
        geopot_k  = diagn.layers[k].dynamics_variables.geopot
        geopot_k1 = diagn.layers[k+1].dynamics_variables.geopot

        for lm in eachharmonic(temp_k,temp_k1,geopot_k,geopot_k1)
            geopot_k½ = geopot_k1[lm] + temp_k1[lm]*Δp_geopot_half[k+1] # 1st half layer integration
            geopot_k[lm] = geopot_k½  + temp_k[lm]*Δp_geopot_full[k]    # 2nd onto full layer
        end      
    end

    # # LAPSERATE CORRECTION IN THE FREE TROPOSPHERE (>nlev)
    # # TODO only for spectral coefficients 1,: ?
    # for k in 2:nlev-1
    #     for j in 1:nx
    #         geopot[1,j,k] = geopot[1,j,k] +
    #                 lapserate_correction[k-1]*(temp[1,j,k+1] - temp[1,j,k-1])
    #     end
    # end
end
