"""
    geopotential!(diagn,progn,B,G)

Compute spectral geopotential `geopot` from spectral temperature `temp`
and spectral surface geopotential `geopot_surf` (orography*gravity).
"""
function geopotential!( diagn::DiagnosticVariables{NF},
                        progn::PrognosticVariables{NF},
                        lf::Int,            # leapfrog step
                        B::Boundaries{NF},  # contains surface geopotential
                        G::Geometry{NF}     # contains precomputed layer-thickness arrays
                        ) where NF          # number format NF

    @unpack geopot_surf = B
    @unpack Δp_geopot_half, Δp_geopot_full, lapserate_correction = G
    @unpack nlev = G

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
            geopot_k[lm] = geopot[i,j,k+1] + xgeop2[k+1]*temp[i,j,k+1] + xgeop1[k]*temp[i,j,k]
            # geopot[i,j,k] = geopot[i,j,k+1] + xgeop2[k+1]*temp[i,j,k+1] + xgeop1[k]*temp[i,j,k]
        end      
    end

    # LAPSERATE CORRECTION IN THE FREE TROPOSPHERE (>nlev)
    # TODO only for spectral coefficients 1,: ?
    for k in 2:nlev-1
        for j in 1:nx
            geopot[1,j,k] = geopot[1,j,k] +
                    lapserate_correction[k-1]*(temp[1,j,k+1] - temp[1,j,k-1])
        end
    end
end
