"""
    Δp_geopot_half, Δp_geopot_full = initialize_geopotential(   σ_levels_full::Vector,
                                                                σ_levels_half::Vector,
                                                                R_dry::Real)

Precomputes """
function initialize_geopotential(   σ_levels_full::Vector,
                                    σ_levels_half::Vector,
                                    R_dry::Real)

    nlev = length(σ_levels_full)    # number of vertical levels
    @assert nlev+1 == length(σ_levels_half) "σ half levels must have length nlev+1"

    Δp_geopot_half = zeros(nlev)    # allocate arrays
    Δp_geopot_full = zeros(nlev)

    # 1. integration onto half levels
    for k in 1:nlev-1               # k is full level index, 1=top, nlev=bottom
        # used for: Φ_{k+1/2} = Φ_{k+1} + R*T_{k+1}*(ln(p_{k+1}) - ln(p_{k+1/2}))
        Δp_geopot_half[k+1] = R_dry*log(σ_levels_full[k+1]/σ_levels_half[k+1])
    end

    # 2. integration onto full levels (same formula but k -> k-1/2)
    for k in 1:nlev
        # used for: Φ_k = Φ_{k+1/2} + R*T_k*(ln(p_{k+1/2}) - ln(p_k))
        Δp_geopot_full[k] = R_dry*log(σ_levels_half[k+1]/σ_levels_full[k])
    end

    return Δp_geopot_half, Δp_geopot_full
end

# function lapserate_correction(  σ_levels_full::Vector,
#                                 σ_levels_half::Vector,
#                                 Δp_geopot_full::Vector)

#     nlev = length(σ_levels_full)
#     @assert nlev+1 == length(σ_levels_half) "σ half levels must have length nlev+1"
#     @assert nlev == length(Δp_geopot_full) "σ half levels must have length nlev"

#     lapserate_corr = zeros(nlev)
#     for k in 2:nlev-1   # only in the free troposphere
#         # TODO reference
#         lapserate_corr[k] = 0.5*Δp_geopot_full[k]*
#                     log(σ_levels_half[k+1]/σ_levels_full[k]) / log(σ_levels_full[k+1]/σ_levels_full[k-1])
#     end

#     return lapserate_corr
# end

# """
#     geopotential!(diagn,B,G)

# Compute spectral geopotential `geopot` from spectral temperature `temp`
# and spectral surface geopotential `geopot_surf` (orography*gravity).
# """
# function geopotential!( diagn::DiagnosticVariables{NF},
#                         B::Boundaries{NF},      # contains surface geopotential
#                         G::Geometry{NF}         # contains precomputed layer-thickness arrays
#                         ) where NF              # number format NF

#     (;geopot_surf) = B.orography           # = orography*gravity
#     (;Δp_geopot_half, Δp_geopot_full) = G  # = R*Δlnp either on half or full levels
#     (;lapserate_corr) = G
#     (;nlev) = G                            # number of vertical levels

#     @boundscheck diagn.nlev == length(Δp_geopot_full) || throw(BoundsError)

#     # for PrimitiveDry virtual temperature = absolute temperature here
#     # note these are not anomalies here as they are only in grid-point fields
    
#     # BOTTOM FULL LAYER
#     temp = diagn.layers[end].dynamics_variables.temp_virt
#     geopot = diagn.layers[end].dynamics_variables.geopot
    
#     @inbounds for lm in eachharmonic(geopot,geopot_surf,temp)
#         geopot[lm] = geopot_surf[lm] + temp[lm]*Δp_geopot_full[end]
#     end

#     # OTHER FULL LAYERS, integrate two half-layers from bottom to top
#     @inbounds for k in nlev-1:-1:1
#         temp_k    = diagn.layers[k].dynamics_variables.temp_virt
#         temp_k1   = diagn.layers[k+1].dynamics_variables.temp_virt
#         geopot_k  = diagn.layers[k].dynamics_variables.geopot
#         geopot_k1 = diagn.layers[k+1].dynamics_variables.geopot

#         for lm in eachharmonic(temp_k,temp_k1,geopot_k,geopot_k1)
#             geopot_k½ = geopot_k1[lm] + temp_k1[lm]*Δp_geopot_half[k+1] # 1st half layer integration
#             geopot_k[lm] = geopot_k½  + temp_k[lm]*Δp_geopot_full[k]    # 2nd onto full layer
#         end      
#     end

#     # # LAPSERATE CORRECTION IN THE FREE TROPOSPHERE (>nlev)
#     # # TODO only for spectral coefficients 1,: ?
#     # lmax,mmax = size(geopot)
#     # for k in 2:nlev-1
#     #     temp_k_above = progn.layers[k-1].timesteps[lf].temp
#     #     temp_k_below = progn.layers[k+1].timesteps[lf].temp
#     #     geopot  = diagn.layers[k].dynamics_variables.geopot

#     #     for l in 1:lmax-1
#     #         geopot[l,l] += lapserate_corr[k]*(temp_k_below[l,l] - temp_k_above[l,l])
#     #     end
#     # end
# end

# function geopotential!( temp::Vector,
#                         G::Geometry)
    
#     (;Δp_geopot_half, Δp_geopot_full, nlev) = G  # = R*Δlnp either on half or full levels
#     geopot = zero(temp)

#     # bottom layer
#     geopot[nlev] = temp[nlev]*Δp_geopot_full[end]

#     # OTHER FULL LAYERS, integrate two half-layers from bottom to top
#     @inbounds for k in nlev-1:-1:1
#         geopot[k] = geopot[k+1] + temp[k+1]*Δp_geopot_half[k+1] + temp[k]*Δp_geopot_full[k]
#     end

#     return geopot
# end

# function geopotential!( diagn::DiagnosticVariablesLayer,
#                         pres::LowerTriangularMatrix,
#                         C::DynamicsConstants)

#     (;gravity) = C
#     (;geopot) = diagn.dynamics_variables

#     geopot .= pres*gravity
# end 

# """
#     virtual_temperature!(   diagn::DiagnosticVariablesLayer,
#                             temp::LowerTriangularMatrix,
#                             M::PrimitiveWet)
                
# Calculates the virtual temperature Tᵥ as

#     Tᵥ = T(1+μq)

# With absolute temperature T, specific humidity q and

#     μ = (1-ξ)/ξ, ξ = R_dry/R_vapour.
    
# In grid-point space and then transforms Tᵥ back into spectral space
# for the geopotential calculation."""
# function virtual_temperature!(  diagn::DiagnosticVariablesLayer,
#                                 ::LowerTriangularMatrix,    # only needed for dispatch compat with DryCore
#                                 model::PrimitiveWet)
    
#     (;temp_grid, humid_grid, temp_virt_grid) = diagn.grid_variables
#     (;temp_virt) = diagn.dynamics_variables
#     μ = model.constants.μ_virt_temp
#     S = model.spectral_transform

#     @inbounds for ij in eachgridpoint(temp_virt_grid, temp_grid, humid_grid)
#         temp_virt_grid[ij] = temp_grid[ij]*(1 + μ*humid_grid[ij])
#     end
#     # spectral!(temp_virt,temp_virt_grid,S)
# end

# function linear_virtual_temperature!(   diagn::DiagnosticVariablesLayer,
#                                         progn::PrognosticLayerTimesteps,
#                                         model::PrimitiveWet,
#                                         lf::Int)
    
#     (;temp_virt) = diagn.dynamics_variables
#     μ = model.constants.μ_virt_temp
#     Tₖ = model.geometry.temp_ref_profile[diagn.k]   
#     (;temp,humid) = progn.timesteps[lf]

#     @. temp_virt = temp + Tₖ*μ*humid
# end


# """
# For the PrimitiveDry temperautre and virtual temperature are the same (humidity=0).
# Just copy over the arrays."""
# function virtual_temperature!(  diagn::DiagnosticVariablesLayer,
#                                 temp::LowerTriangularMatrix,
#                                 ::PrimitiveDry)
    
#     (;temp_grid, temp_virt_grid) = diagn.grid_variables
#     (;temp_virt) = diagn.dynamics_variables

#     copyto!(temp_virt_grid,temp_grid)
#     # copyto!(temp_virt,temp)
# end

# function linear_virtual_temperature!(   diagn::DiagnosticVariablesLayer,
#                                         progn::PrognosticLayerTimesteps,
#                                         ::PrimitiveDry,
#                                         lf::Int)
    
#     (;temp_virt) = diagn.dynamics_variables
#     (;temp) = progn.timesteps[lf]

#     copyto!(temp_virt,temp)
# end
