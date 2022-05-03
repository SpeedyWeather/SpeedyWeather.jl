# =========================================================================
# This file calculates the tendencies of the prognostic variables via the top-level function get_tendencies()
# The general structure is as follows:
#
#get_tendencies()      
#│
#└───get_grid_point_tendencies() ---> returns vor_tend, div_tend,temp_tend,pres_surf_tend,tr_tend + some extras
#│   │   
#│   │   
#│   │
#│   └───get_grid_point_fields() ---> Updates Diag.gridvars 
#|   └───parametrizationd_tendencies()
#|   └───dynamics_tendencies() ---> Returns Diag.Tendencies 
#│        │   
#│        │   
#│        │ 
#|        └───`src/dynamics_tendencies.jl`
#|        └───vor_div_tendency_and_corrections()
#│   
#└───get_spectral_tendencies() ---> returns pres_surf_tend,temp_tend,div_tend
#
#
# Note that get_spectral_tendencies() is quite badly named. It really modifies/updates the **already calculated** spectral tendencies
# =========================================================================


"""
Compute the grid point and spectral tendencies, including an implicit correction for the spectral tendencies. 
"""
function get_tendencies!(   diagn::DiagnosticVariables{NF}, # all diagnostic variables
                            progn::PrognosticVariables{NF}, # all prognostic variables
                            M::ModelSetup{NF},              # struct containing all constants
                            lf2::Int=2                      # leapfrog index 2 (time step used for tendencies)
                            ) where {NF<:AbstractFloat}

    # @unpack α = M.Parameters

    # =========================================================================
    # Computation of grid-point tendencies (converted to spectral at the end of
    # grtend) Diag.Tendencies.GridPoint
    # =========================================================================

    # println("Hello from get_tendencies")
    # get_grid_point_tendencies!(Prog,Diag,l2,M)

    # =========================================================================
    # Computation of spectral tendencies 
    # =========================================================================
    
    # @unpack u_grid, v_grid, vor_grid = diagn.grid_variables
    # @unpack vor_tend = diagn.tendencies

    divergence_uvω!(diagn,M.geospectral)


    # if  α == 0

    #     get_spectral_tendencies!(Prog,Diag,1,M)

    # else
    #     get_spectral_tendencies!(Prog,Diag,l2,M) #Note: l2 currently not used. Needs to be implemented 

    # end

    #if alpha < 0.5 #Coefficient for semi-implicit computations. Previously if alpha = 0?
       # get_spectral_tendencies!(Diag,C)
    #else
        
        #get_spectral_tendencies!(Diag,C)

        # Implicit correction
        #@unpack div_tend,t_tend,ps_tend = Diag.Tendencies #Get some of the tendencies calculated above
        #implicit_terms!(div_tend, t_tend, ps_tend) #TK: no edits have been made (yet!) to this implicit function
   # end


end



"""
Compute the grid point tendencies. These are composed of the dynamics and physical parameterisations. 
"""
function get_grid_point_tendencies!(Prog::PrognosticVariables{NF}, # Prognostic variables
                                    Diag::DiagnosticVariables{NF}, # Diagnostic variables
                                    l2::Int,                       # leapfrog index 2 (time step used for tendencies)
                                    M
                                    ) where {NF<:AbstractFloat}
    
    #1. Convert spectral prognostic variables to grid point space
    get_grid_point_fields!(Prog,Diag,l2,M) #Updates Diag.gridvars
 
    # 2. Parameterised physics tendencies. 
    #Needs to be defined in parameterisation_tendencies.jl. There is a lot of physics here, so would be best as separate, self-contained PR
    parametrization_tendencies!(Prog,Diag,M)


    #3. Dynamics tendencies
    dynamics_tendencies!(Prog,Diag,l2,M) #Takes Diag.gridvars and Diag.ParameterisedTendencies and calculates Diag.Tendencies 



end


"""
Compute grid point fields of the spectral prognostic tendencies
"""
function get_grid_point_fields!(Prog::PrognosticVariables{NF}, # Prognostic variables
                               Diag::DiagnosticVariables{NF}, # Diagnostic variables
                               l2::Int,                       # leapfrog index 2 (time step used for tendencies)
                               M
                               ) where {NF<:AbstractFloat}


    #Unpack spectral prognostic fields
    @unpack vor,div,temp,pres_surf,humid = Prog 

    #Unpack the gridded counterparts. This is what we will be calculating in this subroutine
    @unpack vor_grid,div_grid,temp_grid,pres_surf_grid,humid_grid, u_grid, v_grid = Diag.grid_variables

    #Unpack constants
    @unpack cp = M.Parameters 
    @unpack f_coriolis = M.GeoSpectral.geometry
    

    #Get dimensions of the grid. Can also read from M.Geometry
    nlon,nlat,nlev = size(vor_grid)

    #1. Compute grid-point fields
    #1.1 Update geopotential in spectral space
    #geopotential!(geopot,ϕ0spectral,temp,G)         # geopotential from surface geopotential. Need to get phi0_spectral and G from somewhere

   for k in 1:nlev
        vor_grid[:,:,k]  = gridded(vor[:,:,k])  # vorticity 
        div_grid[:,:,k]  = gridded(div[:,:,k])  # divergence
        temp_grid[:,:,k] = gridded(temp[:,:,k]) # temperature


        #Correct vorticity grid point field
        for j in 1:nlat 
            vor_grid[:,j,k] .+= f_coriolis[j]
        end



        #Calculate zonal velocity u and meridional velocity v in grid-point space,
        #from vorticity and divergence in spectral space.
        #This function is currently not well defined and returns all 0s
        uvspec!(vor[:,:,k], div[:,:,k], u_grid[:,:,k],v_grid[:,:,k])

        #Geopotential. Normalised by cp to avoid overflows. Feature from Paxton/Chantry
       # gridded(geopot[:,:,k]*(1/cp),geopot_grid[:,:,k]) 

   end



    #From Paxton/Chantry: "Truncate variables where the spectral transform is still done at double
    #precision". Or just conversion to per second?
    u_grid = u_grid / 3600.0
    v_grid = v_grid / 3600.0


    # =========================================================================
    # COMMENT: Need to verify l2 indexes here. In previous SPEEDY versions we have an index for
    # physical tendencies and a separate index for dynamical tendencies.
    # =========================================================================


    #Surface pressure
    pres_surf_grid = gridded(pres_surf)


    #Humidity
    #From Paxton/Chantry: "Don't transform the two stratospheric levels where humidity is set to zero
    # because it leads to overflows" 
    for k in 3:nlev
        humid_grid[:,:,k] = gridded(humid[:,:,k])
    end 


end

"""
Compute non-linear tendencies in grid-point space from dynamics and add to physics tendencies. Convert total
gridpoint tendencies to spectral tendencies.
"""
function dynamics_tendencies!(Prog::PrognosticVariables{NF}, # Prognostic variables
                             Diag::DiagnosticVariables{NF}, # Diagnostic variables
                             l2::Int,                       # leapfrog index 2 (time step used for tendencies)
                             M
                             ) where {NF<:AbstractFloat}
    

    # =========================================================================
    # All functions below are defined in `src/dynamics_tendencies.jl`
    # =========================================================================

    #1. Compute tendency of log(surface pressure)
    surface_pressure_tendency!(Prog,Diag,l2,M)               # Calculates pres_surf_tend and associated intermediate fields

    #2. Compute "vertical" velocity
    vertical_velocity_tendency!(Diag,M)                 # Calculates sigma_tend,sigma_m

    # 3. Subtract part of temperature field that is used as reference for implicit terms
    temperature_grid_anomaly!(Diag,M)                        # Calculates temp_grid_anomaly

    # 4. Zonal wind tendency
    zonal_wind_tendency!(Diag,M)                             # Calculates u_tend

    # 5. Meridional wind tendency
    meridional_wind_tendency!(Diag,M)                        # Calculates v_tend

    # 6. Temperature tendency
    temperature_tendency!(Diag,M)                            # Calculates temp_tend

    # 7. Humidity tendency
    humidity_tendency!(Diag,M)                               # Calculates humid_tend



    # =========================================================================
    # Calculate vor_tend,div_tend and then update temp_tend and humid_tend
    # =========================================================================
    vor_div_tendency_and_corrections!(Diag,M)


end

"""
Calculate vor_tend,div_tend and then update temp_tend and humid_tend
"""
function vor_div_tendency_and_corrections!( Diag::DiagnosticVariables{NF},
                                            M
                                            ) where {NF<:AbstractFloat}


    @unpack u_tend,v_tend,vor_tend,div_tend,temp_tend,humid_tend = Diag.tendencies
    @unpack u_grid,v_grid,temp_grid_anomaly,temp_grid,humid_grid = Diag.grid_variables
    _,_,nlev = size(u_grid)

    @unpack L2_velocity_complex = Diag.intermediate_variables



    for k in 1:nlev 

        #1. Calculate vor and div spectral tendencies from u and v grid tendencies 
        vdspec!(u_tend[:,:,k], v_tend[:,:,k], vor_tend[:,:,k], div_tend[:,:,k], true,M.GeoSpectral) 

        #2. Divergence tendency
        ## add -laplacian(0.5*(u**2+v**2)) to divergence tendency
        L2_velocity_complex = spectral(0.50*(u_grid[:,:,k].^2 + v_grid[:,:,k].^2))
        div_tend[:,:,k] = div_tend[:,:,k] - ∇²(L2_velocity_complex,M.GeoSpectral) 

        #3. Temperature tendency
        # add div(vT) to spectral t tendency
        vdspec!(-u_grid[:,:,k].*temp_grid_anomaly[:,:,k], 
                -v_grid[:,:,k].*temp_grid_anomaly[:,:,k], 
                L2_velocity_complex, temp_tend[:,:,k], 
                true,
                M.GeoSpectral)
    
        temp_tend[:,:,k] = temp_tend[:,:,k] + spectral(temp_grid[:,:,k])


        #4. Humidity (tracer) tendency
        vdspec!(-u_grid[:,:,k].*humid_grid[:,:,k], 
                    -v_grid[:,:,k].*humid_grid[:,:,k], 
                    L2_velocity_complex, humid_tend[:,:,k], 
                    true,
                    M.GeoSpectral)

        humid_tend[:,:,k] = humid_tend[:,:,k] + spectral(humid_grid[:,:,k])
         

    end 

 
                                
end 


"""
Compute the spectral tendencies of divergence, temperature
and log_surf.pressure 
"""
function get_spectral_tendencies!(Prog::PrognosticVariables{NF},
                                  Diag::DiagnosticVariables{NF},
                                  l2::Int,                       # leapfrog index 2 (time step used for tendencies)
                                  M
                                 ) where {NF<:AbstractFloat}





    @unpack σ_levels_thick = M.GeoSpectral.geometry # This is dhs in Fortran version
    @unpack div,temp,pres_surf= Prog 
    @unpack pres_surf_tend,temp_tend,div_tend = Diag.tendencies
    @unpack vertical_mean_divergence,sigdtc,dumk,spectral_geopotential = Diag.intermediate_variables
    @unpack tref,σ_levels_half⁻¹_2,tref2,tref3,rgas = M.GeoSpectral.geometry 
    @unpack geopot_surf = M.Boundaries 

    #Get number of levels
    _,_,nlev = size(div)


    #1. Pressure surface tendency 

    #Zero out pre-defined array 
    vertical_mean_divergence .= zero(NF)
    for k in 1:nlev

        vertical_mean_divergence = vertical_mean_divergence + div[:,:,k] * σ_levels_thick[k]

    end 

    
    pres_surf_tend = pres_surf_tend - vertical_mean_divergence
    pres_surf_tend[1,1] = zero(NF) 
    
    

    #2. Temperature tendency 
    sigdtc[:,:,1] .= zero(NF)
    sigdtc[:,:,nlev+1] .= zero(NF) #What exactly is this quantity?
    

    for k in 1:nlev-1
        sigdtc[:,:,k+1] = sigdtc[:,:,k] - σ_levels_thick[k] * (div[:,:,k] - vertical_mean_divergence)
    end 


    dumk[:,:,1] .= zero(NF)
    dumk[:,:,nlev+1] .= zero(NF) #What exactly is this quantity?

    for k in 2:nlev
        dumk[:,:,k] = sigdtc[:,:,k] * (tref[k] - tref[k-1])
    end

    for k in 1:nlev 
        temp_tend[:,:,k] = temp_tend[:,:,k] 
                          -(dumk[:,:,k+1] + dumk[:,:,k])*σ_levels_half⁻¹_2[k]
                          +(sigdtc[:,:,k+1]+sigdtc[:,:,k])*tref3[k]
                          -vertical_mean_divergence*tref2[k]
    end 



    #3. Divergence tendency 

    #First calculate the spectral_geopotential
    geopotential!(spectral_geopotential,geopot_surf,temp, M.GeoSpectral)

    #Then calculate the divergence tendency:
    for k in 1:nlev 
        div_tend[:,:,k] = div_tend[:,:,k] - ∇²(spectral_geopotential[:,:,k] + rgas*tref[k]*pres_surf,M.GeoSpectral) 
    end 


end

function add_tendencies!(   tend::AbstractMatrix{NF},   # tendency to accumulate into
                            term1::AbstractMatrix{NF},  # with term1
                            term2::AbstractMatrix{NF}   # and term2
                            ) where NF                  # number format real or complex

    # term1, term2 can have one more degree l which will be ignored in the loop though
    size_compat1 = size(tend) == size(term1) || (size(term1) .- (1,0)) == size(tend)
    size_compat2 = size(tend) == size(term2) || (size(term2) .- (1,0)) == size(tend)
    @boundscheck size_compat1 || throw(BoundsError)
    @boundscheck size_compat2 || throw(BoundsError)

    lmax,mmax = size(tend) .- 1

    @inbounds for m in 1:mmax+1
        for l in m:lmax+1
            tend[l,m] += (term1[l,m] + term2[l,m])
        end
    end
end

function add_tendencies!(   tend::AbstractMatrix{Complex{NF}},  # tendency to accumulate into
                            term::AbstractMatrix{Complex{NF}}   # with term
                            ) where NF                          # number format real or complex

    # term can have one more degree l which will be ignored in the loop though
    size_compat = size(tend) == size(term) || (size(term) .- (1,0)) == size(tend)
    @boundscheck size_compat || throw(BoundsError)
    lmax,mmax = size(tend) .- 1

    @inbounds for m in 1:mmax+1
        for l in m:lmax+1
            tend[l,m] += term[l,m]
        end
    end
end