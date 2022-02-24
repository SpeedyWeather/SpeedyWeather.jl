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
# We are also currently dealing with tracers rather than e.g. humidity explicitly. For now, lets scrap tracers and just deal with humidity directly?
# Constants struct C will probably really be a model setup var M. See L22 of run_speedy.jl
# =========================================================================


"""
Compute the grid point and spectral tendencies, including an implicit correction for the spectral tendencies. 
"""
function get_tendencies!(Prog::PrognosticVariables{NF}, # Prognostic variables
                         Diag::DiagnosticVariables{NF}, # Diagnostic variables
                         l2::Int,                       # leapfrog index 2 (time step used for tendencies)
                         M,              # struct containing constants
                        ) where {NF<:AbstractFloat}

   # @unpack alpha = C

    # =========================================================================
    # Computation of grid-point tendencies (converted to spectral at the end of
    # grtend) Diag.Tendencies.GridPoint
    # =========================================================================

    println("Hello from get_tendencies")
    get_grid_point_tendencies!(Prog,Diag,l2,M)

    # =========================================================================
    # Computation of spectral tendencies 
    # =========================================================================
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
   # parametrization_tendencies!(Prog,Diag,M)


    #3. Dynamics tendencies
   # dynamics_tendencies!(Prog,Diag,M) #Takes Diag.gridvars and Diag.ParameterisedTendencies and calculates Diag.Tendencies 



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
    @unpack cp,coriol = M.P #model.parameters

    nlat,nlon,nlev = size(vor_grid)


    #1. Compute grid-point fields
    #1.1 Update geopotential in spectral space
    #geopotential!(geopot,ϕ0spectral,temp,G)         # geopotential from surface geopotential. Need to get phi0_spectral and G from somewhere

   for k in 1:nlev
        gridded(vor[:,:,k,l2], vor_grid[:,:,k])  # vorticity 
        gridded(div[:,:,k,l2], div_grid[:,:,k])  # divergence
        gridded(temp[:,:,k,l2],temp_grid[:,:,k]) # temperature


        #Correct vorticity grid point field
        for j in 1:nlat
            for i in 1:nlon
                vor_grid[i,j,k] += coriol[j] 
            end
        end

        #Calculate zonal velocity u and meridional velocity v in grid-point space,
        #from vorticity and divergence in spectral space
        uvspec!(vor[:,:,k,l2], div[:,:,k,l2], u_grid[:,:,k],v_grid[:,:,k])

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
    gridded(pres_surf[:,:,l2], pres_surf_grid)


    #Humidity
    #From Paxton/Chantry: "Don't transform the two stratospheric levels where humidity is set to zero
    # because it leads to overflows" 
    for k in 3:nlev
        gridded(humid[:,:,k,l2,1], humid_grid[:,:,k])
    end 


end

"""
Compute non-linear tendencies in grid-point space from dynamics and add to physics tendencies. Convert total
gridpoint tendencies to spectral tendencies.
"""
function dynamics_tendencies(Prog::PrognosticVariables{NF}, # Prognostic variables
                             Diag::DiagnosticVariables{NF}, # Diagnostic variables
                             l2::Int,                       # leapfrog index 2 (time step used for tendencies)
                             C::Constants{NF}
                             ) where {NF<:AbstractFloat}
    

    # =========================================================================
    # All functions below are defined in `src/dynamics_tendencies.jl`
    # =========================================================================

    #1. Compute tendency of log(surface pressure)
    surface_pressure_tendency!(Prog,Diag,l2,C)               # Calculates pres_surf_tend

    #2. Compute "vertical" velocity
    vertical_velocity_tendency!(Prog,Diag,C)                 # Calculates sigma_tend,sigma_m

    # 3. Subtract part of temperature field that is used as reference for implicit terms
    temperature_grid_anomaly!(Diag,C)                        # Calculates temp_grid_anomaly

    # 4. Zonal wind tendency
    zonal_wind_tendency!(Diag,C)                             # Calculates u_tend

    # 5. Meridional wind tendency
    meridional_wind_tendency!(Diag,C)                        # Calculates v_tend

    # 6. Temperature tendency
    temperature_tendency!(Diag,C)                            # Calculates temp_tend

    # 7. Humidity tendency
    humidity_tendency!(Diag,C)                               # Calculates humid_tend



    # =========================================================================
    # Calculate vor_tend,div_tend and then update temp_tend and tr_tend
    # =========================================================================
    vor_div_tendency_and_corrections!(Diag,C)


end

"""
Convert a set of tendencies in grid point space to spectral space and calculate ...
"""
function vor_div_tendency_and_corrections!( Diag::DiagnosticVariables{NF},
                                            C::Constants{NF}
                                            ) where {NF<:AbstractFloat}

    @unpack u_tend, v_tend, vor_tend,div_tend = Diag.tendencies
    @unpack u_grid,v_grid,temp_grid_anomaly,tr_grid = Diag.gridvars
    @unpack dumc = Diag.miscvars
    
    _,_,nlev = size(u_tend)

    for k in 1:nlev
        #  1. Calculate vor and div spectral tendencies from u and v tendencies
        vdspec!(u_tend[:,:,k], v_tend[:,:,k], vor_tend[:,:,k], div_tend[:,:,k], true) 

        #2. Divergence tendency
        #add -laplacian(0.5*(u**2+v**2)) to divergence tendency
        dumc[:,:,1] = spectral(0.50*(u_grid[:,:,k]^2 + v_grid[:,:,k]^2))
        div_tend[:,:,k] = div_tend[:,:,k] - laplacian(dumc[:,:,1]) 

        
        #3. Temperature tendency
        # add div(vT) to spectral t tendency
        vdspec!(-u_grid[:,:,k]*temp_grid_anomaly[:,:,k], 
               -v_grid[:,:,k]*temp_grid_anomaly[:,:,k], 
               dumc[:,:,1], temp_tend[:,:,k], 
               true)
    
        temp_tend[:,:,k] = temp_tend[:,:,k] + spectral(t_grid[:,:,k])

        #4. Tracer tendency
        for itr in 1:n_trace
            vdspec!(-u_grid[:,:,k]*tr_grid[:,:,k,itr], 
                    -v_grid[:,:,k]*tr_grid[:,:,k,itr], 
                    dumc[:,:,1], tr_tend[:,:,k,itr], 
                    true)
        tr_tend[:,:,k,itr] = tr_tend[:,:,k,itr] + spectral(tr_grid[:,:,k,itr])
        end 

    end



end 


"""
Compute the spectral tendencies of divergence, temperature
and log_surf.pressure 
"""
function get_spectral_tendencies!(Prog::PrognosticVariables{NF},
                                  Diag::DiagnosticVariables{NF},
                                  l2::Int,                       # leapfrog index 2 (time step used for tendencies)
                                  C::Constants{NF}
                                 ) where {NF<:AbstractFloat}


    @unpack div,temp,pres_surf = Prog
    @unpack pres_surf_tend,temp_tend,div_tend = Diag.tendencies
    @unpack d_meanc, sigma_tend_c,dumk = Diag.miscvars
    @unpack dhs,dhsr,temp_ref,tref2,tref3 = C

    _,_,nlev = size(div)

    # 1. Vertical mean divergence and pressure tendency 
    d_meanc[:,:] = 0.0
    for k in 1:nlev
        d_meanc = d_meanc + div[:,:,k,l2]*dhs[k]
    end

    pres_surf_tend = pres_surf_tend - d_meanc
    pres_surf_tend[1,1] = Complex{RealType}(0.0)

    # 2. Sigma-dot "velocity" and temperature tendency
    for k in 1:nlev - 1
        sigma_tend_c[:,:,k+1] = sigma_tend_c[:,:,k] - dhs[k]*(div[:,:,k,j2] - d_meanc)
    end

    for k in 2:nlev
        dumk[:,:,k] = sigma_tend_c[:,:,k]*(temp_ref[k] - temp_ref[k-1])
    end

    for k in 1:nlev
        temp_tend[:,:,k] -= (dumk[:,:,k+1] + dumk[:,:,k])*dhsr[k]
            + tref3[k]*(sigma_tend_c[:,:,k+1] + sigma_tend_c[:,:,k]) - tref2[k]*d_meanc
    end

    # 3. Geopotential and divergence tendency
    geopotential!(geopot,ϕ0spectral,temp,G) # Paxton/Chantry have a geopotential call here. Do we actually need this?       
    for k in 1:nlev
        div_tend[:,:,k] -= ∇²(geopot[:,:,k] + rgas*temp_ref[k]*pres_surf[:,:,l2])
    end
end