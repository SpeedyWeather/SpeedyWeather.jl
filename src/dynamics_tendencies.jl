"""
Compute the spectral tendency of the surface pressure logarithm
"""
function surface_pressure_tendency!(Prog::PrognosticVariables{NF}, # Prognostic variables
                                    Diag::DiagnosticVariables{NF}, # Diagnostic variables
                                    l2::Int,                       # leapfrog index 2 (time step used for tendencies)
                                    M
                                    ) where {NF<:AbstractFloat}


    @unpack pres_surf                            = Prog 
    @unpack pres_surf_tend                       = Diag.tendencies
    @unpack u_grid,v_grid,div_grid               =   Diag.grid_variables
    @unpack u_mean,v_mean,div_mean,
            pres_surf_gradient_spectral_x,
            pres_surf_gradient_spectral_y,
            pres_surf_gradient_grid_x,
            pres_surf_gradient_grid_y            = Diag.intermediate_variables
    @unpack σ_levels_thick                       = M.GeoSpectral.geometry #I think this is dhs

    _,_,nlev = size(u_grid)


    
  
    #Calculate mean fields
    for k in 1:nlev
        u_mean += u_grid[:,:,k]  *σ_levels_thick[k] 
        v_mean += v_grid[:,:,k]  *σ_levels_thick[k]
        div_mean += div_grid[:,:,k]*σ_levels_thick[k]
    end

    #Now use the mean fields
    grad!(pres_surf, pres_surf_gradient_spectral_x, pres_surf_gradient_spectral_y, M.GeoSpectral)
    pres_surf_gradient_grid_x = gridded(pres_surf_gradient_spectral_x*3600, M.GeoSpectral)
    pres_surf_gradient_grid_y = gridded(pres_surf_gradient_spectral_x*3600, M.GeoSpectral) #3600 factor from Paxton/Chantry. I think this is to correct for the underflow rescaling earlier

    pres_surf_tend = spectral(-u_mean.*pres_surf_gradient_grid_x - v_mean.*pres_surf_gradient_grid_y,M.GeoSpectral)
    pres_surf_tend[1,1] = pres_surf_tend[1,1]*0.0 

end

"""
Compute the spectral tendency of the "vertical" velocity
"""
function vertical_velocity_tendency!(Diag::DiagnosticVariables{NF}, # Diagnostic variables
                                     M
                                     ) where {NF<:AbstractFloat}

    @unpack u_grid,v_grid,div_grid = Diag.grid_variables
    @unpack u_mean,v_mean,div_mean,pres_surf_gradient_grid_x,pres_surf_gradient_grid_y,sigma_tend,sigma_m, puv = Diag.intermediate_variables
    @unpack σ_levels_thick = M.GeoSpectral.geometry
    _,_,nlev = size(u_grid)


    for k in 1:nlev
        puv[:,:,k] = (u_grid[:,:,k] - u_mean) .* pres_surf_gradient_grid_x + (v_grid[:,:,k] - v_mean) .* pres_surf_gradient_grid_y
    end

    for k in 1:nlev
        sigma_tend[:,:,k+1] = sigma_tend[:,:,k] - σ_levels_thick[k]*(puv[:,:,k] + div_grid[:,:,k] - div_mean)
        sigma_m[:,:,k+1]    = sigma_m[:,:,k]    - σ_levels_thick[k]*puv[:,:,k]
    end

end



"""
Compute the temperature anomaly in grid point space
"""
function temperature_grid_anomaly!(Diag::DiagnosticVariables{NF}, # Diagnostic variables
                                   M
                                   ) where {NF<:AbstractFloat}

    @unpack temp_grid,temp_grid_anomaly = Diag.grid_variables
    @unpack tref = M.GeoSpectral.geometry #Note that tref is currenrtly not defined correctly 

    _,_,nlev = size(temp_grid)


    for k in 1:nlev
        temp_grid_anomaly[:,:,k] = temp_grid[:,:,k] .- tref[k] 
    end

end





"""
Compute the spectral tendency of the zonal wind
"""
function zonal_wind_tendency!(Diag::DiagnosticVariables{NF}, # Diagnostic variables
                              M
                              )where {NF<:AbstractFloat}
    
    @unpack u_tend = Diag.tendencies
    @unpack u_grid,v_grid,vor_grid,temp_grid_anomaly= Diag.grid_variables
    @unpack sigma_tend,pres_surf_gradient_grid_x,pres_surf_gradient_grid_y,sigma_u = Diag.intermediate_variables
    
    
    @unpack rgas,σ_levels_half⁻¹_2 = M.GeoSpectral.geometry #I think this is dhsr 

    _,_,nlev = size(u_grid)


    #Update px,py

    pres_surf_gradient_grid_x = rgas*pres_surf_gradient_grid_x
    pres_surf_gradient_grid_y = rgas*pres_surf_gradient_grid_y

   

    for k in 2:nlev
        sigma_u[:,:,k] = sigma_tend[:,:,k].*(u_grid[:,:,k] - u_grid[:,:,k-1])
    end


    for k in 1:nlev
        u_tend[:,:,k] = u_tend[:,:,k] + v_grid[:,:,k].*vor_grid[:,:,k] 
                        - temp_grid_anomaly[:,:,k].*pres_surf_gradient_grid_x
                        - (sigma_u[:,:,k+1] + sigma_u[:,:,k])*σ_levels_half⁻¹_2[k]
    end


  

end



"""
Compute the spectral tendency of the meridional wind 
"""
function meridional_wind_tendency!(Diag::DiagnosticVariables{NF}, # Diagnostic variables
                                   M
                                  )where {NF<:AbstractFloat}

    @unpack v_tend = Diag.tendencies
    @unpack vor_grid,u_grid,v_grid,temp_grid_anomaly =Diag.grid_variables
    @unpack sigma_tend,sigma_u,pres_surf_gradient_grid_x,pres_surf_gradient_grid_y = Diag.intermediate_variables
    
    
    @unpack rgas,σ_levels_half⁻¹_2 = M.GeoSpectral.geometry #I think this is dhsr 


    _,_,nlev = size(u_grid)


    for k in 2:nlev
        sigma_u[:,:,k] = sigma_tend[:,:,k].*(v_grid[:,:,k] - v_grid[:,:,k-1])
    end
          
 

    for k in 1:nlev
        v_tend[:,:,k] = v_tend[:,:,k] + u_grid[:,:,k].*vor_grid[:,:,k] 
                        - temp_grid_anomaly[:,:,k].*pres_surf_gradient_grid_y
                       - (sigma_u[:,:,k+1] + sigma_u[:,:,k])*σ_levels_half⁻¹_2[k]
    end




end





"""
Compute the spectral temperature tendency
"""
function temperature_tendency!(Diag::PrognosticVariables{NF}, # Diagnostic variables
                               C::Constants{NF}
                               )where {NF<:AbstractFloat}
    
    @unpack temp_tend = Diag.Tendencies
    @unpack div_grid,temp_grid_anomaly =Diag.gridvars
    @unpack d_mean,arbitrary_array,sigma_tend,sigma_m,puv= Diag.miscvars
    @unpack temp_ref,tref3,dhsr,fsgr,akap = C

    _,_,nlev = size(div_grid)



    for k in 2:nlev
        arbitrary_array[:,:,k] = sigma_tend[:,:,k].*(temp_grid_anomaly[:,:,k] - temp_grid_anomaly[:,:,k-1])
                    + sigma_m[:,:,k].*(temp_ref[k] - temp_ref[k-1])
    end

    for k in 1:nlev
        temp_tend[:,:,k] = temp_tend[:,:,k]
                        + temp_grid_anomaly[:,:,k].*div_grid[:,:,k]
                        - (arbitrary_array[:,:,k+1] + arbitrary_array[:,:,k])*dhsr[k]
                        + fsgr[k]*temp_grid_anomaly[:,:,k].*(sigma_tend[:,:,k+1] + sigma_tend[:,:,k])
                        + tref3[k]*(sigma_m[:,:,k+1] + sigma_m[:,:,k])
                        + akap*(t_grid[:,:,k].*puv[:,:,k] - temp_grid_anomaly[:,:,k].*d_mean)
    end 


end




"""
Compute the humidity tendency
"""
function humidity_tendency!(Diag::PrognosticVariables{NF}, # Diagnostic variables
                          C::Constants{NF}
                          )where {NF<:AbstractFloat}
    

    @unpack div_grid,humid_grid = Diag.gridvars
    @unpack arbitrary_array,sigma_tend= Diag.miscvars
    @unpack dhsr = C


    _,_,nlev = size(div_grid)


    for k in 2:nlev
        arbitrary_array[:,:,k] = sigma_tend[:,:,k].*(humid_grid[:,:,k,itr] - humid_grid[:,:,k-1,itr])
    end
    
    arbitrary_array[:,:,2:3] .= zero(NF)
    
    for k in 1:nlev
        humid_tend[:,:,k,itr] = humid_tend + humid_grid[:,:,k,itr].*div_grid[:,:,k]
                - (arbitrary_array[:,:,k+1] + arbitrary_array[:,:,k])*dhsr[k]
        end

end


