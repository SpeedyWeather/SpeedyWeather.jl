"""
Compute the spectral tendency of the surface pressure logarithm
"""
function surface_pressure_tendency!(Prog::PrognosticVariables{NF}, # Prognostic variables
                                   Diag::PrognosticVariables{NF}, # Diagnostic variables
                                   l2::Int,                       # leapfrog index 2 (time step used for tendencies)
                                   C::Constants{NF}
                                   ) where {NF<:AbstractFloat}


    @unpack logp0= Prog 
    @unpack logp0_tend = Diag.Tendencies

    @unpack u_grid,v_grid,div_grid = Diag.gridvars
    @unpack u_mean,v_mean,d_mean,dumc,px,py = Diag.miscvars
    @unpack nlev, dhs = C 

    #Initialise mean fields. 
    u_mean[:,:] = 0.0
    v_mean[:,:] = 0.0
    d_mean[:,:] = 0.0 #d-mean is calculated here but actually used in vertical_velocity_tendency!()

    #..and calculate values
    for k in 1:nlev
        u_mean += u_grid[:,:,k]*dhs[k] 
        v_mean += v_grid[:,:,k]*dhs[k]
        d_mean += div_grid[:,:,k]*dhs[k]
    end

    #Now use the mean fields
    grad!(logp0[:,:,l2], dumc[:,:,2], dumc[:,:,3])
    px = convert_to_grid(dumc[:,:,2]*3600, scale=true)
    py = convert_to_grid(dumc[:,:,3]*3600, scale=true) #3600 factor from Paxton/Chantry. I think this is to correct for the underflow rescaling earlier

    logp0_tend = convert_to_spectral(-u_mean.*px - v_mean.*py)
    logp0_tend[1,1] = Complex{NF}(zero)

end

"""
Compute the spectral tendency of the "vertical" velocity
"""
function vertical_velocity_tendency!(Prog::PrognosticVariables{NF}, # Prognostic variables
                                     Diag::PrognosticVariables{NF}, # Diagnostic variables
                                     C::Constants{NF}
                                     ) where {NF<:AbstractFloat}

    @unpack u_grid,v_grid,div_grid = Diag.gridvars
    @unpack u_mean,v_mean,d_mean,sigma_tend,sigma_m,puv,px,py = Diag.miscvars
    @unpack nlev,dhs = C 

    #Initialise
    sigma_tend[:,:,:] = 0.0
    sigma_m[:,:,:] = 0.0


    for k in 1:nlev
        puv[:,:,k] = (u_grid[:,:,k] - u_mean)*px + (v_grid[:,:,k] - v_mean)*py
    end

    for k in 1:nlev
        sigma_tend[:,:,k+1] = sigma_tend[:,:,k] - dhs[k]*(puv[:,:,k] + div_grid[:,:,k] - d_mean)
        sigma_m[:,:,k+1] = sigma_m[:,:,k] - dhs[k]*puv[:,:,k]
    end

end



"""
Compute the temperature anomaly in grid point space
"""
function temperature_grid_anomaly!(Diag::PrognosticVariables{NF}, # Diagnostic variables
                                   C::Constants{NF}
                                   ) where {NF<:AbstractFloat}

    @unpack Tabs_grid,Tabs_grid_anomaly = Diag.gridvars
    @unpack nlev,tref = C  #Is tref a constant?

    for k in 1:nlev
        Tabs_grid_anomaly[:,:,k] = Tabs_grid[:,:,k] .- tref[k] 
    end

end





"""
Compute the spectral tendency of the zonal wind
"""
function zonal_wind_tendency!(Diag::PrognosticVariables{NF}, # Diagnostic variables
                              C::Constants{NF}
                              )where {NF<:AbstractFloat}
    
    @unpack u_tend = Diag.Tendencies
    @unpack u_grid,v_grid,vor_grid,Tabs_grid_anomaly= Diag.gridvars
    @unpack sigma_tend,px,py,temp = Diag.miscvars
    @unpack nlev,rgas,dhsr = C 


    #Update px,py
    px = rgas*px
    py = rgas*py

    #Initialise temp
    temp[:,:,:] =0.0

  
    for k in 2:nlev
        temp[:,:,k] = sigma_tend[:,:,k].*(u_grid[:,:,k] - u_grid[:,:,k-1])
    end

    for k in 1:nlev
        u_tend[:,:,k] = u_tend[:,:,k] 
                        + v_grid[:,:,k].*vor_grid[:,:,k] 
                        - Tabs_grid_anomaly[:,:,k]*px
                         - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
    end

end



"""
Compute the spectral tendency of the meridional wind 
"""
function meridional_wind_tendency!(Diag::PrognosticVariables{NF}, # Diagnostic variables
                                   C::Constants{NF}
                                  )where {NF<:AbstractFloat}

    @unpack v_tend = Diag.Tendencies
    @unpack vor_grid,u_grid,v_grid,Tabs_grid_anomaly =Diag.gridvars
    @unpack sigma_tend,temp,px,py = Diag.miscvars
    @unpack nlev,dhsr= C 


    for k in 2:nlev
        temp[:,:,k] = sigma_tend[:,:,k].*(v_grid[:,:,k] - v_grid[:,:,k-1])
    end
                            
    for k in 1:nlev
        v_tend[:,:,k] = v_tend[:,:,k] 
                        -u_grid[:,:,k].*vor_grid[:,:,k] 
                        - Tabs_grid_anomaly[:,:,k]*py
                        - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
    end

end





"""
Compute the spectral temperature tendency
"""
function temperature_tendency!(Diag::PrognosticVariables{NF}, # Diagnostic variables
                               C::Constants{NF}
                               )where {NF<:AbstractFloat}
    
    @unpack Tabs_tend = Diag.Tendencies
    @unpack div_grid,Tabs_grid_anomaly =Diag.gridvars
    @unpack d_mean,temp,sigma_tend,sigma_m,puv= Diag.miscvars
    @unpack tref,tref3,nlev,dhsr,fsgr,akap = C




    for k in 2:nlev
        temp[:,:,k] = sigma_tend[:,:,k].*(Tabs_grid_anomaly[:,:,k] - Tabs_grid_anomaly[:,:,k-1])
                    + sigma_m[:,:,k].*(tref[k] - tref[k-1])
    end

    for k in 1:nlev
        Tabs_tend[:,:,k] = Tabs_tend[:,:,k]
                        + Tabs_grid_anomaly[:,:,k].*div_grid[:,:,k]
                        - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
                        + fsgr[k]*Tabs_grid_anomaly[:,:,k].*(sigma_tend[:,:,k+1] + sigma_tend[:,:,k])
                        + tref3[k]*(sigma_m[:,:,k+1] + sigma_m[:,:,k])
                        + akap*(t_grid[:,:,k].*puv[:,:,k] - Tabs_grid_anomaly[:,:,k].*d_mean)
    end 


end




"""
Compute the tracer tendency
"""
function tracer_tendency!(Diag::PrognosticVariables{NF}, # Diagnostic variables
                          C::Constants{NF}
                          )where {NF<:AbstractFloat}
    

    @unpack div_grid,tr_grid = Diag.gridvars
    @unpack temp,sigma_tend= Diag.miscvars
    @unpack n_trace,dhsr, nlev = C


    for itr in 1:n_trace
        for k in 2:nlev
            temp[:,:,k] = sigma_tend[:,:,k].*(tr_grid[:,:,k,itr] - tr_grid[:,:,k-1,itr])
        end
    
        temp[:,:,2:3] .= zero
    
        for k in 1:nlev
            tr_tend[:,:,k,itr] = tr_tend + tr_grid[:,:,k,itr].*div_grid[:,:,k]
                - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
        end
    end

end


