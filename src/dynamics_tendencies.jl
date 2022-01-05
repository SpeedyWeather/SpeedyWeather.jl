"""
Compute the spectral tendency of the surface pressure logarithm
"""
function surface_pressure_tendency(u_grid::AbstractArray{NF,3},   #IN
                                   v_grid::AbstractArray{NF,3},   #IN
                                   div_grid::AbstractArray{NF,3}, #IN
                                   ps::AbstractArray{NF,4},       #IN
                                   ps_tend::AbstractArray{NF,3},  #OUT
                                   u_mean::AbstractArray{NF,2},   #OUT
                                   v_mean::AbstractArray{NF,2},   #OUT
                                   d_mean::AbstractArray{NF,2},   #OUT
                                   dumc::AbstractArray{NF,3},     #OUT
                                   px::AbstractArray{NF,2},       #OUT
                                   py::AbstractArray{NF,2},       #OUT
                                   C::Constants{NF}               #IN
                                   )where {NF<:AbstractFloat}


    @unpack nlev, dhs = C 

    #Now do some calculation
    for k in 1:nlev
        u_mean += u_grid[:,:,k]*dhs[k] 
        v_mean += v_grid[:,:,k]*dhs[k]
        d_mean += div_grid[:,:,k]*dhs[k]
    end

    
    grad!(ps[:,:,j2], @view(dumc[:,:,1]), @view(dumc[:,:,2]))
    px = convert_to_grid(dumc[:,:,1], scale=true)
    py = convert_to_grid(dumc[:,:,2], scale=true)

    ps_tend = convert_to_spectral(-umean.*px - vmean.*py)
    ps_tend[1,1] = Complex{NF}(zero)

end

"""
Compute the spectral tendency of the "vertical" velocity
"""
function vertical_velocity_tendency(
                                   u_grid::AbstractArray{NF,3}, #IN
                                   v_grid::AbstractArray{NF,3}, #IN
                                   div_grid::AbstractArray{NF,3}, #IN
                                   u_mean::AbstractArray{NF,2}, #IN
                                   px::AbstractArray{NF,2}, #IN
                                   v_mean::AbstractArray{NF,2}, #IN
                                   py::AbstractArray{NF,2},#IN
                                   d_mean::AbstractArray{NF,2},#IN,
                                   puv::AbstractArray{NF,3},#OUT,
                                   sigma_tend::AbstractArray{NF,3},#OUT,
                                   sigma_m::AbstractArray{NF,3},#OUT,
                                   C::Constants{NF} 
                                   )where {NF<:AbstractFloat}

    #Get constants
    @unpack nlev,dhs = C #

    #Declare empty arrays
    
    for k in 1:nlev
        puv[:,:,k] = (u_grid[:,:,k] - u_mean).*px + (v_grid[:,:,k] - v_mean).*py
    end

    for k in 1:nlev
        sigma_tend[:,:,k+1] = sigma_tend[:,:,k] - dhs[k]*(puv[:,:,k] + div_grid[:,:,k] - d_mean)
        sigma_m[:,:,k+1] = sigma_m[:,:,k] - dhs[k]*puv[:,:,k]
    end

end


"""
Compute the spectral tendency of the zonal wind
"""
function zonal_wind_tendency(
                                       u_grid::AbstractArray{NF,3}, #IN
                                       v_grid::AbstractArray{NF,3}, #IN
                                       vor_grid::AbstractArray{NF,3}, #IN
                                       t_grid_anom::AbstractArray{NF,3}, #IN
                                       px::AbstractArray{NF,2}, #IN
                                       sigma_tend::AbstractArray{NF,3},#IN,
                                       u_tend::AbstractArray{NF,3},#OUT
                                       M
                                       )where {NF<:AbstractFloat}

    @unpack nlon, nlat,nlev,dhsr = M 

    #Temporary array
    temp = Array{NF,3}(undef,nlon, nlat,nlev+1) 

    for k in 2:nlev
        temp[:,:,k] = sigma_tend[:,:,k].*(u_grid[:,:,k] - u_grid[:,:,k-1])
    end

    for k in 1:nlev
        u_tend[:,:,k] = u_tend[:,:,k] + v_grid[:,:,k].*vor_grid[:,:,k] - t_grid_anom[:,:,k]*px
            - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
    end

end


"""
Compute the spectral tendency of the meridional wind 
"""
function meridional_wind_tendency(
                                       u_grid::AbstractArray{NF,3}, #IN
                                       v_grid::AbstractArray{NF,3}, #IN
                                       vor_grid::AbstractArray{NF,3}, #IN
                                       t_grid_anom::AbstractArray{NF,3}, #IN
                                       px::AbstractArray{NF,2}, #IN
                                       sigma_tend::AbstractArray{NF,3},#IN,
                                       v_tend::AbstractArray{NF,3},#OUT
                                       M
                                       )where {NF<:AbstractFloat}
    @unpack nlon, nlat,nlev,dhsr = M 

    #Temporary array
    temp = Array{NF,3}(undef,nlon, nlat,nlev+1) 
    for k in 2:nlev
        temp[:,:,k] = sigma_tend[:,:,k].*(v_grid[:,:,k] - v_grid[:,:,k-1])
      end
                                
    for k in 1:nlev
        v_tend[:,:,k] = v_tend[:,:,k] -u_grid[:,:,k].*vor_grid[:,:,k] - t_grid_anom[:,:,k]*R.*py
            - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
    end

end


"""
Compute the spectral temperature tendency
"""
function temperature_tendency(
                                        sigma_tend::AbstractArray{NF,3},  #IN,
                                        t_grid_anom::AbstractArray{NF,3}, #IN
                                        sigma_m::AbstractArray{NF,3},     #IN
                                        puv::AbstractArray{NF,3},         #IN
                                        div_grid::AbstractArray{NF,3},    #IN
                                        t_grid::AbstractArray{NF,3},      #IN
                                        t_tend::AbstractArray{NF,3},      #IN/OUT
                                        M
                                        )where {NF<:AbstractFloat}
    
    
    
    @unpack nlon, nlat,nlev,dhsr,tref,tref3,fsgr,akap = M 


    #Temporary array
    temp = Array{NF,3}(undef,nlon, nlat,nlev+1) 


    for k in 2:nlev
        temp[:,:,k] = sigma_tend[:,:,k].*(t_grid_anom[:,:,k] - t_grid_anom[:,:,k-1])
            + sigma_m[:,:,k].*(tref[k] - tref[k-1])
    end

    for k in 1:nlev
        t_tend[:,:,k] = t_tend[:,:,k]+ t_grid_anom[:,:,k].*div_grid[:,:,k]
            - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
            + fsgr[k]*t_grid_anom[:,:,k].*(sigma_tend[:,:,k+1] + sigma_tend[:,:,k])
            + tref3[k]*(sigma_m[:,:,k+1] + sigma_m[:,:,k])
            + akap*(t_grid[:,:,k].*puv[:,:,k] - t_grid_anom[:,:,k].*dmean)
    end 


end




"""
Compute the tracer tendency
"""
function tracer_tendency(
                             sigma_tend::AbstractArray{NF,3},  #IN,
                             tr_grid:::AbstractArray{NF,3}, 
                             div_grid::AbstractArray{NF,3},
                             tr_tend::AbstractArray{NF,3}
                             M
                            )where {NF<:AbstractFloat}
    
    
    
    @unpack nlon, nlat,nlev,n_trace,dhsr = M 

    #Temporary array
    temp = Array{NF,3}(undef,nlon, nlat,nlev+1) 

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


