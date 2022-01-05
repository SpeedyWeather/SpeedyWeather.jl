
# """Struct holding the spectral prognostic variables in grid point space"""
# struct PrognosticVariablesGridded{NF<:AbstractFloat}
#     vor_grid = Array{NF,3}(undef,nlon, nlat, nlev) #Gridpoint field of vorticity
#     div_grid = similar(vor_grid)                   #Gridpoint field of divergence
#     t_grid = similar(vor_grid)                     #Gridpoint field of temperature
#     ps_grid=Array{NF,2}(undef,nlon, nlat)          #Gridpoint field of surface pressure logarithm
#     tr_grid = similar(vor_grid)                    #Gridpoint field of tracers
#     u_grid = similar(vor_grid)                     #Gridpoint field of zonal velocity
#     v_grid = similar(vor_grid)                     #Gridpoint field of meridional velocity
#     q_grid = similar(vor_grid)                     #Gridpoint field of specific_humidity
#     phi_grid = similar(vor_grid)                   #Gridpoint field of geopotential
# end


# """Struct holding the spectral prognostic variables tendencies"""
# struct PrognosticVariableTendencies{NF<:AbstractFloat}
#     vor_tend = Array{NF,3}(undef,mx, nx, nlev) #spectral tendency of vorticity
#     div_tend = similar(vor_tend)               #spectral tendency of divergence
#     t_tend   = similar(vor_tend)               #spectral tendency of temperature
#     ps_tend  = similar(vor_tend)               #spectral tendency of surface pressure logarithm
#     tr_tend  = similar(vor_tend)               #spectral tendency of tracers
# end



"""
Compute the grid point and spectral tendencies, including an implicit correction for the spectral tendencies. 
"""
function get_tendencies!(Prog::PrognosticVariables{NF}, # Prognostic variables
                         Diag::PrognosticVariables{NF}, # Diagnostic variables
                         C::Constants{NF}
                        ) where {NF<:AbstractFloat}

    @unpack alpha = C

    # =========================================================================
    # Computation of grid-point tendencies (converted to spectral at the end of
    # grtend) Diag.Tendencies.GridPoint
    # =========================================================================

    get_grid_point_tendencies!(Prog,Diag,C)

    # =========================================================================
    # Computation of spectral tendencies 
    # =========================================================================
    @unpack div_tend,t_tend,ps_tend = Diag.Tendencies #Get some of the tendencies calculated above
    if alpha < 0.5 #Coefficient for semi-implicit computations. Previously if alpha = 0?
        get_spectral_tendencies!(div_tend, t_tend, ps_tend,C)
    else
        get_spectral_tendencies!(div_tend, t_tend, ps_tend,C)

        # Implicit correction
        implicit_terms!(div_tend, t_tend, ps_tend) #TK: no edits have been made (yet!) to this implicit function
    end

end



"""
Compute the grid point tendencies. These are composed of the dynamics and physical parameterisations. 
"""
function get_grid_point_tendencies!(Prog::PrognosticVariables{NF}, # Prognostic variables
                                    Diag::PrognosticVariables{NF}, # Diagnostic variables
                                    C::Constants{NF}
                                    ) where {NF<:AbstractFloat}
    
    #1. Convert spectral prognostic variables to grid point space
    get_grid_point_fields!(Prog,C) #Updates Prog.Gridded
 
    # 2. Parameterised physics tendencies. 
    #Needs to be defined in parameterisation_tendencies.jl. There is a lot of physics here, so would be best as separate, self-contained PR
    parametrization_tendencies!(Prog,Diag,C)


    #3. Dynamics tendencies
    dynamics_tendencies!(Prog,Diag,C) #Takes Prog.Grid and Diag.ParameterisedTendencies and calculates Diag.Tendencies in GP space



end


"""
Compute grid point fields of the spectral prognostic tendencies
"""
function get_grid_point_fields(Prog::PrognosticVariables{NF}, # Prognostic variables
                               C::Constants{NF}
                               ) where {NF<:AbstractFloat}


    #Unpack spectral prognostic fields
    @unpack vor,div,t,ps,tr = Prog

    #Unpack the gridded counterparts. This is what we will be calculating in this subroutine
    @unpack vor_grid,div_grid,t_grid,ps_grid,tr_grid = Prog.Gridded

    #Unpack constants
    @unpack cp, j2, n_trace,nlev,nlat,nlon,coriol = C


    #1. Compute grid-point fields
    #1.1 Update geopotential in spectral space
    geopotential!(phi, phi_surface, t)

   
   for k in 1:nlev

    
        convert_to_grid!(vor[:,:,k,j2],vor_grid[:,:,k]) # vorticity 
        convert_to_grid!(div[:,:,k,j2],div_grid[:,:,k]) # divergence
        convert_to_grid!(t[:,:,k,j2], t_grid[:,:,k])    # temperature

        for itr in 1:n_trace #tracers
            convert_to_grid!(tr[:,:,k,j2,itr],tr_grid[:,:,k,itr] )
        end


        #Calculate zonal velocity u and meridional velocity v in grid-point space,
        #from vorticity and divergence in spectral space
        uvspec!(vor[:,:,k,j2], div[:,:,k,j2], u_grid[:,:,k],v_grid[:,:,k])

        #Normalise geopotential by cp to avoid overflows in physics
        convert_to_grid(phi[:,:,k]*(1/cp), phi_grid[:,:,k]) 

       

        #Correct vorticity grid point field
        for j in 1:nlat
            for i in 1:nlon
                vor_grid[i,j,k] += coriol[j] 
            end
        end

   end



    #Truncate variables where the spectral transform is still done at double
    #precision. Convert to per second
    u_grid = u_grid / 3600.0
    v_grid = v_grid / 3600.0

    #Surface pressure spectral transform to grid
    convert_to_grid!(ps[:,:,j1], ps_grid)


    #Don't transform the two stratospheric levels where humidity is set to zero
    # because it leads to overflows 
    for k in 3:nlev
        convert_to_grid!(tr[:,:,k,j1,1], tr_grid[:,:,k])
    end 


end


"""
Compute non-linear tendencies in grid-point space from dynamics and add to physics tendencies. Convert total
gridpoint tendencies to spectral tendencies.
"""
function dynamics_tendencies(Prog::PrognosticVariables{NF}, # Prognostic variables
                             Diag::PrognosticVariables{NF}, # Diagnostic variables
                             C::Constants{NF}
                             ) where {NF<:AbstractFloat}
    
    #Unpack spectral prognostic fields
    @unpack vor,div,t,ps,tr = Prog

    #Unpack the gridded counterparts
    @unpack vor_grid,div_grid,t_grid,ps_grid,tr_grid = Prog.Gridded

    #Unpack the parameterised physics tendencies
    @unpack u_tend, v_tend, t_tend, q_tend = Diag.ParameterisedTendencies
    
    #Unpack the spectral prognostic tendencies. These are what we will be calculating in this subroutine
    @unpack vor_tend, div_tend, t_tend, ps_tend, tr_tend  = Diag.Tendencies

    #Unpack constants
    @unpack cp, j2, n_trace,nlev,nlat,nlon,coriol = C


    #Declare some empty arrays which are shared between subroutines

    #Filled in 1. dynamics_tendencies_surface_pressure()
    u_mean = Array{NF,2}(undef,nlon, nlat) 
    v_mean = similar(u_mean)
    d_mean = similar(u_mean)
    dumc = Array{NF,3}(mx,nx,3) 
    px = Array{NF,2}(udef,ix,il)
    py = similar(px)

    #Filled in 2. dynamics_tendencies_surface_pressure()
    puv = Array{NF,3}(undef,nlon, nlat,nlev) 
    sigma_m = Array{NF,3}(undef,nlon, nlat,nlev+1) 
    sigma_tend = Array{NF,3}(undef,nlon, nlat,nlev+1) 


    #1. Compute tendency of log(surface pressure)
    surface_pressure_tendency(u_grid,v_grid,div_grid,ps, #IN
                              ps_tend,u_mean,v_mean,d_mean,dumc,px,py,#OUT. Tendencies + other terms used later for other tendency calulcations
                              C)
 
    # 2. Compute "vertical" velocity
    vertical_velocity_tendency(u_grid, v_grid, div_grid,u_mean,px,v_mean,d_mean, #IN
                               puv,sigma_tend,sigma_m, #OUT
                               C)
   

    # 3. Subtract part of temperature field that is used as reference for implicit terms
    t_grid_anom = Array{NF,3}(undef,nlon, nlat,nlev) 
    
    for k in 1:nlev
        t_grid_anom[:,:,k] = t_grid[:,:,k] .- tref[k] 
    end

    # 4. Zonal wind tendency
    zonal_wind_tendency(u_grid,v_grid,vor_grid,t_grid_anom,px,sigma_tend, # IN
                        u_tend, #OUT
                        C)


    # 5. Meridional wind tendency
    meridional_wind_tendency(u_grid,v_grid,vor_grid,t_grid_anom,px,sigma_tend, # IN
                             v_tend, #OUT
                             C)

    # 6. Temperature tendency
    temperature_tendency(sigma_tend,t_grid_anom, sigma_m, puv, div_grid, t_grid, #IN
                         t_tend, #IN/OUT
                         C)

    
    # 7. Tracer tendency
    tracer_tendency(sigma_tend, tr_grid,div_grid,#IN
                    tr_tend, #IN/OUT
                    C)



    # =========================================================================
    # Convert tendencies to spectral space
    # =========================================================================


    for k in 1:nlev
        #  Convert u and v tendencies to vor and div spectral tendencies
        #  vdspec takes a grid u and a grid v and converts them to spectral vor and div

        vdspec(utend[:,:,k], vtend[:,:,k], vor_tend[:,:,k], div_tend[:,:,k], true) 

        #Divergence tendency
        #add -laplacian(0.5*(u**2+v**2)) to divergence tendency
        div_tend[:,:,k] = div_tend[:,:,k] - laplacian(convert_to_spectral(half*(u_grid[:,:,k]**two + v_grid[:,:,k]**two))) 

        #Temperature tendency
        #and add div(vT) to spectral t tendency
        vdspec(-u_grid[:,:,k]*t_grid_anom[:,:,k], 
               -v_grid[:,:,k]*t_grid_anom[:,:,k], 
               dumc(:,:,1), t_tend[:,:,k], 
               true)
    
        t_tend(:,:,k) = t_tend(:,:,k) + convert_to_spectral(t_grid[:,:,k])

        #tracer tendency
        do itr = 1, n_trace
            call vdspec(-u_grid[:,:,k]*tr_grid(:,:,k,itr), 
                        -v_grid[:,:,k]*tr_grid[:,:,k,itr], 
                        dumc[:,:,1], tr_tend[:,:,k,itr], 
                        true)
        tr_tend[:,:,k,itr] = tr_tend[:,:,k,itr] + convert_to_spectral(tr_grid[:,:,k,itr])
    end do

    end do




"""
Compute the spectral tendencies of divergence, temperature
and log_surf.pressure 
"""
function get_spectral_tendencies!(div_tend::AbstractArray{NF,3}, #spectral tendency of divergence
                                 t_tend::AbstractArray{NF,3}, #spectral tendency of temperature
                                 ps_tend::AbstractArray{NF,3},#spectral tendency of surface pressure logarithm
                                 C::Constants{NF}
                                 ) where {NF<:AbstractFloat}


    
   #Get any constants you might need
   @unpack mx, nx, nlev, dhs,j2,tref,tref2,tref3,rgas = C


   #Declare local arrays   
   dmeanc = Array{NF,2}(mx,nx) 
   sigma_tend_c = Array{NF,3}(mx, nx, nlev+1)
   dumk = Array{NF,3}(mx, nx, nlev+1)


    # Vertical mean divergence... 
    for k in 1:nlev
        dmeanc = dmeanc + div[:,:,k,j2]*dhs[k]
    end

    #...and and pressure tendency
    ps_tend = ps_tend - dmeanc
    ps_tend[1,1] = Complex{RealType}(zero)

    # Sigma-dot "velocity" and temperature tendency
    for k in 1:nlev - 1
        sigma_tend_c[:,:,k+1] = sigma_tend_c[:,:,k] - dhs[k]*(div[:,:,k,j2] - dmeanc)
    end

    for k in 2:nlev
        dumk[:,:,k] = sigma_tend_c[:,:,k]*(tref[k] - tref[k-1])
    end

    for k in 1:nlev
        t_tend[:,:,k] -= (dumk[:,:,k+1] + dumk[:,:,k])*dhsr[k]
            + tref3[k]*(sigma_tend_c[:,:,k+1] + sigma_tend_c[:,:,k]) - tref2[k]*dmeanc
    end

    # Geopotential and divergence tendency
    geopotential!(phi, phi_surface, t)
    for k in 1:nlev
        div_tend[:,:,k] -= ∇²(phi[:,:,k] + rgas*tref[k]*ps[:,:,j2])
    end
end