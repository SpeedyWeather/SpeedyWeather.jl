
"""
Compute the grid point and spectral tendencies, including an implicit correction for the spectral tendencies. 
"""
function get_tendencies!(
                         #---Diagnostic Variables, IN
                         vor::AbstractArray{NF,4}, #Vorticity.   IN
                         div::AbstractArray{NF,4}, #Divergence.  IN
                         t::AbstractArray{NF,4},   #Temperature. IN
                         ps::AbstractArray{NF,4},  #Surface pressure (logarithm). IN
                         tr::AbstractArray{NF,4},  #Tracers. IN

                         
                        #---Prognositic Variables, IN/OUT
                        phi ::AbstractArray{NF,3} #Atmos. spectral geopotential. IN/OUT
                        phi_surface::AbstractArray{NF,2} #Surface geopotential. IN(?)

                         #Diagnostic variable tendencies
                         vor_tend::AbstractArray{NF,3}, #Vorticity tendency. OUT
                         div_tend::AbstractArray{NF,3}, #Divergence tendency. OUT
                         t_tend::AbstractArray{NF,3},   #Temperature tendency. OUT
                         ps_tend::AbstractArray{NF,3},  #Surface pressure tendency. OUT
                         tr_tend::AbstractArray{NF,3},  #Tracers tendency. OUT
                         M:: ModelSetup #Struct that holds all constants. IN
                         ) where {NF<:AbstractFloat}



    # =========================================================================
    # Computation of grid-point tendencies (converted to spectral at the end of
    # grtend)
    # =========================================================================

    get_grid_point_tendencies!(vor,div,t,ps,tr,
                               phi, phi_surface, 
                               vor_tend, div_tend, t_tend, ps_tend,tr_tend,
                               M)

    # =========================================================================
    # Computation of spectral tendencies
    # =========================================================================

    if α < 0.5 #Coefficient for semi-implicit computations. Previously if alpha = 0?
        get_spectral_tendencies!(D_tend, Tₐ_tend, pₛ_tend, j2)
    else
        get_spectral_tendencies!(D_tend, Tₐ_tend, pₛ_tend, 1)

        # Implicit correction
        implicit_terms!(D_tend, Tₐ_tend, pₛ_tend)
    end
end



"""
Compute the grid point tendencies. These are composed of the dynamics and physical parameterisations. 
"""
function get_grid_point_tendencies!(vor::AbstractArray{NF,4}, #Vorticity.   IN
                                    div::AbstractArray{NF,4}, #Divergence.  IN
                                    t::AbstractArray{NF,4},   #Temperature. IN
                                    ps::AbstractArray{NF,4},  #Surface pressure (logarithm). IN
                                    tr::AbstractArray{NF,4},  #Tracers. IN
                                    phi ::AbstractArray{NF,3} #Atmos. spectral geopotential. IN/OUT
                                    phi_surface::AbstractArray{NF,2} #Surface geopotential. IN(?)
                                    vor_tend::AbstractArray{NF,3}, #Vorticity tendency. OUT
                                    div_tend::AbstractArray{NF,3}, #Divergence tendency. OUT
                                    t_tend::AbstractArray{NF,3},   #Temperature tendency. OUT
                                    ps_tend::AbstractArray{NF,3},  #Surface pressure tendency. OUT
                                    tr_tend::AbstractArray{NF,3},  #Tracers tendency. OUT
                                    M:: ModelSetup #Struct that holds all constants. IN
                                    ) where {NF<:AbstractFloat}


    ### Declare/Allocate local variables. Can we pre-allocate here?

    #Unpack any constants you might need
    @unpack nlon, nlat, nlev = M #

    #Grid point fields
    vor_grid = Array{NF,3}(undef,nlon, nlat, nlev) #Gridpoint field of vorticity
    div_grid = similar(vor_grid)                   #Gridpoint field of divergence
    t_grid = similar(vor_grid)                     #Gridpoint field of temperature
    ps_grid=Array{NF,2}(undef,nlon, nlat)          #Gridpoint field of surface pressure logarithm
    tr_grid = similar(vor_grid)                    #Gridpoint field of tracers
    u_grid = similar(vor_grid)                     #Gridpoint field of zonal velocity
    v_grid = similar(vor_grid)                     #Gridpoint field of meridional velocity
    q_grid = similar(vor_grid)                     #Gridpoint field of specific_humidity
    phi_grid = similar(vor_grid)                   #Gridpoint field of geopotential






    #1. Convert spectral prognostic variables to grid point space
    get_grid_point_fields(vor,div,t,ps,tr
                          phi,phi_surface 
                          vor_grid,div_grid,t_grid,ps_grid, tr_grid, u_grid,v_grid,q_grid,phi_grid,
                          M)

   # 2. Parameterised physics tendencies. 
   #COMMENT: NEEDS DEFINING in parameterisation_tendencies.jl. There is a lot of physics here
   parametrization_tendencies(u_tend, v_tend, t_tend, q_tend)
    #utend: u-wind tendency (gp)
    #vtend: v-wind tendency (gp)
    #ttend: temp. tendency (gp)
    #qtend: spec. hum. tendency (gp)


   #3. Dynamics tendencies
   dynamics_tendencies(vor,div,t,ps,tr, #Prognostic variables
                       vor_grid,div_grid,t_grid,ps_grid, tr_grid, u_grid,v_grid,q_grid,phi_grid, #IN. Grid point fields from 1.
                       u_tend, v_tend, t_tend, q_tend, #IN. Parameterised physics tendencies from 2.
                       vor_tend, div_tend,t_tend,ps_tend,tr_tend # OUT. Spectral tendencies,
                       M
                       )


end


"""
Compute grid point fields of the spectral prognostic tendencies
"""
function get_grid_point_fields(vor::AbstractArray{NF,4}, #Vorticity. IN
                               div::AbstractArray{NF,4}, #Divergence. IN
                               t::AbstractArray{NF,4},   #Temperature IN
                               ps::AbstractArray{NF,4},  #Surface pressure (logarithm). IN
                               tr::AbstractArray{NF,4},  #Tracers. IN
                               phi ::AbstractArray{NF,3}       #Atmos. spectral geopotential# #IN/OUT
                               phi_surface::AbstractArray{NF,2} #Surface geopotential. IN(?)
                               vor_grid::AbstractArray{NF,3}, #Gridpoint field of vorticity. OUT
                               div_grid::AbstractArray{NF,3}, #Gridpoint field of vorticity. OUT
                               t_grid::AbstractArray{NF,3},   #Gridpoint field of vorticity. OUT
                               ps_grid::AbstractArray{NF,2},  #Gridpoint field of vorticity. OUT
                               tr_grid::AbstractArray{NF,3},  #Gridpoint field of vorticity. OUT
                               u_grid::AbstractArray{NF,3},   #Gridpoint field of vorticity. OUT
                               v_grid::AbstractArray{NF,3},   #Gridpoint field of vorticity. OUT 
                               q_grid::AbstractArray{NF,3},   #Gridpoint field of vorticity. OUT
                               phi_grid::AbstractArray{NF,3}, #Gridpoint field of vorticity. OUT
                               M:: ModelSetup) # Struct that holds all constants.


    #Unpack any constants you might need
    @unpack cp, j2, n_trace,nlev,nlat,nlon,coriol = M #


    #1. Compute grid-point fields
    #1.1 Update geopotential in spectral space
    geopotential!(phi, phi_surface, t)
   
   for k in 1:nlev

        
        #Calculate u and v in grid-point space from vorticity and
        #divergence in spectral space
        uvspec!(vor[:,:,k,j2], div[:,:,k,j2], u_grid[:,:,k],v_grid[:,:,k]
    
        #Temperature in grid point space
        convert_to_grid(t[:,:,k,j2], t_grid[:,:,k])  

        #Normalise geopotential by cp to avoid overflows in physics
        convert_to_grid(phi[:,:,k]*(1/cp), phi_grid[:,:,k]) 

        #Vorticity and divergence grids. I think there is a conversion here to 'per hour' to reduce underflows 
        vor_grid[:,:,k] = convert_to_grid(vor[:,:,k,j2])
        div_grid[:,:,k] = convert_to_grid(div[:,:,k,j2])

        for itr in 1:n_trace
            tr_grid[:,:,k,itr] = convert_to_grid(tr[:,:,k,j2,itr])
        end

        for j in 1:nlat
            for i in 1:nlon
                vor_grid[i,j,k] += coriol[j] 
            end
        end

   end



    #Truncate variables where the spectral transform is still done at double
    #precision
    u_grid = u_grid / 3600.0
    v_grid = v_grid / 3600.0

    #Surface pressure spectral transform to grid
    convert_to_grid(ps(:,:,j1), ps_grid)


    #Don't transform the two stratospheric levels where humidity is set to zero
    # because it leads to overflows 
    for k in 3:nlev
        convert_to_grid(tr(:,:,k,j1,1), q_grid(:.:,k), )
    end 

end






"""
Compute non-linear tendencies in grid-point space from dynamics and add to physics tendencies. Convert total
gridpoint tendencies to spectral tendencies
"""
function dynamics_tendencies(
                             vor,div,t,ps,tr, #Prognostics
                             vor_grid::AbstractArray{NF,3},
                             div_grid::AbstractArray{NF,3},
                             t_grid::AbstractArray{NF,3},
                             tr_grid:::AbstractArray{NF,3}, 
                             u_grid::AbstractArray{NF,3},
                             v_grid::AbstractArray{NF,3},
                             q_grid::AbstractArray{NF,3}, 
                             ps_grid::AbstractArray{NF,2}
                             phi_grid::AbstractArray{NF,3},
                             u_tend, v_tend, t_tend, q_tend, #Need to define after parametrization_tendencies() setup
                             vor_tend::AbstractArray{NF,3},
                             div_tend::AbstractArray{NF,3},
                             t_tend::AbstractArray{NF,3},
                             ps_tend::AbstractArray{NF,3},
                             tr_tend::AbstractArray{NF,3},
                             M
                             ) where {NF<:AbstractFloat}

    #Get any constants you might need
    @unpack tref = M #

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
                                         M)
 
    # 2. Compute "vertical" velocity
    vertical_velocity_tendency(u_grid, v_grid, div_grid,u_mean,px,v_mean,d_mean, #IN
                                         puv,sigma_tend,sigma_m #OUT
                                         M)
   

    # 3. Subtract part of temperature field that is used as reference for implicit terms
    t_grid_anom = Array{NF,3}(undef,nlon, nlat,nlev) 
    
    for k in 1:nlev
        t_grid_anom[:,:,k] = t_grid[:,:,k] .- tref[k] 
    end

    # 4. Zonal wind tendency
    zonal_wind_tendency(u_grid,v_grid,vor_grid,t_grid_anom,px,sigma_tend, # IN
                                   u_tend, #OUT
                                   M)


    # 5. Meridional wind tendency
    meridional_wind_tendency(u_grid,v_grid,vor_grid,t_grid_anom,px,sigma_tend, # IN
                                        v_tend, #OUT
                                        M)

    # 6. Temperature tendency
    temperature_tendency(sigma_tend,t_grid_anom, sigma_m, ,puv, div_grid, t_grid, #IN
                                    t_tend, #IN/OUT
                                    M)

    
    # 7. Tracer tendency
    tracer_tendency(sigma_tend, tr_grid,div_grid,#IN
                    tr_tend, #IN/OUT
                    M)



    # =========================================================================
    # Convert tendencies to spectral space
    # =========================================================================


    for k in 1:nlev
        #  Convert u and v tendencies to vor and div spectral tendencies
        #  vdspec takes a grid u and a grid v and converts them to spectral vor and div

        vdspec(utend[:,:,k], vtend[:,:,k], vor_tend[:,:,k], div_tend[:,:,k], true) #COMMENT: is vdspec working OK?

        #Divergence tendency
        #add -laplacian(0.5*(u**2+v**2)) to divergence tendency
        div_tend[:,:,k] = div_tend[:,:,k] - laplacian(convert_to_spectral(half*(u_grid[:,:,k]**two + v_grid[:,:,k]**two))) #COMMENT: need to define ints/floats in consistent way

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
Compute the spectral tendencies.  
"""
function get_spectral_tendencies!(
                                 div_tend::AbstractArray{NF,3}, 
                                 t_tend::AbstractArray{NF,3}, 
                                 ps_tend::AbstractArray{NF,3}
                                 
                                 div #Prognostic variable, divergence
                                 
                                 
                                 
                                 ) where {NF<:AbstractFloat}


    
   #Get any constants you might need
   @unpack mx, nx, nlev, dhs,j2,tref,tref2,tref3 = M #


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
    get_geopotential(Tₐ[:,:,:,j2])
    for k in 1:nlev
        D_tend[:,:,k] -= ∇²(ϕ[:,:,k] + R*tref[k]*pₛ[:,:,j2])
    end
end