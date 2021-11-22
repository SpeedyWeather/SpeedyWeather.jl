
"""
Compute the grid point and spectral tendencies, including an implicit correction for the spectral tendencies. 
"""
function get_tendencies!(vor::AbstractArray{NF,4},
                         div::AbstractArray{NF,4}, 
                         t::AbstractArray{NF,4},
                         ps::AbstractArray{NF,4}, 
                         tr::AbstractArray{NF,4},
                         phi::AbstractArray{NF,3},
                         vor_tend::AbstractArray{NF,3}, 
                         div_tend::AbstractArray{NF,3},
                         t_tend::AbstractArray{NF,3},
                         ps_tend::AbstractArray{NF,3},
                         tr_tend::AbstractArray{NF,3},
                         j1::int, #COMMENT: Maybe j1/j2 should be read in from struct?
                         j2::int) where {NF<:AbstractFloat}



    # =========================================================================
    # Computation of grid-point tendencies (converted to spectral at the end of
    # grtend)
    # =========================================================================

    get_grid_point_tendencies!(vor,div,t,ps,tr,phi, 
                               vor_tend, div_tend, t_tend, ps_tend,tr_tend,
                               j1, j2)

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
function get_gridpoint_tendencies!(vor::AbstractArray{NF,4},
                                   div::AbstractArray{NF,4}, 
                                   t::AbstractArray{NF,4},
                                   ps::AbstractArray{NF,4}, 
                                   tr::AbstractArray{NF,4},
                                   phi::AbstractArray{NF,4},
                                   vor_tend::AbstractArray{NF,3}, 
                                   div_tend::AbstractArray{NF,3},
                                   t_tend::AbstractArray{NF,3},
                                   ps_tend::AbstractArray{NF,3},
                                   tr_tend::AbstractArray{NF,3},
                                   j1::int,
                                   j2::int) where {NF<:AbstractFloat}

    ### Declare/Allocate variables. Can we pre-allocate here?

    #Grid point fields
    vor_grid = Array{NF,3}(undef,nlon, nlat, nlev) #Gridpoint field of vorticity
    div_grid = similar(vor_grid)                   #Gridpoint field of divergence
    t_grid = similar(vor_grid)                     #Gridpoint field of temperature
    tr_grid = similar(vor_grid)                    #Gridpoint field of tracers
    u_grid = similar(vor_grid)                     #Gridpoint field of zonal velocity
    v_grid = similar(vor_grid)                     #Gridpoint field of meridional velocity
    q_grid = similar(vor_grid)                     #Gridpoint field of specific_humidity
    ps_grid=Array{NF,2}(undef,nlon, nlat)          #Gridpoint field of surface pressure logarithm
    phi_grid = similar(vor_grid)                   #Gridpoint field of geopotential


    #1. Convert spectral prognostic variables to grid point space
    get_grid_point_fields(vor,div,t,phi,ps,tr,
                          vor_grid,div_grid,t_grid,tr_grid, u_grid,v_grid,q_grid, ps_grid,phi_grid)


   # 2. Parameterised physics tendencies. 
   #COMMENT: NEEDS DEFINING in phypar.jl. There is a lot of physics here
   phypar(u_tend, div_tend, t_tend, q_tend)
    #utend: u-wind tendency (gp)
    #vtend: v-wind tendency (gp)
    #ttend: temp. tendency (gp)
    #qtend: spec. hum. tendency (gp)


   #3. Dynamics tendencies
   dynamics_tendencies(u_tend, div_tend, t_tend, q_tend, #IN, Parameterised physics tendencies
                       vor_grid,div_grid,t_grid,tr_grid, u_grid,v_grid,q_grid, ps_grid,phi_grid, #IN, Grid point fields
                       j2, #IN, timestep parameter
                       vor_tend, div_tend,t_tend,ps_tend,tr_tend,#OUT, tendencies
                       ps #IN, other. Typically defined globally in the fortran version.
                       )


end


"""
Compute grid point fields of the spectral prognostic tendencies
"""
function get_grid_point_fields(vor::AbstractArray{NF,4}, #IN
                               div::AbstractArray{NF,4}, #IN
                               t::AbstractArray{NF,4},   #IN
                               phi::AbstractArray{NF,4}, #IN
                               ps::AbstractArray{NF,4},  #IN
                               tr::AbstractArray{NF,4},  #IN
                               vor_grid::AbstractArray{NF,3},#OUT
                               div_grid::AbstractArray{NF,3},#OUT
                               t_grid::AbstractArray{NF,3},  #OUT
                               tr_grid::AbstractArray{NF,3}, #OUT
                               u_grid::AbstractArray{NF,3},  #OUT
                               v_grid::AbstractArray{NF,3},  #OUT 
                               q_grid::AbstractArray{NF,3},  #OUT 
                               ps_grid::AbstractArray{NF,2}, #OUT
                               phi_grid::AbstractArray{NF,3}) #OUT



    #1. Compute grid-point fields
    #1.1 Update geopotential in spectral space
    #COMMENT: THIS NEEDS TO BE CONFIGURED PROPERLY W.R.T ARGUMENTS.
    geopotential!(phi)

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
                vor_grid[i,j,k] += coriol[j] #COMMENT: coriol will need to be defined, perhaps in a struct?
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
function dynamics_tendencies(u_tend, div_tend, t_tend, q_tend,#COMMENT. Need to declare physics tendencies types once phypar is completed
                             vor_grid::AbstractArray{NF,3},
                             div_grid::AbstractArray{NF,3},
                             t_grid::AbstractArray{NF,3},
                             tr_gri:::AbstractArray{NF,3}, 
                             u_grid::AbstractArray{NF,3},
                             v_grid::AbstractArray{NF,3},
                             q_grid::AbstractArray{NF,3}, 
                             ps_grid::AbstractArray{NF,2}
                             phi_grid::AbstractArray{NF,3},
                             j2 :: ::int
                             vor_tend::AbstractArray{NF,3},
                             div_tend::AbstractArray{NF,3},
                             t_tend::AbstractArray{NF,3},
                             ps_tend::AbstractArray{NF,3},
                             tr_tend::AbstractArray{NF,3},
                             ps::AbstractArray{NF,4}) where {NF<:AbstractFloat}

#COMMENT: this is an ugly subroutine. Better to quantise this into smaller subroutines for readability and consistency

    #Declare empty mean arrays which we will fill
    umean = Array{NF,2}(undef,nlon, nlat) 
    vmean = similar(umean)
    dmean = similar(umean)
    dumc = Array{NF,3}(mx,nx,3) #COMMENT: mx, nx will have to be defined, again in a struct like nlon, nlat


    for k in 1:nlev
        umean += u_grid[:,:,k]*dhs[k] #COMMENT: dhs needs to be defined. 
        vmean += v_grid[:,:,k]*dhs[k]
        dmean += div_grid[:,:,k]*dhs[k]
    end

    # 1. Compute tendency of log(surface pressure)
    grad!(pₛ[:,:,j2], @view(dumc[:,:,1]), @view(dumc[:,:,2])) #COMMENT: grad! needs to be defined
    px = convert_to_grid(dumc[:,:,1], scale=true)
    py = convert_to_grid(dumc[:,:,2], scale=true)

    ps_tend = convert_to_spectral(-umean.*px - vmean.*py)
    ps_tend[1,1] = Complex{RealType}(zero) #COMMENT: how to do this?

    # 2. Compute "vertical" velocity
    sigma_tend = Array{NF,3}(undef,nlon, nlat,nlev+1) 
    sigma_m = Array{NF,3}(undef,nlon, nlat,nlev+1) 
    puv = Array{NF,3}(undef,nlon, nlat,nlev) 

    for k in 1:nlev
        puv[:,:,k] = (u_grid[:,:,k] - umean).*px + (v_grid[:,:,k] - vmean).*py
    end

    for k in 1:nlev
        sigma_tend[:,:,k+1] = sigma_tend[:,:,k] - dhs[k]*(puv[:,:,k] + D_grid[:,:,k] - dmean)
        sigma_m[:,:,k+1] = sigma_m[:,:,k] - dhs[k]*puv[:,:,k]
    end

    # 3. Subtract part of temperature field that is used as reference for implicit terms
    t_grid_anom = Array{NF,3}(undef,nlon, nlat,nlev) 
    for k in 1:nlev
        t_grid_anom[:,:,k] = t_grid[:,:,k] .- tref[k] #COMMENT: tref needs to be defined
    end

    # 4. Zonal wind tendency
    temp = Array{NF,3}(undef,nlon, nlat,nlev+1) 

    for k in 2:nlev
        temp[:,:,k] = sigma_tend[:,:,k].*(u_grid[:,:,k] - u_grid[:,:,k-1])
    end

    for k in 1:nlev
        u_tend[:,:,k] = u_tend[:,:,k] + v_grid[:,:,k].*vor_grid[:,:,k] - t_grid_anom[:,:,k]*R.*px
            - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
    end

    # 5. Meridional wind tendency
    for k in 2:nlev
        temp[:,:,k] = sigma_tend[:,:,k].*(v_grid[:,:,k] - v_grid[:,:,k-1])
    end

    for k in 1:nlev
        v_tend[:,:,k] = v_tend[:,:,k] -u_grid[:,:,k].*vor_grid[:,:,k] - t_grid_anom[:,:,k]*R.*py
            - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
    end

    # 6. Temperature tendency
    for k in 2:nlev
        temp[:,:,k] = sigma_tend[:,:,k].*(t_grid_anom[:,:,k] - t_grid_anom[:,:,k-1])
            + sigma_m[:,:,k].*(tref[k] - tref[k-1])
    end

    for k in 1:nlev
        t_tend[:,:,k] = t_tend[:,:,k]+ t_grid_anom[:,:,k].*div_grid[:,:,k]
            - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
            + fsgr[k]*t_grid_anom[:,:,k].*(σ_tend[:,:,k+1] + sigma_tend[:,:,k])
            + tref3[k]*(σ_m[:,:,k+1] + σ_m[:,:,k])
            + akap*(Tₐ_grid[:,:,k].*puv[:,:,k] - t_grid_anom[:,:,k].*dmean)
    end #COMMENT: lots of extras to be defined here

    # 7. Tracer tendency
    for itr in 1:n_trace
        for k in 2:nlev
            temp[:,:,k] = sigma_tend[:,:,k].*(tr_grid[:,:,k,itr] - tr_grid[:,:,k-1,itr])
        end

        temp[:,:,2:3] .= zero

        for k in 1:nlev
            tr_tend[:,:,k,itr] = tr_tend + tr_grid[:,:,k,itr].*D_grid[:,:,k]
                - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
        end
    end



    # =========================================================================
    # Convert tendencies to spectral space
    # =========================================================================

    for k in 1:nlev
        #  Convert u and v tendencies to vor and div spectral tendencies
        #  vdspec takes a grid u and a grid v and converts them to
        #  spectral vor and div
        vdspec!(u_tend[:,:,k], v_tend[:,:,k],
            @view(vor_tend[:,:,k]), @view(div_tend[:,:,k]), true)

        # Divergence tendency
        # add -lapl(0.5*(u**2+v**2)) to div tendency
        div_tend[:,:,k] = div_tend[:,:,k]
            - ∇²(convert_to_spectral(half*(u_grid[:,:,k].^two + v_grid[:,:,k].^two)))

        # Temperature tendency
        # and add div(vT) to spectral t tendency
        vdspec!(-u_grid[:,:,k].*t_grid_anom[:,:,k], -v_grid[:,:,k].*t_grid_anom[:,:,k],
            dumc[:,:,1], @view(t_tend[:,:,k]), true)
            t_tend[:,:,k] += convert_to_spectral(t_grid[:,:,k])

        # tracer tendency
        for itr in 1:n_trace
            vdspec!(-u_grid[:,:,k].*tr_grid[:,:,k,itr], -v_grid[:,:,k].*tr_grid[:,:,k,itr],
                dumc[:,:,1], @view(tr_tend[:,:,k,itr]), true)
            tr_tend[:,:,k,itr] += convert_to_spectral(tr_tend_grid[:,:,k,itr])
        end
    end
end




"""
Compute the spectral tendencies.  
"""
function get_spectral_tendencies!(D_tend, Tₐ_tend, pₛ_tend, j2)
    # Vertical mean div and pressure tendency
    dmeanc = zeros(Complex{RealType}, mx, nx)
    for k in 1:nlev
        dmeanc = dmeanc + D[:,:,k,j2]*dhs[k]
    end

    pₛ_tend = pₛ_tend - dmeanc
    pₛ_tend[1,1] = Complex{RealType}(zero)

    # Sigma-dot "velocity" and temperature tendency
    σ_tendc = zeros(Complex{RealType}, mx, nx, nlev+1)

    for k in 1:nlev - 1
        σ_tendc[:,:,k+1] = σ_tendc[:,:,k] - dhs[k]*(D[:,:,k,j2] - dmeanc)
    end

    dumk = zeros(Complex{RealType}, mx, nx, nlev+1)

    for k in 2:nlev
        dumk[:,:,k] = σ_tendc[:,:,k]*(tref[k] - tref[k-1])
    end

    for k in 1:nlev
        Tₐ_tend[:,:,k] -= (dumk[:,:,k+1] + dumk[:,:,k])*dhsr[k]
            + tref3[k]*(σ_tendc[:,:,k+1] + σ_tendc[:,:,k]) - tref2[k]*dmeanc
    end

    # Geopotential and divergence tendency
    get_geopotential(Tₐ[:,:,:,j2])
    for k in 1:nlev
        D_tend[:,:,k] -= ∇²(ϕ[:,:,k] + R*tref[k]*pₛ[:,:,j2])
    end
end