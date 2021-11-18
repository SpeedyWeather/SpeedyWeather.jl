
"""
Compute the grid point and spectral tendencies 
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
                         j1::int,
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


    #1. Convert spectral prognostic to grid point space
    get_grid_point_fields(
                          
                          vor_grid,div_grid,t_grid,tr_grid, u_grid,v_grid,q_grid, ps_grid,phi_grid)


   # 2. Parameterised physics tendencies. NEEDS DEFINING in phypar.jl
   phypar(ξ_tend, D_tend, Tₐ_tend, tr_tend)


   ompute physical parametrization tendencies for u, v, t, q
    !
    ! Output arguments:
    !     utend: u-wind tendency (gp)
    !     vtend: v-wind tendency (gp)
    !     ttend: temp. tendency (gp)
    !     qtend: spec. hum. tendency (gp)


   #3. Dynamics tendencies
   dynamics_tendencies(,
   )

end


"""
Compute grid point fields of the spectral prognostic tendencies
"""
function get_grid_point_fields(
                               vor,div,t,phi,ps,tr
                               vor_grid,div_grid,t_grid,tr_grid, u_grid,v_grid, q_grid, ps_grid,phi_grid)



    #1. Compute grid-point fields
    #1.1 Update geopotential in spectral space
    #COMMENT: THIS NEEDS TO BE CONFIGURED PROPERLY W.R.T ARGUMENTS.
    geopotential!(phi)

    #1.2 Grid point variables for physics tendencies
    #Calculate u and v in grid-point space from vorticity and
    #divergence in spectral space
    for j in 1:nlev
        uvspec!(vor[:,:,k,j2], div[:,:,k,j2], u_grid[:,:,k],v_grid[:,:,k]
    end
    
    #Truncate variables where the spectral transform is still done at double
    #precision
    u_grid = u_grid / 3600.0
    v_grid = v_grid / 3600.0



    for k in 1:nlev
        gridded(t[:,:,k,j2], t_grid[:,:,k])   #Transform from spectral to gridpoint space
    end

    #Normalise geopotential by cp to avoid overflows in physics
    for k in 1:nlev
        gridded(phi[:,:,k]*(1/cp), phi_grid[:,:,k])   #Transform from spectral to gridpoint space
    end


    #Don't transform the two stratospheric levels where humidity is set to zero
    # because it leads to overflows 
    for k in 3:nlev
        gridded(tr(:,:,k,j1,1), q_grid(:.:,k), )
    end 

    #Surface pressure spectral transform to grid
    gridded(ps(:,:,j1), ps_grid)


    #1.3 Grid-point variables for dynamics tendencies
    #Set units of vorticity and divergence to 'per hour' to reduce underflows 

   for k in 1:nlev
    vor_grid[:,:,k] = gridded(vor[:,:,k,j2])
    div_grid[:,:,k] = gridded(div[:,:,k,j2])

    for itr in 1:n_trace
        tr_grid[:,:,k,itr] = gridded(tr[:,:,k,j2,itr])
    end

    for j in 1:nlat
        for i in 1:nlon
            vor_grid[i,j,k] += coriol[j] #add coriolo
        end
    end

   end



end

"""
Compute non-linear tendencies in grid-point space from ynamics and add to physics tendencies. Convert total
gridpoint tendencies to spectral tendencies
"""
function dynamics_tendencies()


    umean = zeros(RealType, nlon, nlat)
    vmean = zeros(RealType, nlon, nlat)
    dmean = zeros(RealType, nlon, nlat)

    for k in 1:nlev
        umean += u_grid[:,:,k]*dhs[k]
        vmean += v_grid[:,:,k]*dhs[k]
        dmean += D_grid[:,:,k]*dhs[k]
    end

    # Compute tendency of log(surface pressure)
    # pₛ(1,1,j2)=zero
    grad!(pₛ[:,:,j2], @view(dumc[:,:,1]), @view(dumc[:,:,2]))
    px = spec_to_grid(dumc[:,:,1], scale=true)
    py = spec_to_grid(dumc[:,:,2], scale=true)

    pₛ_tend = grid_to_spec(-umean.*px - vmean.*py)
    pₛ_tend[1,1] = Complex{RealType}(zero)

    # Compute "vertical" velocity
    σ_tend = zeros(RealType, nlon, nlat, nlev+1)
    σ_m = zeros(RealType, nlon, nlat, nlev+1)

    # (The following combination of terms is utilized later in the
    #     temperature equation)
    puv = zeros(RealType, nlon, nlat, nlev)
    for k in 1:nlev
        puv[:,:,k] = (u_grid[:,:,k] - umean).*px + (v_grid[:,:,k] - vmean).*py
    end

    for k in 1:nlev
        σ_tend[:,:,k+1] = σ_tend[:,:,k] - dhs[k]*(puv[:,:,k] + D_grid[:,:,k] - dmean)
        σ_m[:,:,k+1] = σ_m[:,:,k] - dhs[k]*puv[:,:,k]
    end

    # Subtract part of temperature field that is used as reference for implicit terms
    t_grid_anom = zeros(RealType, nlon, nlat, nlev)

    for k in 1:nlev
        t_grid_anom[:,:,k] = Tₐ_grid[:,:,k] .- tref[k]
    end

    # Zonal wind tendency
    temp = zeros(RealType, nlon, nlat, nlev+1)
    u_tend = zeros(RealType, nlon, nlat, nlev)

    for k in 2:nlev
        temp[:,:,k] = σ_tend[:,:,k].*(u_grid[:,:,k] - u_grid[:,:,k-1])
    end

    for k in 1:nlev
        u_tend[:,:,k] = v_grid[:,:,k].*ξ_grid[:,:,k] - t_grid_anom[:,:,k]*R.*px
            - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
    end

    # Meridional wind tendency
    v_tend = zeros(RealType, nlon, nlat, nlev)

    for k in 2:nlev
        temp[:,:,k] = σ_tend[:,:,k].*(v_grid[:,:,k] - v_grid[:,:,k-1])
    end

    for k in 1:nlev
        v_tend[:,:,k] = -u_grid[:,:,k].*ξ_grid[:,:,k] - t_grid_anom[:,:,k]*R.*py
            - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
    end

    # Temperature tendency
    Tₐ_tend_grid = zeros(RealType, nlon, nlat, nlev)

    for k in 2:nlev
        temp[:,:,k] = σ_tend[:,:,k].*(t_grid_anom[:,:,k] - t_grid_anom[:,:,k-1])
            + σ_m[:,:,k].*(tref[k] - tref[k-1])
    end

    for k in 1:nlev
        Tₐ_tend_grid[:,:,k] = t_grid_anom[:,:,k].*D_grid[:,:,k]
            - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
            + fsgr[k]*t_grid_anom[:,:,k].*(σ_tend[:,:,k+1] + σ_tend[:,:,k])
            + tref3[k]*(σ_m[:,:,k+1] + σ_m[:,:,k])
            + akap*(Tₐ_grid[:,:,k].*puv[:,:,k] - t_grid_anom[:,:,k].*dmean)
    end

    # Tracer tendency
    tr_tend_grid = zeros(RealType, nlon, nlat, nlev, n_trace)
    for itr in 1:n_trace
        for k in 2:nlev
            temp[:,:,k] = σ_tend[:,:,k].*(tr_grid[:,:,k,itr] - tr_grid[:,:,k-1,itr])
        end

        temp[:,:,2:3] .= zero

        for k in 1:nlev
            tr_tend_grid[:,:,k,itr] = tr_grid[:,:,k,itr].*D_grid[:,:,k]
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
            @view(ξ_tend[:,:,k]), @view(D_tend[:,:,k]), true)

        # Divergence tendency
        # add -lapl(0.5*(u**2+v**2)) to div tendency
        D_tend[:,:,k] = D_tend[:,:,k]
            - ∇²(grid_to_spec(half*(u_grid[:,:,k].^two + v_grid[:,:,k].^two)))

        # Temperature tendency
        # and add div(vT) to spectral t tendency
        vdspec!(-u_grid[:,:,k].*t_grid_anom[:,:,k], -v_grid[:,:,k].*t_grid_anom[:,:,k],
            dumc[:,:,1], @view(Tₐ_tend[:,:,k]), true)
        Tₐ_tend[:,:,k] += grid_to_spec(Tₐ_tend_grid[:,:,k])

        # tracer tendency
        for itr in 1:n_trace
            vdspec!(-u_grid[:,:,k].*tr_grid[:,:,k,itr], -v_grid[:,:,k].*tr_grid[:,:,k,itr],
                dumc[:,:,1], @view(tr_tend[:,:,k,itr]), true)
            tr_tend[:,:,k,itr] += grid_to_spec(tr_tend_grid[:,:,k,itr])
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