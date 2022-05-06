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
    @unpack u_grid,v_grid,div_grid               = Diag.grid_variables
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
    pres_surf_gradient_grid_x = gridded(pres_surf_gradient_spectral_x*3600)
    pres_surf_gradient_grid_y = gridded(pres_surf_gradient_spectral_x*3600) #3600 factor from Paxton/Chantry. I think this is to correct for the underflow rescaling earlier

    pres_surf_tend = spectral(-u_mean.*pres_surf_gradient_grid_x - v_mean.*pres_surf_gradient_grid_y)
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
    @unpack tref = M.GeoSpectral.geometry #Note that tref is currently not defined correctly 

    _,_,nlev = size(temp_grid)


    for k in 1:nlev
        temp_grid_anomaly[:,:,k] = temp_grid[:,:,k] .- tref[k] #+ 0K correction?
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
function temperature_tendency!(Diag::DiagnosticVariables{NF}, # Diagnostic variables
                               M
                               )where {NF<:AbstractFloat}
    

    @unpack temp_tend = Diag.tendencies
    @unpack div_grid,temp_grid,temp_grid_anomaly = Diag.grid_variables
    @unpack sigma_u,sigma_tend,sigma_m,puv,div_mean = Diag.intermediate_variables
    @unpack tref,σ_levels_half⁻¹_2,fsgr,tref3 = M.GeoSpectral.geometry #Note that tref is currenrtly not defined correctly 
    @unpack akap = M.Parameters

    _,_,nlev = size(div_grid)



    for k in 2:nlev
        sigma_u[:,:,k] = sigma_tend[:,:,k].*(temp_grid_anomaly[:,:,k] - temp_grid_anomaly[:,:,k-1])
                    + sigma_m[:,:,k].*(tref[k] - tref[k-1])
    end

    for k in 1:nlev
        temp_tend[:,:,k] = temp_tend[:,:,k]
                        + temp_grid_anomaly[:,:,k].*div_grid[:,:,k]
                        - (sigma_u[:,:,k+1] + sigma_u[:,:,k])*σ_levels_half⁻¹_2[k]
                        + fsgr[k]*temp_grid_anomaly[:,:,k].*(sigma_tend[:,:,k+1] + sigma_tend[:,:,k])
                        + tref3[k]*(sigma_m[:,:,k+1] + sigma_m[:,:,k])
                        + akap*(temp_grid[:,:,k].*puv[:,:,k] - temp_grid_anomaly[:,:,k].*div_mean)
    end 


end

"""
Compute the humidity tendency
"""
function humidity_tendency!(Diag::DiagnosticVariables{NF}, # Diagnostic variables
                            M
                            )where {NF<:AbstractFloat}
    
    @unpack humid_tend = Diag.tendencies
    @unpack div_grid,humid_grid = Diag.grid_variables
    @unpack sigma_u,sigma_tend,= Diag.intermediate_variables

    @unpack σ_levels_half⁻¹_2 = M.GeoSpectral.geometry 

    _,_,nlev = size(div_grid)


    for k in 2:nlev
        sigma_u[:,:,k] = sigma_tend[:,:,k].*(humid_grid[:,:,k] - humid_grid[:,:,k-1])
    end
    

    # From Paxton/Chantry: dyngrtend.f90. Unsure if we need this here since we are dealing solely with humidity,
        # !spj for moisture, vertical advection is not possible between top
        # !spj two layers
        # !kuch three layers
        # !if(iinewtrace==1)then
        # do k=2,3
        #     temp(:,:,k)=0.0_dp # temp is equivalent to sigma_u. i.e a temporary array that is reused for calculations 
        # enddo
        # !endif


    for k in 1:nlev
        humid_tend[:,:,k] = humid_tend[:,:,k]
                            + humid_grid[:,:,k].*div_grid[:,:,k]
                            - (sigma_u[:,:,k+1] + sigma_u[:,:,k])*σ_levels_half⁻¹_2[k]
        end 

end

"""Spectral tendency of ∇⋅(uv*ω) from vector uv=(u,v) in grid space and absolute vorticity ω.
Step 1 (grid space): Add Coriolis f to the relative vorticity ζ (=`vor_grid`) to obtain abs vorticity ω.
Step 2 (grid space): Multiply u,v with abs vorticity ω.
Step 3 (grid space): Unscale with coslat, cosine of latitude, as the gradients will include a coslat term.
Step 4 (spectral space): convert uω/coslat, vω/coslat from grid to spectral space
Step 5 (spectral space): Compute gradients ∂/∂lon(uω/coslat) and ∂/∂lat(vω/coslat)
Step 6 (spectral space): Add ∂/∂lon(uω/coslat)+∂/∂θ(vω/coslat) and return.
"""
function vorticity_advection!(  D::DiagnosticVariables{NF}, # all diagnostic variables   
                                G::GeoSpectral{NF}          # struct with geometry and spectral transform
                                ) where {NF<:AbstractFloat}                                   

    S = G.spectral_transform
    @unpack u_grid, v_grid, vor_grid = D.grid_variables
    @unpack uω, vω, uω_grid,vω_grid = D.intermediate_variables
    @unpack ∂uω_∂lon,∂vω_∂lat = D.intermediate_variables
    @unpack vor_tend = D.tendencies

    # STEP 1-3: Abs vorticity, velocity times abs vort
    vorticity_fluxes!(uω_grid,vω_grid,u_grid,v_grid,vor_grid,G.geometry)

    spectral!(uω,uω_grid,S)                     # STEP 4: to spectral space
    spectral!(vω,vω_grid,S)

    gradient_longitude!(∂uω_∂lon,-uω)           # STEP 5: spectral gradients
    gradient_latitude!( ∂vω_∂lat,vω,S)

    add_tendencies!(vor_tend,∂uω_∂lon,∂vω_∂lat) # STEP 6: Add tendencies
end

function vorticity_fluxes!( uω::AbstractMatrix{NF},    # Output: u*(vor+coriolis) in grid space
                            vω::AbstractMatrix{NF},    # Output: v*(vor+coriolis) in grid space
                            u::AbstractMatrix{NF},     # Input: zonal velocity in grid space
                            v::AbstractMatrix{NF},     # Input: meridional velocity in grid space
                            vor::AbstractMatrix{NF},   # Input: relative vorticity in grid space       
                            G::Geometry{NF}            # struct with precomputed geometry arrays
                            ) where {NF<:AbstractFloat}# number format NF

    nlon,nlat = size(u)
    @boundscheck size(u) == size(v) || throw(BoundsError)
    @boundscheck size(u) == size(vor) || throw(BoundsError)
    @boundscheck size(u) == size(uω) || throw(BoundsError)
    @boundscheck size(u) == size(vω) || throw(BoundsError)

    @unpack f_coriolis = G

    @inbounds for j in 1:nlat
        for i in 1:nlon
            ω = vor[i,j] + f_coriolis[j]    # = relative vorticity + coriolis
            uω[i,j] = ω*u[i,j]              # = u(vor+f)
            vω[i,j] = ω*v[i,j]              # = v(vor+f)
        end
    end
end

function gridded!(  diagn::DiagnosticVariables{NF}, # all diagnostic variables
                    progn::PrognosticVariables{NF}, # all prognostic variables
                    M::ModelSetup,                  # everything that's constant
                    lf::Int=1                       # leapfrog index
                    ) where NF
    
    @unpack vor = progn                             # relative vorticity
    @unpack vor_grid, u_grid, v_grid = diagn.grid_variables
    @unpack stream_function, coslat_u, coslat_v = diagn.intermediate_variables
    
    G = M.geospectral.geometry
    S = M.geospectral.spectral_transform

    vor_lf = view(vor,:,:,lf,:)     # pick leapfrog index without memory allocation
    gridded!(vor_grid,vor_lf,S)     # get vorticity on grid from spectral vor
    ∇⁻²!(stream_function,vor_lf,S)  # invert Laplacian ∇² for stream function
    
    gradient_longitude!(coslat_v, stream_function)
    gradient_latitude!( coslat_u, stream_function, S)

    gridded!(u_grid,coslat_u,S)     # get u,v on grid from spectral
    gridded!(v_grid,coslat_v,S)

    unscale_coslat!(u_grid,G)       # undo the coslat scaling from gradients     
    unscale_coslat!(v_grid,G)

    return nothing
end