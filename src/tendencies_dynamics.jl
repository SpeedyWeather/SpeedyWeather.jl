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

    gradient_longitude!(∂uω_∂lon,uω)            # STEP 5: spectral gradients
    gradient_latitude!( ∂vω_∂lat,vω,S)

    flipsign!(∂uω_∂lon)                         # because ∂ζ/∂t = -∇⋅(uv*(ζ+f))       
    flipsign!(∂vω_∂lat)

    add_tendencies!(vor_tend,∂uω_∂lon,∂vω_∂lat) # STEP 6: Add tendencies
end

function curl_vorticity_fluxes!(D::DiagnosticVariables{NF}, # all diagnostic variables   
                                G::GeoSpectral{NF}          # struct with geometry and spectral transform
                                ) where {NF<:AbstractFloat}                                   

    S = G.spectral_transform
    @unpack uω, vω = D.intermediate_variables
    @unpack ∂uω_∂lat, ∂vω_∂lon = D.intermediate_variables
    @unpack div_tend = D.tendencies

    gradient_longitude!(∂vω_∂lon,vω)                    # 1st component of ∇×(uv(ζ+f))
    gradient_latitude!( ∂uω_∂lat,uω,S,flipsign=true)    # 2nd component of ∇×(uv(ζ+f))
    
    add_tendencies!(div_tend,∂vω_∂lon,∂uω_∂lat)         # evaluate after bernoulli_potential!
end

function vorticity_fluxes!( uω::AbstractMatrix{NF},     # Output: u*(vor+coriolis) in grid space
                            vω::AbstractMatrix{NF},     # Output: v*(vor+coriolis) in grid space
                            u::AbstractMatrix{NF},      # Input: zonal velocity in grid space
                            v::AbstractMatrix{NF},      # Input: meridional velocity in grid space
                            vor::AbstractMatrix{NF},    # Input: relative vorticity in grid space       
                            G::Geometry{NF}             # struct with precomputed geometry arrays
                            ) where {NF<:AbstractFloat} # number format NF

    nlon,nlat = size(u)
    @boundscheck size(u) == size(v) || throw(BoundsError)
    @boundscheck size(u) == size(vor) || throw(BoundsError)
    @boundscheck size(u) == size(uω) || throw(BoundsError)
    @boundscheck size(u) == size(vω) || throw(BoundsError)

    @unpack f_coriolis = G
    @boundscheck length(f_coriolis) == nlat || throw(BoundsError)

    @inbounds for j in 1:nlat
        for i in 1:nlon
            ω = vor[i,j] + f_coriolis[j]    # = relative vorticity + coriolis
            uω[i,j] = ω*u[i,j]              # = u(vor+f)
            vω[i,j] = ω*v[i,j]              # = v(vor+f)
        end
    end
end

function volume_fluxes!(D::DiagnosticVariables{NF}, # all diagnostic variables   
                        G::GeoSpectral{NF},         # struct with geometry and spectral transform
                        B::Boundaries,              # used for orography
                        H₀::Real                    # layer thickness
                        ) where {NF<:AbstractFloat}                                   

    @unpack pres_tend = D.tendencies
    @unpack u_grid, v_grid, pres_grid = D.grid_variables
    @unpack uh, vh, uh_grid, vh_grid = D.intermediate_variables
    @unpack ∂uh_∂lon, ∂vh_∂lat = D.intermediate_variables
    @unpack orography = B
    S = G.spectral_transform

    @unpack nlon, nlat = G.geometry
    @boundscheck size(pres_grid) == (nlon,nlat) || throw(BoundsError)
    @boundscheck size(pres_grid) == size(orography) || throw(BoundsError)
    @boundscheck size(uh_grid) == size(u_grid) || throw(BoundsError)
    @boundscheck size(vh_grid) == size(v_grid) || throw(BoundsError)

    H₀ = convert(NF,H₀)

    # compute (uh,vh) on the grid
    # pres_grid is η, the interface displacement
    # layer thickness h = η + H, H is the layer thickness at rest
    # H = H₀ - orography, H₀ is the layer thickness without mountains

    @inbounds for j in 1:nlat
        for i in 1:nlon
            h = pres_grid[i,j] + H₀ - orography[i,j]    # h = η + H₀ - orography
            uh_grid[i,j] = u_grid[i,j]*h                # = uh
            vh_grid[i,j] = v_grid[i,j]*h                # = vh
        end
    end

    spectral!(uh,uh_grid,S)
    spectral!(vh,vh_grid,S)

    gradient_longitude!(∂uh_∂lon,uh,  flipsign=true)    # 1st component of -∇⋅(uh)
    gradient_latitude!( ∂vh_∂lat,vh,S,flipsign=true)    # 2nd component of -∇⋅(uh)

    ∂uh_∂lon_surf = view(∂uh_∂lon,:,:,1)                # create views of surface layer
    ∂vh_∂lat_surf = view(∂vh_∂lat,:,:,1)

    add_tendencies!(pres_tend,∂uh_∂lon_surf,∂vh_∂lat_surf)
end

"""
    bernoulli_potential!(   B::AbstractMatrix{NF},  # Output: Bernoulli potential B = 1/2*(u^2+v^2)+Φ
                            u::AbstractMatrix{NF},  # zonal velocity
                            v::AbstractMatrix{NF},  # meridional velocity
                            η::AbstractMatrix{NF},  # interface displacement
                            g::Real                 # gravity
                            ) where {NF<:AbstractFloat}

Computes the Bernoulli potential 1/2*(u^2 + v^2) + g*η."""
function bernoulli_potential!(  B::AbstractMatrix{NF},  # Output: Bernoulli potential B = 1/2*(u^2+v^2)+Φ
                                u::AbstractMatrix{NF},  # zonal velocity
                                v::AbstractMatrix{NF},  # meridional velocity
                                η::AbstractMatrix{NF},  # interface displacement
                                g::Real                 # gravity
                                ) where {NF<:AbstractFloat}
    
    nlon,nlat = size(B)
    @boundscheck (nlon,nlat) == size(u) || throw(BoundsError)
    @boundscheck (nlon,nlat) == size(v) || throw(BoundsError)
    @boundscheck (nlon,nlat) == size(η) || throw(BoundsError)

    one_half = convert(NF,0.5)
    gravity = convert(NF,g)

    @inbounds for j in 1:nlat
        for i in 1:nlon
            B[i,j] = one_half*(u[i,j]^2 + v[i,j]^2) + gravity*η[i,j]
        end
    end
end

function bernoulli_potential!(  D::DiagnosticVariables{NF}, # all diagnostic variables   
                                G::GeoSpectral{NF},         # struct with geometry and spectral transform
                                g::Real                     # gravity
                                ) where {NF<:AbstractFloat}   
    
    @unpack u_grid,v_grid,pres_grid = D.grid_variables
    @unpack bernoulli, bernoulli_grid = D.intermediate_variables
    @unpack div_tend = D.tendencies
    S = G.spectral_transform

    bernoulli_potential!(bernoulli_grid,u_grid,v_grid,pres_grid,g)  # = 1/2(u^2 + v^2) + gη on grid
    spectral!(bernoulli,bernoulli_grid,S)                           # to spectral space

    # write directly in div_tend, ie bernoulli potential has to be the first tendency
    ∇²!(div_tend,bernoulli,S)                                       # = ∇²(1/2(u^2 + v^2) + gη)
    flipsign!(div_tend)                                              # = -∇²(1/2(u^2 + v^2) + gη) on RHS
end

"""
    gridded!(   diagn::DiagnosticVariables{NF}, # all diagnostic variables
                progn::PrognosticVariables{NF}, # all prognostic variables
                M::BarotropicModel,             # everything that's constant
                lf::Int=1                       # leapfrog index
                ) where NF

Propagate the spectral state of the prognostic variables `progn` to the
diagnostic variables in `diagn`. Updates grid vorticity, spectral stream function
and spectral and grid velocities u,v."""
function gridded!(  diagn::DiagnosticVariables{NF}, # all diagnostic variables
                    progn::PrognosticVariables{NF}, # all prognostic variables
                    M::BarotropicModel,             # everything that's constant
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
    gradient_latitude!( coslat_u, stream_function, S, flipsign=true)

    gridded!(u_grid,coslat_u,S)     # get u,v on grid from spectral
    gridded!(v_grid,coslat_v,S)

    unscale_coslat!(u_grid,G)       # undo the coslat scaling from gradients     
    unscale_coslat!(v_grid,G)

    return nothing
end

"""
    gridded!(   diagn::DiagnosticVariables{NF}, # all diagnostic variables
                progn::PrognosticVariables{NF}, # all prognostic variables
                M::ShallowWaterModel,           # everything that's constant
                lf::Int=1                       # leapfrog index
                ) where NF

Propagate the spectral state of the prognostic variables `progn` to the
diagnostic variables in `diagn`. Updates grid vorticity, spectral stream function
and spectral and grid velocities u,v."""
function gridded!(  diagn::DiagnosticVariables{NF}, # all diagnostic variables
                    progn::PrognosticVariables{NF}, # all prognostic variables
                    M::ShallowWaterModel,           # everything that's constant
                    lf::Int=1                       # leapfrog index
                    ) where NF
    
    @unpack vor, div, pres = progn                  # relative vorticity, divergence, pressure
    @unpack vor_grid, div_grid, u_grid, v_grid, pres_grid = diagn.grid_variables
    @unpack stream_function, coslat_u, coslat_v, velocity_potential = diagn.intermediate_variables
    
    G = M.geospectral.geometry
    S = M.geospectral.spectral_transform

    vor_lf = view(vor,:,:,lf,:)         # pick leapfrog index without memory allocation
    div_lf = view(div,:,:,lf,:)
    pres_lf = view(pres,:,:,lf)       

    gridded!(vor_grid,vor_lf,S)         # get vorticity on grid from spectral vor
    gridded!(div_grid,div_lf,S)         # get divergence on grid from spectral div
    gridded!(pres_grid,pres_lf,S)       # get pressure on grid from spectral pres

    ∇⁻²!(stream_function,vor_lf,S)      # invert Laplacian ∇² for stream function
    ∇⁻²!(velocity_potential,div_lf,S)   # invert Laplacian ∇² for velocity potential

    # contribution from stream function (non-divergent component)
    gradient_longitude!(coslat_v, stream_function)
    gradient_latitude!( coslat_u, stream_function, S, flipsign=true)

    # add contribution from velocity potential (non-rotational component)
    gradient_longitude!(coslat_u, velocity_potential,    add=true)
    gradient_latitude!( coslat_v, velocity_potential, S, add=true)

    gridded!(u_grid,coslat_u,S)     # get u,v on grid from spectral
    gridded!(v_grid,coslat_v,S)

    unscale_coslat!(u_grid,G)       # undo the coslat scaling from gradients     
    unscale_coslat!(v_grid,G)

    return nothing
end