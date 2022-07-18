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

"""
    vorticity_flux_divergence!( D::DiagnosticVariables{NF}, # all diagnostic variables   
                                G::GeoSpectral{NF}          # struct with geometry and spectral transform
                                ) where {NF<:AbstractFloat}

Compute the vorticity advection as the (negative) divergence of the vorticity fluxes -∇⋅(uv*(ζ+f)).
First, compute the uv*(ζ+f), then transform to spectral space and take the divergence and flip the sign."""
function vorticity_flux_divergence!(D::DiagnosticVariables{NF}, # all diagnostic variables   
                                    G::GeoSpectral{NF}          # struct with geometry and spectral transform
                                    ) where {NF<:AbstractFloat}                                   

    S = G.spectral_transform
    @unpack U_grid, V_grid, vor_grid = D.grid_variables
    @unpack uω_coslat⁻¹, vω_coslat⁻¹, uω_coslat⁻¹_grid, vω_coslat⁻¹_grid = D.intermediate_variables
    @unpack vor_tend = D.tendencies

    # STEP 1-3: Abs vorticity, velocity times abs vort
    vorticity_fluxes!(uω_coslat⁻¹_grid,vω_coslat⁻¹_grid,U_grid,V_grid,vor_grid,G.geometry)

    spectral!(uω_coslat⁻¹,uω_coslat⁻¹_grid,S)
    spectral!(vω_coslat⁻¹,vω_coslat⁻¹_grid,S)

    # flipsign as RHS is negative ∂ζ/∂t = -∇⋅(uv*(ζ+f)), write directly into tendency
    divergence!(vor_tend,uω_coslat⁻¹,vω_coslat⁻¹,S,flipsign=true)
end

"""
    vorticity_flux_curl!(   D::DiagnosticVariables{NF}, # all diagnostic variables   
                            G::GeoSpectral{NF}          # struct with geometry and spectral transform
                            ) where {NF<:AbstractFloat}
    
Compute the curl of the vorticity fluxes ∇×(uω,vω) and add to divergence tendency."""
function vorticity_flux_curl!(  D::DiagnosticVariables{NF}, # all diagnostic variables   
                                G::GeoSpectral{NF}          # struct with geometry and spectral transform
                                ) where {NF<:AbstractFloat}                                   

    S = G.spectral_transform
    @unpack uω_coslat⁻¹, vω_coslat⁻¹ = D.intermediate_variables
    @unpack div_tend = D.tendencies

    # write directly into divergence tendency, add as Bernoulli potential is already in div_tend
    curl!(div_tend,uω_coslat⁻¹,vω_coslat⁻¹,S,add=true)
end

"""
    vorticity_fluxes!(  uω_coslat⁻¹::AbstractMatrix{NF},    # Output: u*(ζ+f)/coslat in grid space
                        vω_coslat⁻¹::AbstractMatrix{NF},    # Output: v*(ζ+f)/coslat in grid space
                        U::AbstractMatrix{NF},              # Input: u*coslat in grid space
                        V::AbstractMatrix{NF},              # Input: v*coslat in grid space
                        vor::AbstractMatrix{NF},            # Input: relative vorticity ζ in grid space       
                        G::Geometry{NF}                     # struct with precomputed geometry arrays
                        ) where {NF<:AbstractFloat}         # number format NF

Compute the vorticity fluxes (u,v)*(ζ+f)/coslat in grid-point space from U,V and vorticity ζ."""
function vorticity_fluxes!( uω_coslat⁻¹::AbstractMatrix{NF},    # Output: u*(ζ+f)/coslat in grid space
                            vω_coslat⁻¹::AbstractMatrix{NF},    # Output: v*(ζ+f)/coslat in grid space
                            U::AbstractMatrix{NF},              # Input: u*coslat in grid space
                            V::AbstractMatrix{NF},              # Input: v*coslat in grid space
                            vor::AbstractMatrix{NF},            # Input: relative vorticity ζ in grid space       
                            G::Geometry{NF}                     # struct with precomputed geometry arrays
                            ) where {NF<:AbstractFloat}         # number format NF

    nlon,nlat = size(U)
    @boundscheck size(U) == size(V) || throw(BoundsError)
    @boundscheck size(U) == size(vor) || throw(BoundsError)
    @boundscheck size(U) == size(uω_coslat⁻¹) || throw(BoundsError)
    @boundscheck size(U) == size(vω_coslat⁻¹) || throw(BoundsError)

    @unpack f_coriolis, coslat⁻² = G
    @boundscheck length(f_coriolis) == nlat || throw(BoundsError)
    @boundscheck length(coslat⁻²) == nlat || throw(BoundsError)

    @inbounds for j in 1:nlat
        for i in 1:nlon
            # ω = relative vorticity + coriolis and unscale with coslat²
            ω = coslat⁻²[j]*(vor[i,j] + f_coriolis[j])
            uω_coslat⁻¹[i,j] = ω*U[i,j]              # = u(ζ+f)/coslat
            vω_coslat⁻¹[i,j] = ω*V[i,j]              # = v(ζ+f)/coslat
        end
    end
end

"""
    volume_fluxes!( D::DiagnosticVariables{NF}, # all diagnostic variables   
                    G::GeoSpectral{NF},         # struct with geometry and spectral transform
                    B::Boundaries,              # contains orography
                    H₀::Real                    # layer thickness
                    ) where {NF<:AbstractFloat}   

Computes the (negative) divergence of the volume fluxes `uh,vh` for the continuity equation, -∇⋅(uh,vh)"""
function volume_flux_divergence!(   D::DiagnosticVariables{NF}, # all diagnostic variables   
                                    G::GeoSpectral{NF},         # struct with geometry and spectral transform
                                    B::Boundaries,              # contains orography
                                    H₀::Real                    # layer thickness
                                    ) where {NF<:AbstractFloat}                                   

    @unpack pres_tend = D.tendencies
    @unpack U_grid, V_grid, pres_grid = D.grid_variables
    @unpack uh_coslat⁻¹, vh_coslat⁻¹, uh_coslat⁻¹_grid, vh_coslat⁻¹_grid = D.intermediate_variables
    @unpack orography = B
    S = G.spectral_transform
    @unpack coslat⁻² = G.geometry

    @unpack nlon, nlat = G.geometry
    @boundscheck size(pres_grid) == (nlon,nlat) || throw(BoundsError)
    @boundscheck size(pres_grid) == size(orography) || throw(BoundsError)
    @boundscheck size(uh_coslat⁻¹_grid) == size(U_grid) || throw(BoundsError)
    @boundscheck size(vh_coslat⁻¹_grid) == size(V_grid) || throw(BoundsError)
    @boundscheck length(coslat⁻²) == nlat || throw(BoundsError) 

    H₀ = convert(NF,H₀)

    # compute (uh,vh) on the grid
    # pres_grid is η, the interface displacement
    # layer thickness h = η + H, H is the layer thickness at rest
    # H = H₀ - orography, H₀ is the layer thickness without mountains

    k = 1
    @inbounds for j in 1:nlat
        for i in 1:nlon
            # h = (η + H₀ - orography)/coslat²
            h = coslat⁻²[j]*(pres_grid[i,j,k] + H₀ - orography[i,j])
            uh_coslat⁻¹_grid[i,j,k] = U_grid[i,j,k]*h      # = uh/coslat
            vh_coslat⁻¹_grid[i,j,k] = V_grid[i,j,k]*h      # = vh/coslat
        end
    end

    spectral!(uh_coslat⁻¹,uh_coslat⁻¹_grid,S)       # to spectral space
    spectral!(vh_coslat⁻¹,vh_coslat⁻¹_grid,S)

    uh_coslat⁻¹_surf = view(uh_coslat⁻¹,:,:,1)      # create views of surface layer
    vh_coslat⁻¹_surf = view(vh_coslat⁻¹,:,:,1)

    # compute divergence of volume fluxes and flip sign as ∂η/∂ = -∇⋅(uh,vh)
    divergence!(pres_tend,uh_coslat⁻¹_surf,vh_coslat⁻¹_surf,S,flipsign=true)
end

"""
    bernoulli_potential!(   B::AbstractMatrix{NF},  # Output: Bernoulli potential B = 1/2*(u^2+v^2)+g*η
                            u::AbstractMatrix{NF},  # zonal velocity
                            v::AbstractMatrix{NF},  # meridional velocity
                            η::AbstractMatrix{NF},  # interface displacement
                            g::Real                 # gravity
                            ) where {NF<:AbstractFloat}

Computes the Bernoulli potential 1/2*(u^2 + v^2) + g*η in grid-point space."""
function bernoulli_potential!(  B::AbstractMatrix{NF},  # Output: Bernoulli potential B = 1/2*(u^2+v^2)+Φ
                                U::AbstractMatrix{NF},  # zonal velocity *coslat
                                V::AbstractMatrix{NF},  # meridional velocity *coslat
                                η::AbstractMatrix{NF},  # interface displacement
                                g::Real,                # gravity
                                G::Geometry{NF}         # used for precomputed cos²(lat)
                                ) where {NF<:AbstractFloat}
    
    @unpack coslat⁻² = G
    nlon, nlat = size(B)
    @boundscheck nlat == length(coslat⁻²) || throw(BoundsError)
    @boundscheck size(B) == size(U) || throw(BoundsError)
    @boundscheck size(B) == size(V) || throw(BoundsError)
    @boundscheck size(B) == size(η) || throw(BoundsError)

    one_half = convert(NF,0.5)                          # convert to number format NF
    gravity = convert(NF,g)

    @inbounds for j in 1:nlat
        one_half_coslat⁻² = one_half*coslat⁻²[j]
        for i in 1:nlon
            B[i,j] = one_half_coslat⁻²*(U[i,j]^2 + V[i,j]^2) + gravity*η[i,j]
        end
    end
end

"""
    bernoulli_potential!(   D::DiagnosticVariables{NF}, # all diagnostic variables   
                            GS::GeoSpectral{NF},        # struct with geometry and spectral transform
                            g::Real                     # gravity
                            ) where {NF<:AbstractFloat}

Computes the Laplace operator ∇² of the Bernoulli potential `B` in spectral space. First, computes the Bernoulli potential
on the grid, then transforms to spectral space and takes the Laplace operator."""
function bernoulli_potential!(  D::DiagnosticVariables{NF}, # all diagnostic variables   
                                GS::GeoSpectral{NF},        # struct with geometry and spectral transform
                                g::Real                     # gravity
                                ) where {NF<:AbstractFloat}   
    
    @unpack U_grid,V_grid,pres_grid = D.grid_variables
    @unpack bernoulli, bernoulli_grid = D.intermediate_variables
    @unpack div_tend = D.tendencies
    S = GS.spectral_transform
    G = GS.geometry

    bernoulli_potential!(bernoulli_grid,U_grid,V_grid,pres_grid,g,G)# = 1/2(u^2 + v^2) + gη on grid
    spectral!(bernoulli,bernoulli_grid,S)                           # to spectral space

    # write directly in div_tend, ie bernoulli potential has to be the first tendency
    ∇²!(div_tend,bernoulli,S)                                       # = ∇²(1/2(u^2 + v^2) + gη)
    flipsign!(div_tend)                                             # = -∇²(1/2(u^2 + v^2) + gη) on RHS
end

"""
    gridded!(   diagn::DiagnosticVariables{NF}, # all diagnostic variables
                progn::PrognosticVariables{NF}, # all prognostic variables
                M::BarotropicModel,             # everything that's constant
                lf::Int=1                       # leapfrog index
                ) where NF

Propagate the spectral state of the prognostic variables `progn` to the
diagnostic variables in `diagn` for the barotropic vorticity model.
Updates grid vorticity, spectral stream function and spectral and grid velocities u,v."""
function gridded!(  diagn::DiagnosticVariables{NF}, # all diagnostic variables
                    progn::PrognosticVariables{NF}, # all prognostic variables
                    M::BarotropicModel,             # everything that's constant
                    lf::Int=1                       # leapfrog index
                    ) where NF
    
    @unpack vor = progn                             # relative vorticity
    @unpack vor_grid, U_grid, V_grid = diagn.grid_variables
    @unpack u_coslat, v_coslat = diagn.intermediate_variables
    S = M.geospectral.spectral_transform

    vor_lf = view(vor,:,:,lf,:)     # pick leapfrog index without memory allocation
    gridded!(vor_grid,vor_lf,S)     # get vorticity on grid from spectral vor
    
    # get spectral U,V from spectral vorticity via stream function Ψ
    # U = u*coslat = -coslat*∂Ψ/∂lat
    # V = v*coslat = ∂Ψ/∂lon, radius omitted in both cases
    UV_from_vor!(u_coslat,v_coslat,vor_lf,S)

    # transform to U,V on grid (U,V = u,v*coslat)
    gridded!(U_grid,u_coslat,S)
    gridded!(V_grid,v_coslat,S)

    return nothing
end

"""
    gridded!(   diagn::DiagnosticVariables{NF}, # all diagnostic variables
                progn::PrognosticVariables{NF}, # all prognostic variables
                M::ShallowWaterModel,           # everything that's constant
                lf::Int=1                       # leapfrog index
                ) where NF

Propagate the spectral state of the prognostic variables `progn` to the
diagnostic variables in `diagn` for the shallow water model. Updates grid vorticity,
grid divergence, grid interface displacement (`pres_grid`) and the velocities
U,V (scaled by cos(lat))."""
function gridded!(  diagn::DiagnosticVariables{NF}, # all diagnostic variables
                    progn::PrognosticVariables{NF}, # all prognostic variables
                    M::ShallowWaterModel,           # everything that's constant
                    lf::Int=1                       # leapfrog index
                    ) where NF
    
    @unpack vor, div, pres = progn                  # relative vorticity, divergence, pressure
    @unpack vor_grid, div_grid, U_grid, V_grid, pres_grid = diagn.grid_variables
    @unpack u_coslat, v_coslat = diagn.intermediate_variables
    S = M.geospectral.spectral_transform

    vor_lf = view(vor,:,:,lf,:)         # pick leapfrog index without memory allocation
    div_lf = view(div,:,:,lf,:)
    pres_lf = view(pres,:,:,lf)       

    gridded!(vor_grid,vor_lf,S)         # get vorticity on grid from spectral vor
    gridded!(div_grid,div_lf,S)         # get divergence on grid from spectral div
    gridded!(pres_grid,pres_lf,S)       # get pressure on grid from spectral pres

    # get spectral U,V from vorticity and divergence via stream function Ψ and vel potential ϕ
    # U = u*coslat = -coslat*∂Ψ/∂lat + ∂ϕ/dlon
    # V = v*coslat =  coslat*∂ϕ/∂lat + ∂Ψ/dlon
    UV_from_vordiv!(u_coslat,v_coslat,vor_lf,div_lf,S)

    # transform to U,V on grid (U,V = u,v*coslat)
    gridded!(U_grid,u_coslat,S)
    gridded!(V_grid,v_coslat,S)

    return nothing
end