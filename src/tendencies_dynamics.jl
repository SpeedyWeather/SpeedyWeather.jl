"""
    surface_pressure_tendency!( Prog::PrognosticVariables,
                                Diag::DiagnosticVariables,
                                lf::Int,
                                M::PrimitiveEquationModel)

Computes the tendency of the logarithm of surface pressure as

    -(ū*px + v̄*py) - D̄

with ū,v̄ being the vertically averaged velocities; px, py the gradients
of the logarithm of surface pressure ln(p_s) and D̄ the vertically averaged divergence.
1. Calculate ∇ln(p_s) in spectral space, convert to grid.
2. Multiply ū,v̄ with ∇ln(p_s) in grid-point space, convert to spectral.
3. D̄ is subtracted in spectral space.
4. Set tendency of the l=m=0 mode to 0 for better mass conservation."""
function surface_pressure_tendency!(progn::PrognosticVariables{NF},
                                    diagn::DiagnosticVariables{NF},
                                    lf::Int,                      # leapfrog index
                                    model::PrimitiveEquationModel
                                    ) where {NF<:AbstractFloat}

    # CALCULATE ∇lnp_s
    pres = progn.pres.leapfrog[lf]
    @unpack dpres_dlon, dpres_dlat = diagn.surface
    @unpack dpres_dlon_grid, dpres_dlat_grid = diagn.surface
    
    ∇!(dpres_dlon,dpres_dlat,pres,model.spectral_transform)
    gridded!(dpres_dlon_grid,dpres_dlon,model.spectral_transform)
    gridded!(dpres_dlat_grid,dpres_dlat,model.spectral_transform)

    # TENDENCY: -(ū,v̄)⋅∇lnp_s
    # vertical averages need to be computed first!
    @unpack pres_tend, pres_tend_grid = diagn.surface
    @unpack U_mean, V_mean, div_mean = diagn.surface
    @unpack coslat⁻² = model.geometry

    # precompute ring indices
    rings = eachring(pres_tend_grid,dpres_dlon_grid,dpres_dlat_grid,U_mean,V_mean)

    @inbounds for (j,ring) in enumerate(rings)
        coslat⁻²j = coslat⁻²[j]
        for ij in ring
            # -(ū,v̄)⋅∇lnp_s only, do -D̄ in spectral space
            pres_tend_grid[ij] = -(U_mean[ij]*dpres_dlon_grid[ij] +
                                    V_mean[ij]*dpres_dlat_grid[ij])*coslat⁻²j
        end
    end

    spectral!(pres_tend,pres_tend_grid,model.spectral_transform)

    # now do the -D̄ term
    @inbounds for lm in eachharmonic(pres_tend,div_mean)
        pres_tend[lm] -= div_mean[lm]
    end

    pres_tend[1] = zero(NF)     # for mass conservation
    return nothing
end

"""
    vertical_averages!(Diag::DiagnosticVariables,G::Geometry)

Calculates the vertically averaged (weighted by the thickness of the σ level)
velocities (*coslat) and divergence. E.g.

    U_mean = ∑_k=1^nlev Δσ_k * U_k

U,V are averaged in grid-point space, divergence in spectral space.
"""
function vertical_averages!(progn::PrognosticVariables{NF},
                            diagn::DiagnosticVariables{NF},
                            lf::Int,            # leapfrog index
                            G::Geometry{NF}) where NF
    
    @unpack σ_levels_thick, nlev = G
    @unpack U_mean, V_mean, div_mean = diagn.surface

    @boundscheck nlev == diagn.nlev || throw(BoundsError)

    fill!(U_mean,0)     # reset accumulators from previous vertical average
    fill!(V_mean,0)
    fill!(div_mean,0)

    for k in 1:nlev
        Δσ_k = σ_levels_thick[k]
        U = diagn.layers[k].grid_variables.U_grid
        V = diagn.layers[k].grid_variables.V_grid
        div = progn.layers[k].leapfrog[lf].div

        # UV in grid-point space
        @inbounds for ij in eachgridpoint(diagn.surface)
            U_mean[ij] += U[ij]*Δσ_k
            V_mean[ij] += V[ij]*Δσ_k
        end

        # divergence in spectral space
        @inbounds for lm in eachharmonic(div,div_mean)
            div_mean[lm] += div[lm]*Δσ_k
        end
    end
end
        
"""
Compute the spectral tendency of the "vertical" velocity
"""
function vertical_velocity!(Diag::DiagnosticVariables{NF},
                            M::PrimitiveEquationModel
                            ) where {NF<:AbstractFloat}


    @unpack dpres_dlon_grid, dpres_dlat_grid = Diag.surface
    @unpack nlev, σ_levels_thick = M.geometry
    @unpack U_mean, V_mean, div_mean = Diag.surface
    
    @boundscheck nlev == Diag.nlev || throw(BoundsError)

    # make sure integration starts with 0
    fill!(Diag.layers[1].dynamics_variables.σ_tend,0)
    fill!(Diag.layers[1].dynamics_variables.σ_m,0)

    @inbounds for k in 1:nlev     # top to bottom, bottom layer separate

        U = Diag.layers[k].grid_variables.U_grid
        V = Diag.layers[k].grid_variables.V_grid
        D = Diag.layers[k].grid_variables.div_grid

        # σ_tend & σ_m sit on half layers below (k+1/2), but its 0 at
        # k=1/2 and nlev+1/2, don't explicitly store k=1/2
        σ_tend = Diag.layers[k].dynamics_variables.σ_tend   # actually on half levels  
        σ_m =  Diag.layers[k].dynamics_variables.σ_tend     # actually on half levels
        uv∇p = Diag.layers[k].dynamics_variables.uv∇p       # on full levels

        # next layer below
        kmax = min(k+1,nlev)    # to avoid access to k = nlev+1
        σ_tend_below = Diag.layers[kmax].dynamics_variables.σ_tend
        σ_m_below = Diag.layers[kmax].dynamics_variables.σ_m

        # TODO check whether coslat unscaling is needed
        for ij in eachgridpoint(U,V,D,U_mean,V_mean,div_mean,dpres_dlon_grid,dpres_dlat_grid)
            uv∇p_ij = (U[ij]-U_mean[ij])*dpres_dlon_grid[ij] + (V[ij]-V_mean[ij])*dpres_dlat_grid[ij]
            uv∇p[ij] = uv∇p_ij
        
            # integration from the top: σ_tend[k] = σ_tend[k-1] - σ_levels_thick...
            # here achieved via -= and the copy into the respective array in the layer below
            σ_tend[ij] -= σ_levels_thick[k]*(uv∇p_ij + D[ij] - div_mean[ij])
            σ_m[ij] -= σ_levels_thick[k]*uv∇p_ij

            # copy into layer below for vertical integration
            # for k = nlev, σ_tend_below == σ_tend, so nothing actually happens here
            σ_tend_below[ij] = σ_tend[ij]
            σ_m_below[ij] = σ_m[ij]
        end
    end
end

function vertical_advection!(   diagn::DiagnosticVariables,
                                model::PrimitiveEquationModel)
    
    @unpack σ_levels_thick⁻¹_half, nlev = model.geometry
    @boundscheck nlev == diagn.nlev || throw(BoundsError)

    # ALL LAYERS (but use indexing tricks to avoid out of bounds access for top/bottom)
    @inbounds for k in 1:nlev       
        # for k==1 "above" term is 0, for k==nlev "below" term is zero
        # avoid out-of-bounds indexing with k_above, k_below as follows
        k_above = k == 1 ? nlev : k-1   # wrap around to access M_nlev+1/2 = 0 (which zeros that term)
        k_below = min(k+1,nlev)         # just saturate, because M_nlev+1/2 = 0 (which zeros that term)
        
        # mass fluxes, M_1/2 = M_nlev+1/2 = 0, but k=1/2 isn't explicitly stored
        σ_tend_above = diagn.layers[k_above].dynamics_variables.σ_tend
        σ_tend_below = diagn.layers[k].dynamics_variables.σ_tend

        # zonal wind
        u_tend = diagn.layers[k].tendencies.u_tend_grid
        U_above = diagn.layers[k_above].grid_variables.U_grid
        U = diagn.layers[k].grid_variables.U_grid
        U_below = diagn.layers[k_below].grid_variables.U_grid

        # meridional wind
        v_tend = diagn.layers[k].tendencies.v_tend_grid
        V_above = diagn.layers[k_above].grid_variables.V_grid
        V = diagn.layers[k].grid_variables.V_grid
        V_below = diagn.layers[k_below].grid_variables.V_grid

        # temperature
        T_tend = diagn.layers[k].tendencies.temp_tend_grid
        T_above = diagn.layers[k_above].grid_variables.temp_grid
        T = diagn.layers[k].grid_variables.temp_grid
        T_below = diagn.layers[k_below].grid_variables.temp_grid

        # humidity
        q_tend = diagn.layers[k].tendencies.humid_tend_grid
        q_above = diagn.layers[k_above].grid_variables.humid_grid
        q = diagn.layers[k].grid_variables.humid_grid
        q_below = diagn.layers[k_below].grid_variables.humid_grid

        Δσk = σ_levels_thick⁻¹_half[k]      # = 1/(2Δσ_k), for convenience

        # TODO check whether coslat unscaling is needed
        @inbounds for ij in eachgridpoint(u_tend,v_tend)
            u_tend[ij] = (σ_tend_above[ij]*(U_above[ij] - U[ij]) +
                            σ_tend_below[ij]*(U[ij] - U_below[ij]))*Δσk
            v_tend[ij] = (σ_tend_above[ij]*(V_above[ij] - V[ij]) +
                            σ_tend_below[ij]*(V[ij] - V_below[ij]))*Δσk
            T_tend[ij] = (σ_tend_above[ij]*(T_above[ij] - T[ij]) +
                            σ_tend_below[ij]*(T[ij] - T_below[ij]))*Δσk
        end

        if model.parameters.dry_core != true    # then also compute vertical advection of humidity
            @inbounds for ij in eachgridpoint(q_tend)
                q_tend[ij] = (σ_tend_above[ij]*(q_above[ij] - q[ij]) +
                                σ_tend_below[ij]*(q[ij] - q_below[ij]))*Δσk
            end
        end
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

function uv_tendencies!(diagn::DiagnosticVariablesLayer,
                        surf::SurfaceVariables,
                        model::PrimitiveEquationModel)
    
    @unpack f_coriolis, coslat⁻² = model.geometry
    @unpack R_dry = model.constants

    @unpack u_tend_grid, v_tend_grid = diagn.tendencies   # already contains vertical advection
    U = diagn.grid_variables.U_grid             # U = u*coslat
    V = diagn.grid_variables.V_grid             # V = v*coslat
    vor = diagn.grid_variables.vor_grid         # relative vorticity
    dpres_dx = surf.dpres_dlon_grid             # zonal gradient of logarithm of surface pressure
    dpres_dy = surf.dpres_dlat_grid             # meridional gradient thereof
    Tᵥ = diagn.grid_variables.temp_virt_grid    # virtual temperature

    # precompute ring indices and boundscheck
    rings = eachring(u_tend_grid,v_tend_grid,U,V,vor,dpres_dx,dpres_dy,Tᵥ)

    @inbounds for (j,ring) in enumerate(rings)
        coslat⁻²j = coslat⁻²[j]
        f = f_coriolis[j]
        for ij in ring
            ω = vor[ij] + f         # absolute vorticity
            RTᵥ = R_dry*Tᵥ[ij]      # gas constant (dry air) times virtual temperature
            u_tend_grid[ij] = (u_tend_grid[ij] + V[ij]*ω - RTᵥ*dpres_dx[ij])*coslat⁻²j
            v_tend_grid[ij] = (v_tend_grid[ij] - U[ij]*ω - RTᵥ*dpres_dy[ij])*coslat⁻²j
        end
    end

    # convert to spectral space
    @unpack u_tend, v_tend = diagn.tendencies   # spectral fields
    S = model.spectral_transform

    spectral!(u_tend,u_tend_grid,S)
    spectral!(v_tend,v_tend_grid,S)

    return nothing
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
function humidity_tendency!(diagn::DiagnosticVariablesLayer,
                            model::PrimitiveEquationModel)

    model.parameters.dry_core && return nothing     # escape immediately for no humidity
    
    @unpack humid_tend, humid_tend_grid = diagn.tendencies
    @unpack div_grid, humid_grid = diagn.grid_variables
    
    # +q*div term of the advection operator
    @inbounds for ij in eachgridpoint(humid_tend_grid,humid_grid,div_grid)
        # add as tend already contains parameterizations + vertical advection
        humid_tend_grid[ij] += humid_grid[ij]*div_grid[ij]
    end

    spectral!(humid_tend,humid_tend_grid,model.spectral_transform)

    # -∇⋅((u,v)*humid)
    @unpack U_grid, V_grid = diagn.grid_variables
    # reuse general work arrays a,b,a_grid,b_grid
    uq = diagn.dynamics_variables.a             # = u*humid in spectral
    vq = diagn.dynamics_variables.b             # = v*humid in spectral
    uq_grid = diagn.dynamics_variables.a_grid   # = u*humid on grid
    vq_grid = diagn.dynamics_variables.b_grid   # = v*humid on grid

    @inbounds for ij in eachgridpoint(humid_tend_grid,humid_grid,U_grid,V_grid)
        uq_grid[ij] = U_grid[ij]*humid_grid[ij]
        vq_grid[ij] = V_grid[ij]*humid_grid[ij]
    end

    spectral!(uq,uq_grid,model.spectral_transform)
    spectral!(vq,vq_grid,model.spectral_transform)

    divergence!(humid_tend,uq,vq,model.spectral_transform,add=true,flipsign=true)
end

"""
    vorticity_flux_divcurl!(    D::DiagnosticVariables{NF}, # all diagnostic variables   
                                G::GeoSpectral{NF}          # struct with geometry and spectral transform
                                ) where {NF<:AbstractFloat}

1) Compute the vorticity advection as the (negative) divergence of the vorticity fluxes -∇⋅(uv*(ζ+f)).
First, compute the uv*(ζ+f), then transform to spectral space and take the divergence and flip the sign.
2) Compute the curl of the vorticity fluxes ∇×(uω,vω) and store in divergence tendency."""
function vorticity_flux_divcurl!(   diagn::DiagnosticVariablesLayer,
                                    G::Geometry,
                                    S::SpectralTransform;
                                    div::Bool=true,         # calculate divergence of vor flux?
                                    curl::Bool=true         # calculate curl of vor flux?
                                    )

    @unpack U_grid, V_grid, vor_grid = diagn.grid_variables
    @unpack vor_tend, div_tend = diagn.tendencies

    uω_coslat⁻¹ = diagn.dynamics_variables.a            # reuse work arrays a,b
    vω_coslat⁻¹ = diagn.dynamics_variables.b
    uω_coslat⁻¹_grid = diagn.dynamics_variables.a_grid
    vω_coslat⁻¹_grid = diagn.dynamics_variables.b_grid

    # STEP 1-3: Abs vorticity, velocity times abs vort
    vorticity_fluxes!(uω_coslat⁻¹_grid,vω_coslat⁻¹_grid,U_grid,V_grid,vor_grid,G)

    spectral!(uω_coslat⁻¹,uω_coslat⁻¹_grid,S)
    spectral!(vω_coslat⁻¹,vω_coslat⁻¹_grid,S)

    # flipsign as RHS is negative ∂ζ/∂t = -∇⋅(uv*(ζ+f)), write directly into tendency
    div && divergence!(vor_tend,uω_coslat⁻¹,vω_coslat⁻¹,S,flipsign=true)

    # = ∇×(uω,vω) = ∇×(uv*(ζ+f)), write directly into tendency
    # curl not needed for BarotropicModel
    curl && curl!(div_tend,uω_coslat⁻¹,vω_coslat⁻¹,S)               
end

"""
    vorticity_fluxes!(  uω_coslat⁻¹::AbstractGrid{NF},      # Output: u*(ζ+f)/coslat
                        vω_coslat⁻¹::AbstractGrid{NF},      # Output: v*(ζ+f)/coslat
                        U::AbstractGrid{NF},                # Input: u*coslat
                        V::AbstractGrid{NF},                # Input: v*coslat
                        vor::AbstractGrid{NF},              # Input: relative vorticity ζ
                        G::Geometry{NF}                     # struct with precomputed geometry arrays
                        ) where {NF<:AbstractFloat}         # number format NF

Compute the vorticity fluxes (u,v)*(ζ+f)/coslat in grid-point space from U,V and vorticity ζ."""
function vorticity_fluxes!( uω_coslat⁻¹::AbstractGrid{NF},  # Output: u*(ζ+f)/coslat
                            vω_coslat⁻¹::AbstractGrid{NF},  # Output: v*(ζ+f)/coslat
                            U::AbstractGrid{NF},            # Input: u*coslat
                            V::AbstractGrid{NF},            # Input: v*coslat
                            vor::AbstractGrid{NF},          # Input: relative vorticity ζ
                            G::Geometry{NF}                 # struct with precomputed geometry arrays
                            ) where {NF<:AbstractFloat}     # number format NF

    nlat = get_nlat(U)
    @unpack f_coriolis, coslat⁻² = G
    @boundscheck length(f_coriolis) == nlat || throw(BoundsError)
    @boundscheck length(coslat⁻²) == nlat || throw(BoundsError)

    rings = eachring(uω_coslat⁻¹,vω_coslat⁻¹,U,V,vor)       # precompute ring indices

    @inbounds for (j,ring) in enumerate(rings)
        coslat⁻²j = coslat⁻²[j]
        f = f_coriolis[j]
        for ij in ring
            # ω = relative vorticity + coriolis and unscale with coslat²
            ω = coslat⁻²j*(vor[ij] + f)
            uω_coslat⁻¹[ij] = ω*U[ij]              # = u(ζ+f)/coslat
            vω_coslat⁻¹[ij] = ω*V[ij]              # = v(ζ+f)/coslat
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
function bernoulli_potential!(  diagn::DiagnosticVariablesLayer,
                                surf::SurfaceVariables,
                                G::Geometry,            
                                S::SpectralTransform,
                                g::Real,                            # gravity
                                )
    
    @unpack U_grid,V_grid = diagn.grid_variables
    @unpack pres_grid = surf
    @unpack bernoulli, bernoulli_grid = diagn.dynamics_variables
    @unpack div_tend = diagn.tendencies

    bernoulli_potential!(bernoulli_grid,U_grid,V_grid,pres_grid,g,G)# = 1/2(u^2 + v^2) + gη on grid
    spectral!(bernoulli,bernoulli_grid,S)                           # to spectral space
    ∇²!(div_tend,bernoulli,S,add=true,flipsign=true)                # add -∇²(1/2(u^2 + v^2) + gη)
end

"""
    bernoulli_potential!(   B::AbstractGrid,    # Output: Bernoulli potential B = 1/2*(u^2+v^2)+g*η
                            U::AbstractGrid,    # zonal velocity *coslat
                            V::AbstractGrid,    # meridional velocity *coslat
                            η::AbstractGrid,    # interface displacement
                            g::Real,            # gravity
                            G::Geometry)

Computes the Bernoulli potential 1/2*(u^2 + v^2) + g*η in grid-point space. This is the
ShallowWater variant that adds the interface displacement η."""
function bernoulli_potential!(  B::AbstractGrid{NF},    # Output: Bernoulli potential B = 1/2*(u^2+v^2)+Φ
                                U::AbstractGrid{NF},    # zonal velocity *coslat
                                V::AbstractGrid{NF},    # meridional velocity *coslat
                                η::AbstractGrid{NF},    # interface displacement
                                g::Real,                # gravity
                                G::Geometry{NF}         # used for precomputed cos²(lat)
                                ) where {NF<:AbstractFloat}
    
    @unpack coslat⁻² = G
    @boundscheck length(coslat⁻²) == get_nlat(U) || throw(BoundsError)

    one_half = convert(NF,0.5)                      # convert to number format NF
    gravity = convert(NF,g)

    rings = eachring(B,U,V,η)

    @inbounds for (j,ring) in enumerate(rings)
        one_half_coslat⁻² = one_half*coslat⁻²[j]
        for ij in ring
            B[ij] = one_half_coslat⁻²*(U[ij]^2 + V[ij]^2) + gravity*η[ij]
        end
    end
end

"""
    bernoulli_potential!(   diagn::DiagnosticVariables, 
                            G::Geometry,
                            S::SpectralTransform)

Computes the Laplace operator ∇² of the Bernoulli potential `B` in spectral space.
    (1) computes the kinetic energy KE=1/2(u^2+v^2) on the grid
    (2) transforms KE to spectral space
    (3) adds geopotential for the bernoulli potential in spectral space
    (4) takes the Laplace operator.
    
This version is used for the PrimitiveEquation model"""
function bernoulli_potential!(  diagn::DiagnosticVariablesLayer,
                                G::Geometry,            
                                S::SpectralTransform,
                                )
    
    @unpack U_grid,V_grid = diagn.grid_variables
    @unpack bernoulli, bernoulli_grid, geopot = diagn.dynamics_variables
    @unpack div_tend = diagn.tendencies

    bernoulli_potential!(bernoulli_grid,U_grid,V_grid,G)    # = 1/2(u^2 + v^2) on grid
    spectral!(bernoulli,bernoulli_grid,S)                   # to spectral space
    add_tendencies!(bernoulli,geopot)                       # add geopotential Φ
    ∇²!(div_tend,bernoulli,S,add=true,flipsign=true)        # add -∇²(1/2(u^2 + v^2) + ϕ)
end

"""
    bernoulli_potential!(   B::AbstractGrid,    # Output: Bernoulli potential B = 1/2*(u^2+v^2)+g*η
                            u::AbstractGrid,    # zonal velocity
                            v::AbstractGrid,    # meridional velocity
                            η::AbstractGrid,    # interface displacement
                            g::Real,            # gravity
                            G::Geometry)

Computes the Bernoulli potential 1/2*(u^2 + v^2), excluding the geopotential, in grid-point space.
This is the PrimitiveEquation-variant where the geopotential is added later in spectral space."""
function bernoulli_potential!(  B::AbstractGrid{NF},    # Output: Bernoulli potential B = 1/2*(u^2+v^2)
                                U::AbstractGrid{NF},    # zonal velocity *coslat
                                V::AbstractGrid{NF},    # meridional velocity *coslat
                                G::Geometry{NF}         # used for precomputed cos²(lat)
                                ) where {NF<:AbstractFloat}
    
    @unpack coslat⁻² = G
    @boundscheck length(coslat⁻²) == get_nlat(U) || throw(BoundsError)

    one_half = convert(NF,0.5)                      # convert to number format NF
    rings = eachring(B,U,V)

    @inbounds for (j,ring) in enumerate(rings)
        one_half_coslat⁻² = one_half*coslat⁻²[j]
        for ij in ring
            B[ij] = one_half_coslat⁻²*(U[ij]^2 + V[ij]^2)
        end
    end
end

function volume_fluxes!(    uh_coslat⁻¹::AbstractGrid{NF},  # Output: zonal volume flux uh/coslat
                            vh_coslat⁻¹::AbstractGrid{NF},  # Output: meridional volume flux vh/coslat
                            U::AbstractGrid{NF},            # U = u*coslat, zonal velocity
                            V::AbstractGrid{NF},            # V = v*coslat, meridional velocity
                            η::AbstractGrid{NF},            # interface displacement
                            orography::AbstractGrid{NF},    # orography
                            H₀::Real,                       # layer thickness at rest
                            G::Geometry{NF},
                            ) where {NF<:AbstractFloat}                                   

    @unpack coslat⁻² = G
    @boundscheck length(coslat⁻²) == get_nlat(η) || throw(BoundsError) 

    H₀ = convert(NF,H₀)

    # compute (uh,vh) on the grid
    # pres_grid is η, the interface displacement
    # layer thickness h = η + H, H is the layer thickness at rest
    # H = H₀ - orography, H₀ is the layer thickness without mountains

    rings = eachring(uh_coslat⁻¹,vh_coslat⁻¹,U,V,η,orography)   # precompute ring indices

    @inbounds for (j,ring) in enumerate(rings)
        coslat⁻²j = coslat⁻²[j]
        for ij in ring
            h = coslat⁻²j*(η[ij] + H₀ - orography[ij])
            uh_coslat⁻¹[ij] = U[ij]*h       # = uh/coslat
            vh_coslat⁻¹[ij] = V[ij]*h       # = vh/coslat
        end
    end
end

"""
    volume_fluxes!( D::DiagnosticVariables{NF},
                    G::Geometry{NF},
                    S::SpectralTransform{NF},
                    B::Boundaries,
                    H₀::Real                    # layer thickness
                    ) where {NF<:AbstractFloat}   

Computes the (negative) divergence of the volume fluxes `uh,vh` for the continuity equation, -∇⋅(uh,vh)"""
function volume_flux_divergence!(   diagn::DiagnosticVariablesLayer,
                                    surface::SurfaceVariables,
                                    G::Geometry,
                                    S::SpectralTransform,
                                    B::Boundaries,              # contains orography
                                    H₀::Real                    # layer thickness
                                    )                           

    @unpack pres_grid, pres_tend = surface
    @unpack U_grid, V_grid = diagn.grid_variables
    @unpack orography = B

    uh_coslat⁻¹ = diagn.dynamics_variables.a            # reuse work arrays a,b
    vh_coslat⁻¹ = diagn.dynamics_variables.b
    uh_coslat⁻¹_grid = diagn.dynamics_variables.a_grid
    vh_coslat⁻¹_grid = diagn.dynamics_variables.b_grid

    volume_fluxes!(uh_coslat⁻¹_grid,vh_coslat⁻¹_grid,U_grid,V_grid,pres_grid,orography,H₀,G)
    
    spectral!(uh_coslat⁻¹,uh_coslat⁻¹_grid,S)
    spectral!(vh_coslat⁻¹,vh_coslat⁻¹_grid,S)

    # compute divergence of volume fluxes and flip sign as ∂η/∂ = -∇⋅(uh,vh)
    divergence!(pres_tend,uh_coslat⁻¹,vh_coslat⁻¹,S,flipsign=true)
end

function interface_relaxation!( η::LowerTriangularMatrix{Complex{NF}},
                                surface::SurfaceVariables{NF},
                                time::DateTime,         # time of relaxation
                                M::ShallowWaterModel,   # contains η⁰, which η is relaxed to
                                ) where NF    

    @unpack pres_tend = surface
    @unpack seasonal_cycle, equinox, tropic_cancer = M.parameters
    A = M.parameters.interface_relax_amplitude

    s = 45/23.5     # heuristic conversion to Legendre polynomials
    θ = seasonal_cycle ? s*tropic_cancer*sin(Dates.days(time - equinox)/365.25*2π) : 0
    η2 = convert(NF,A*(2sind(θ)))           # l=1,m=0 harmonic
    η3 = convert(NF,A*(0.2-1.5cosd(θ)))     # l=2,m=0 harmonic

    τ⁻¹ = inv(M.constants.interface_relax_time)
    pres_tend[2] += τ⁻¹*(η2-η[2])
    pres_tend[3] += τ⁻¹*(η3-η[3])
end

function gridded!(  diagn::DiagnosticVariables,     # all diagnostic variables
                    progn::PrognosticVariables,     # all prognostic variables
                    lf::Int,                        # leapfrog index
                    M::ModelSetup,
                    )

    # all variables on layers
    for (progn_layer,diagn_layer) in zip(progn.layers,diagn.layers)
        gridded!(diagn_layer,progn_layer,lf,M)
    end

    # surface only for ShallowWaterModel or PrimitiveEquationModel
    S = M.spectral_transform
    M isa BarotropicModel ? nothing : gridded!(diagn.surface.pres_grid,progn.pres.leapfrog[lf],S)

    return nothing
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
function gridded!(  diagn::DiagnosticVariablesLayer,   
                    progn::PrognosticVariablesLeapfrog,
                    lf::Int,                            # leapfrog index
                    M::BarotropicModel,
                    )
    
    @unpack vor_grid, U_grid, V_grid = diagn.grid_variables
    u_coslat = diagn.dynamics_variables.a   # use work array a,b for u,v*coslat
    v_coslat = diagn.dynamics_variables.b
    S = M.spectral_transform

    vor_lf = progn.leapfrog[lf].vor     # relative vorticity at leapfrog step lf
    gridded!(vor_grid,vor_lf,S)         # get vorticity on grid from spectral vor
    
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
                lf::Int=1                       # leapfrog index
                M::ShallowWaterModel,           # everything that's constant
                ) where NF

Propagate the spectral state of the prognostic variables `progn` to the
diagnostic variables in `diagn` for the shallow water model. Updates grid vorticity,
grid divergence, grid interface displacement (`pres_grid`) and the velocities
U,V (scaled by cos(lat))."""
function gridded!(  diagn::DiagnosticVariablesLayer,
                    progn::PrognosticVariablesLeapfrog,
                    lf::Int,                            # leapfrog index
                    M::ShallowWaterModel,               # everything that's constant
                    )
    
    @unpack vor_grid, div_grid, U_grid, V_grid = diagn.grid_variables
    u_coslat = diagn.dynamics_variables.a   # use work array a,b for u,v*coslat
    v_coslat = diagn.dynamics_variables.b
    S = M.spectral_transform

    vor_lf = progn.leapfrog[lf].vor     # pick leapfrog index without memory allocation
    div_lf = progn.leapfrog[lf].div   

    # get spectral U,V from vorticity and divergence via stream function Ψ and vel potential ϕ
    # U = u*coslat = -coslat*∂Ψ/∂lat + ∂ϕ/dlon
    # V = v*coslat =  coslat*∂ϕ/∂lat + ∂Ψ/dlon
    UV_from_vordiv!(u_coslat,v_coslat,vor_lf,div_lf,S)

    gridded!(vor_grid,vor_lf,S)         # get vorticity on grid from spectral vor
    gridded!(div_grid,div_lf,S)         # get divergence on grid from spectral div

    # transform to U,V on grid (U,V = u,v*coslat)
    gridded!(U_grid,u_coslat,S)
    gridded!(V_grid,v_coslat,S)

    return nothing
end

function gridded!(  diagn::DiagnosticVariablesLayer,
                    progn::PrognosticVariablesLeapfrog,
                    lf::Int,                            # leapfrog index
                    model::PrimitiveEquationModel,      # everything that's constant
                    )
    
    @unpack vor_grid, div_grid, U_grid, V_grid = diagn.grid_variables
    @unpack temp_grid, humid_grid = diagn.grid_variables
    u_coslat = diagn.dynamics_variables.a       # reuse work arrays a,b
    v_coslat = diagn.dynamics_variables.b

    @unpack dry_core = model.parameters
    S = model.spectral_transform

    vor_lf = progn.leapfrog[lf].vor     # pick leapfrog index without memory allocation
    div_lf = progn.leapfrog[lf].div
    temp_lf = progn.leapfrog[lf].temp
    humid_lf = progn.leapfrog[lf].humid

    # get spectral U,V from vorticity and divergence via stream function Ψ and vel potential ϕ
    # U = u*coslat = -coslat*∂Ψ/∂lat + ∂ϕ/dlon
    # V = v*coslat =  coslat*∂ϕ/∂lat + ∂Ψ/dlon
    UV_from_vordiv!(u_coslat,v_coslat,vor_lf,div_lf,S)

    gridded!(vor_grid,vor_lf,S)         # get vorticity on grid from spectral vor
    gridded!(div_grid,div_lf,S)         # get divergence on grid from spectral div
    gridded!(temp_grid,temp_lf,S)       # (absolute) temperature
    dry_core || gridded!(humid_grid,humid_lf,S)     # specific humidity (wet core only)

    # transform to U,V on grid (U,V = u,v*coslat)
    gridded!(U_grid,u_coslat,S)
    gridded!(V_grid,v_coslat,S)

    return nothing
end