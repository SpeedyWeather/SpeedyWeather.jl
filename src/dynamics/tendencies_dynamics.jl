function pressure_gradients!(   diagn::DiagnosticVariables,
                                progn::PrognosticVariables,
                                lf::Integer,        # leapfrog index
                                S::SpectralTransform)
    
    pres = progn.pres.leapfrog[lf]                      # log of surface pressure
    ∇lnp_x_spec = diagn.layers[1].dynamics_variables.a  # reuse work arrays for gradients
    ∇lnp_y_spec = diagn.layers[1].dynamics_variables.b  # in spectral space
    @unpack ∇lnp_x, ∇lnp_y = diagn.surface              # but store in grid space

    ∇!(∇lnp_x_spec,∇lnp_y_spec,pres,S)                  # CALCULATE ∇ln(pₛ)
    gridded!(∇lnp_x,∇lnp_x_spec,S)                      # transform to grid: zonal gradient
    gridded!(∇lnp_y,∇lnp_y_spec,S)                      # meridional gradient
end

function thickness_weighted_divergence!(diagn::DiagnosticVariablesLayer,
                                        surf::SurfaceVariables,
                                        G::Geometry,
                                        )

    @unpack ∇lnp_x, ∇lnp_y = surf   # zonal, meridional gradient of log surface pressure
    @unpack u_grid, v_grid, div_grid = diagn.grid_variables
    @unpack uv∇lnp, div_weighted = diagn.dynamics_variables
    @unpack coslat⁻¹ = G
    Δσₖ = G.σ_levels_thick[diagn.k]

    rings = eachring(uv∇lnp,u_grid,v_grid,div_grid,∇lnp_x,∇lnp_y)

    @inbounds for (j,ring) in enumerate(rings)
        coslat⁻¹j = coslat⁻¹[j]
        for ij in ring
            uv∇lnp_ij = coslat⁻¹j*(u_grid[ij]*∇lnp_x[ij] + v_grid[ij]*∇lnp_y[ij])
            uv∇lnp[ij] = uv∇lnp_ij
            div_weighted = Δσₖ*(uv∇lnp_ij + div_grid[ij])
        end
    end
end

"""Convert absolute and virtual temperature to anomalies wrt to the reference profile"""
function temperature_anomaly!(  diagn::DiagnosticVariablesLayer,
                                model::PrimitiveEquation)
                    
    Tₖ = model.geometry.temp_ref_profile[diagn.k]   # reference temperature at this level k
    @unpack temp_grid, temp_virt_grid = diagn.grid_variables

    @. temp_grid -= Tₖ          # absolute temperature -> anomaly
    @. temp_virt_grid -= Tₖ     # virtual temperature -> anomaly
end

"""
    vertical_averages!(Diag::DiagnosticVariables,G::Geometry)

Calculates the vertically averaged (weighted by the thickness of the σ level)
velocities (*coslat) and divergence. E.g.

    u_mean = ∑_k=1^nlev Δσ_k * u_k

u,v are averaged in grid-point space, divergence in spectral space.
"""
function vertical_averages!(diagn::DiagnosticVariables{NF},
                            progn::PrognosticVariables{NF},
                            lf::Int,    # leapfrog index for D̄_spec
                            G::Geometry{NF}) where NF
    
    @unpack σ_levels_thick, nlev = G
    ū = diagn.surface.u_mean_grid       # rename for convenience
    v̄ = diagn.surface.v_mean_grid
    D̄ = diagn.surface.div_mean_grid
    D̄_spec = diagn.surface.div_mean

    @boundscheck nlev == diagn.nlev || throw(BoundsError)

    fill!(ū,0)     # reset accumulators from previous vertical average
    fill!(v̄,0)
    fill!(D̄,0)
    fill!(D̄_spec,0)
    # fill!(diagn.layers[1].dynamics_variables.div_sum_above,0)

    @inbounds for k in 1:nlev

        # arrays for layer-thickness weighted column averages
        Δσₖ = σ_levels_thick[k]
        u = diagn.layers[k].grid_variables.u_grid
        v = diagn.layers[k].grid_variables.v_grid
        D = diagn.layers[k].grid_variables.div_grid
        D_spec = progn.layers[k].leapfrog[lf].div
        
        # arrays for sum of divergences for level k from level r=1 to k-1
        k_above = max(1,k-1)
        D_weighted_above = diagn.layers[k_above].dynamics_variables.div_weighted
        D̄ᵣ_above = diagn.layers[k_above].dynamics_variables.div_sum_above
        D̄ᵣ = diagn.layers[k].dynamics_variables.div_sum_above
        
        # u,v,D in grid-point space, with thickness weighting Δσₖ
        @inbounds for ij in eachgridpoint(diagn.surface)
            ū[ij] += u[ij]*Δσₖ
            v̄[ij] += v[ij]*Δσₖ
            D̄[ij] += D[ij]*Δσₖ
            D̄ᵣ[ij] = D̄ᵣ_above[ij] + D_weighted_above[ij]
        end
        
        # above code will incorrectly start the summation at k=1
        k == 1 && fill!(D̄ᵣ,0)   # set to 0 to remove D_weighted at k=1

        # but also divergence in spectral space
        @inbounds for lm in eachharmonic(D̄_spec,D_spec)
            D̄_spec[lm] += D_spec[lm]*Δσₖ
        end
    end
end
        
"""
    surface_pressure_tendency!( Prog::PrognosticVariables,
                                Diag::DiagnosticVariables,
                                lf::Int,
                                M::PrimitiveEquation)

Computes the tendency of the logarithm of surface pressure as

    -(ū*px + v̄*py) - D̄

with ū,v̄ being the vertically averaged velocities; px, py the gradients
of the logarithm of surface pressure ln(p_s) and D̄ the vertically averaged divergence.
1. Calculate ∇ln(p_s) in spectral space, convert to grid.
2. Multiply ū,v̄ with ∇ln(p_s) in grid-point space, convert to spectral.
3. D̄ is subtracted in spectral space.
4. Set tendency of the l=m=0 mode to 0 for better mass conservation."""
function surface_pressure_tendency!(surf::SurfaceVariables{NF},
                                    model::PrimitiveEquation
                                    ) where {NF<:AbstractFloat}

    @unpack pres_tend, pres_tend_grid, ∇lnp_x, ∇lnp_y = surf
    @unpack coslat⁻¹ = model.geometry
    
    # vertical averages need to be computed first!
    ū = surf.u_mean_grid       # rename for convenience
    v̄ = surf.v_mean_grid
    D̄ = surf.div_mean_grid
    D̄_spec = surf.div_mean

    # precompute ring indices
    rings = eachring(pres_tend_grid,∇lnp_x,∇lnp_y,ū,v̄,D̄)

    @inbounds for (j,ring) in enumerate(rings)
        coslat⁻¹j = coslat⁻¹[j]
        for ij in ring
            # -(ū,v̄)⋅∇lnp_s in grid-point space
            # -D̄ is missing here as it's subtracted in spectral space for pres_tend
            # and subtracted in grid-point space for vertical velocity therein
            pres_tend_grid[ij] = -coslat⁻¹j*(ū[ij]*∇lnp_x[ij] +
                                            v̄[ij]*∇lnp_y[ij])
        end
    end

    spectral!(pres_tend,pres_tend_grid,model.spectral_transform)

    # the -D̄ term in spectral
    # for semi-implicit D̄ is calc at time step i-1, i.e. leapfrog lf=1 in vertical_averages!
    pres_tend .-= D̄_spec

    pres_tend[1] = zero(NF)     # for mass conservation
    return nothing
end

function vertical_velocity!(diagn::DiagnosticVariablesLayer,
                            surf::SurfaceVariables,
                            model::PrimitiveEquation)

    @unpack k = diagn                                   # vertical level
    σ̇ = diagn.dynamics_variables.σ_tend                 # vertical mass flux M = pₛσ̇ at k+1/2
    D̄_above = diagn.dynamics_variables.div_sum_above    # sum of thickness-weighted div from level 1:k
    ∂lnpₛ_∂t = surf.pres_tend_grid                      # calculated in surface_pressure_tendency! (excl -D̄)
    D̄ = surf.div_mean_grid                              # vertical avrgd div to be subtract from ∂lnpₛ_∂t
    σk_half = model.geometry.σ_levels_half[k+1]         # σ at k+1/2
    
    # mass flux σ̇ is zero at k=1/2 (not explicitly stored) and k=nlev+1/2 (stored in layer k)
    # set to zero for bottom layer then, and exit immediately
    k == model.geometry.nlev && (fill!(σ̇,0); return nothing)

    # Hoskins and Simmons, 1975 between eq (5) and (6)
    @inbounds for ij in eachgridpoint(σ̇,D̄_above,∂lnpₛ_∂t)
        σ̇[ij] = -D̄_above[ij] - σk_half*(∂lnpₛ_∂t[ij] - D̄[ij])
    end
end

function vertical_advection!(   diagn::DiagnosticVariables,
                                model::PrimitiveEquation)
    
    wet_core = model isa PrimitiveWetCore
    @unpack σ_levels_thick, nlev = model.geometry
    @boundscheck nlev == diagn.nlev || throw(BoundsError)

    # set the k=1 level to zero in the beginning
    u_tend_top = diagn.layers[1].tendencies.u_tend_grid
    v_tend_top = diagn.layers[1].tendencies.v_tend_grid
    temp_tend_top = diagn.layers[1].tendencies.temp_tend_grid
    humid_tend_top = diagn.layers[1].tendencies.humid_tend_grid
    fill!(u_tend_top,0)
    fill!(v_tend_top,0)
    fill!(temp_tend_top,0)
    wet_core && fill!(humid_tend_top,0)

    # ALL LAYERS (but use indexing tricks to avoid out of bounds access for top/bottom)
    @inbounds for k in 1:nlev       
        # for k==1 "above" term is 0, for k==nlev "below" term is zero
        # avoid out-of-bounds indexing with k_above, k_below as follows
        k_below = min(k+1,nlev)         # just saturate, because M_nlev+1/2 = 0 (which zeros that term)
        
        # mass fluxes, M_1/2 = M_nlev+1/2 = 0, but k=1/2 isn't explicitly stored
        σ_tend = diagn.layers[k].dynamics_variables.σ_tend
        
        # layer thickness Δσ on level k
        Δσₖ = σ_levels_thick[k]
        
        u_tend_k = diagn.layers[k].tendencies.u_tend_grid
        u_tend_below = diagn.layers[k_below].tendencies.u_tend_grid
        u = diagn.layers[k].grid_variables.u_grid
        u_below = diagn.layers[k_below].grid_variables.u_grid

        _vertical_advection!(u_tend_below,u_tend_k,σ_tend,u_below,u,Δσₖ)

        v_tend_k = diagn.layers[k].tendencies.v_tend_grid
        v_tend_below = diagn.layers[k_below].tendencies.v_tend_grid
        v = diagn.layers[k].grid_variables.v_grid
        v_below = diagn.layers[k_below].grid_variables.v_grid

        _vertical_advection!(v_tend_below,v_tend_k,σ_tend,v_below,v,Δσₖ)

        T_tend_k = diagn.layers[k].tendencies.temp_tend_grid
        T_tend_below = diagn.layers[k_below].tendencies.temp_tend_grid
        T = diagn.layers[k].grid_variables.temp_grid
        T_below = diagn.layers[k_below].grid_variables.temp_grid

        _vertical_advection!(T_tend_below,T_tend_k,σ_tend,T_below,T,Δσₖ)

        if wet_core
            q_tend_k = diagn.layers[k].tendencies.humid_tend_grid
            q_tend_below = diagn.layers[k_below].tendencies.humid_tend_grid
            q = diagn.layers[k].grid_variables.humid_grid
            q_below = diagn.layers[k_below].grid_variables.humid_grid

            _vertical_advection!(q_tend_below,q_tend_k,σ_tend,q_below,q,Δσₖ)
        end
    end
end

function _vertical_advection!(  ξ_tend_below::Grid,     # tendency of quantity ξ at k+1
                                ξ_tend_k::Grid,         # tendency of quantity ξ at k
                                σ_tend::Grid,           # vertical velocity at k+1/2
                                ξ_below::Grid,          # quantity ξ at k+1
                                ξ::Grid,                # quantity ξ at k
                                Δσₖ::NF                 # layer thickness on σ levels
                                ) where {NF<:AbstractFloat,Grid<:AbstractGrid{NF}}
    Δσₖ2⁻¹ = 1/2Δσₖ                                                 # precompute     
    @inbounds for ij in eachgridpoint(ξ,ξ_tend_k,σ_tend)
        ξ_tend_below[ij] = σ_tend[ij] * (ξ_below[ij] - ξ[ij])       # coslat⁻¹ scaling not here
        ξ_tend_k[ij] = Δσₖ2⁻¹ * (ξ_tend_k[ij] - ξ_tend_below[ij])   # but in vordiv_tendencies!
    end
end


function vordiv_tendencies!(diagn::DiagnosticVariablesLayer,
                            surf::SurfaceVariables,
                            model::PrimitiveEquation)
    
    @unpack f_coriolis, coslat⁻¹, temp_ref_profile = model.geometry
    @unpack R_dry = model.constants

    @unpack u_tend_grid, v_tend_grid = diagn.tendencies   # already contains vertical advection
    u = diagn.grid_variables.u_grid             # velocity
    v = diagn.grid_variables.v_grid             # velocity
    vor = diagn.grid_variables.vor_grid         # relative vorticity
    ∇lnp_x = surf.∇lnp_x                        # zonal gradient of logarithm of surface pressure
    ∇lnp_y = surf.∇lnp_y                        # meridional gradient thereof
    Tᵥ = diagn.grid_variables.temp_virt_grid    # virtual temperature (anomaly!)

    # precompute ring indices and boundscheck
    rings = eachring(u_tend_grid,v_tend_grid,u,v,vor,∇lnp_x,∇lnp_y,Tᵥ)

    @inbounds for (j,ring) in enumerate(rings)
        coslat⁻¹j = coslat⁻¹[j]
        f = f_coriolis[j]
        for ij in ring
            ω = vor[ij] + f                     # absolute vorticity
            RTᵥ = R_dry*Tᵥ[ij]                  # dry gas constant * virtual temperature anomaly
            u_tend_grid[ij] = (u_tend_grid[ij] + v[ij]*ω - RTᵥ*∇lnp_x[ij])*coslat⁻¹j
            v_tend_grid[ij] = (v_tend_grid[ij] - u[ij]*ω - RTᵥ*∇lnp_y[ij])*coslat⁻¹j
        end
    end

    # divergence and curl of that u,v_tend vector for vor,div tendencies
    @unpack vor_tend, div_tend = diagn.tendencies
    u_tend = diagn.dynamics_variables.a
    v_tend = diagn.dynamics_variables.b
    S = model.spectral_transform

    spectral!(u_tend,u_tend_grid,S)
    spectral!(v_tend,v_tend_grid,S)

    curl!(vor_tend,u_tend,v_tend,S)             # ∂ζ/∂t = ∇×(u_tend,v_tend)
    divergence!(div_tend,u_tend,v_tend,S)       # ∂D/∂t = ∇⋅(u_tend,v_tend)
end

"""
Compute the temperature tendency
"""
function temperature_tendency!( diagn::DiagnosticVariablesLayer,
                                surf::SurfaceVariables,
                                model::PrimitiveEquation)

    @unpack temp_tend, temp_tend_grid = diagn.tendencies
    @unpack div_grid, temp_grid = diagn.grid_variables
    @unpack div_sum_above, div_weighted, uv∇lnp = diagn.dynamics_variables
    @unpack κ = model.constants
    Tᵥ = diagn.grid_variables.temp_virt_grid
    Tₖ = model.geometry.temp_ref_profile[diagn.k]

    ∂lnpₛ_∂t = surf.pres_tend_grid                      # calculated in surface_pressure_tendency! (excl -D̄)
    D̄ = surf.div_mean_grid                              # vertical avrgd div to be subtract from ∂lnpₛ_∂t
    
    @unpack k = diagn           # model level 
    σ_lnp_A = model.geometry.σ_lnp_A[k]
    σ_lnp_B = model.geometry.σ_lnp_B[k]

    # +T*div term of the advection operator
    @inbounds for ij in eachgridpoint(temp_tend_grid,temp_grid,div_grid)

        Dlnp_Dt_ij = σ_lnp_A*div_sum_above[ij] + σ_lnp_B*div_weighted[ij] + uv∇lnp[ij]

        # += as tend already contains parameterizations + vertical advection
        temp_tend_grid[ij] += temp_grid[ij]*div_grid[ij] +      # +TD term of hori advection
            κ*(Tᵥ[ij]+Tₖ)*Dlnp_Dt_ij                            # +κTᵥ*Dlnp/Dt, adiabatic term
    end

    spectral!(temp_tend,temp_tend_grid,model.spectral_transform)

    # now add the -∇⋅((u,v)*T') term
    flux_divergence!(temp_tend,temp_grid,diagn,model,add=true,flipsign=true)
end

function humidity_tendency!(diagn::DiagnosticVariablesLayer,
                            model::PrimitiveWetCore)

    @unpack humid_tend, humid_tend_grid = diagn.tendencies
    @unpack humid_grid = diagn.grid_variables

    horizontal_advection!(humid_tend,humid_tend_grid,humid_grid,diagn,model)
end

# no humidity tendency for dry core
humidity_tendency!(::DiagnosticVariablesLayer,::PrimitiveDryCore) = nothing

function horizontal_advection!( A_tend::LowerTriangularMatrix{Complex{NF}}, # Ouput: tendency to write into
                                A_tend_grid::AbstractGrid{NF},              # Input: tendency incl prev terms
                                A_grid::AbstractGrid{NF},                   # Input: grid field to be advected
                                diagn::DiagnosticVariablesLayer{NF},        
                                model::ModelSetup;
                                add::Bool=true) where NF                    # add/overwrite A_tend_grid?

    @unpack div_grid = diagn.grid_variables
    
    @inline kernel(a,b,c) = add ? a+b*c : b*c

    # +A*div term of the advection operator
    @inbounds for ij in eachgridpoint(A_tend_grid,A_grid,div_grid)
        # add as tend already contains parameterizations + vertical advection
        A_tend_grid[ij] = kernel(A_tend_grid[ij],A_grid[ij],div_grid[ij])
    end

    spectral!(A_tend,A_tend_grid,model.spectral_transform)  # for +A*div in spectral space
    
    # now add the -∇⋅((u,v)*A) term
    flux_divergence!(A_tend,A_grid,diagn,model,add=true,flipsign=true)
end

"""Computes -∇⋅((u,v)*A)"""
function flux_divergence!(  A_tend::LowerTriangularMatrix{Complex{NF}}, # Ouput: tendency to write into
                            A_grid::AbstractGrid{NF},                   # Input: grid field to be advected
                            diagn::DiagnosticVariablesLayer{NF},        
                            model::ModelSetup;
                            add::Bool=true,                 # add result to A_tend or overwrite for false
                            flipsign::Bool=true) where NF   # compute -∇⋅((u,v)*A) (true) or ∇⋅((u,v)*A)? 

    @unpack u_grid, v_grid = diagn.grid_variables
    @unpack coslat⁻¹ = model.geometry

    # reuse general work arrays a,b,a_grid,b_grid
    uA = diagn.dynamics_variables.a             # = u*A in spectral
    vA = diagn.dynamics_variables.b             # = v*A in spectral
    uA_grid = diagn.dynamics_variables.a_grid   # = u*A on grid
    vA_grid = diagn.dynamics_variables.b_grid   # = v*A on grid

    rings = eachring(uA_grid,vA_grid,u_grid,v_grid,A_grid)  # precompute ring indices

    @inbounds for (j,ring) in enumerate(rings)
        coslat⁻¹j = coslat⁻¹[j]
        for ij in ring
            Acoslat⁻¹j = A_grid[ij]*coslat⁻¹j
            uA_grid[ij] = u_grid[ij]*Acoslat⁻¹j
            vA_grid[ij] = v_grid[ij]*Acoslat⁻¹j
        end
    end

    spectral!(uA,uA_grid,model.spectral_transform)
    spectral!(vA,vA_grid,model.spectral_transform)

    divergence!(A_tend,uA,vA,model.spectral_transform;add,flipsign)
    return uA,vA        # return for curl calculation (ShallowWater)
end

"""
function vorticity_flux_divcurl!(   diagn::DiagnosticVariablesLayer,
                                    model::ModelSetup;
                                    curl::Bool=true)    # calculate curl of vor flux?

1) Compute the vorticity advection as the (negative) divergence of the vorticity fluxes -∇⋅(uv*(ζ+f)).
First, compute the uv*(ζ+f), then transform to spectral space and take the divergence and flip the sign.
2) Compute the curl of the vorticity fluxes ∇×(uω,vω) and store in divergence tendency."""
function vorticity_flux_divcurl!(   diagn::DiagnosticVariablesLayer,
                                    model::ModelSetup;
                                    curl::Bool=true)    # calculate curl of vor flux?

    G = model.geometry
    S = model.spectral_transform

    @unpack u_grid, v_grid, vor_grid = diagn.grid_variables
    @unpack vor_tend, div_tend = diagn.tendencies

    # add the planetary vorticity f to relative vorticity ζ = absolute vorticity ω
    absolute_vorticity!(vor_grid,G)

    # now do -∇⋅(uω,vω) and store in vor_tend
    uω,vω = flux_divergence!(vor_tend,vor_grid,diagn,model,add=false,flipsign=true)

    # = ∇×(uω,vω) = ∇×(uv*(ζ+f)), write directly into tendency
    # curl not needed for BarotropicModel
    curl && curl!(div_tend,uω,vω,S,add=false,flipsign=false)               
end

function absolute_vorticity!(   vor::AbstractGrid,
                                G::Geometry)

    @unpack f_coriolis = G
    @boundscheck length(f_coriolis) == get_nlat(vor) || throw(BoundsError)

    rings = eachring(vor)
    @inbounds for (j,ring) in enumerate(rings)
        f = f_coriolis[j]
        for ij in ring
            vor[ij] += f
        end
    end
end

"""
    bernoulli_potential!(   diagn::DiagnosticVariablesLayer, 
                            G::Geometry,
                            S::SpectralTransform)

Computes the Laplace operator ∇² of the Bernoulli potential `B` in spectral space.
    (1) computes the kinetic energy KE=1/2(u^2+v^2) on the grid
    (2) transforms KE to spectral space
    (3) adds geopotential for the bernoulli potential in spectral space
    (4) takes the Laplace operator.
    
This version is used for both ShallowWater and PrimitiveEquation, only the geopotential
calculation in geopotential! differs."""
function bernoulli_potential!(  diagn::DiagnosticVariablesLayer{NF},     
                                S::SpectralTransform,
                                ) where NF
    
    @unpack u_grid,v_grid = diagn.grid_variables
    @unpack bernoulli, bernoulli_grid, geopot = diagn.dynamics_variables
    @unpack div_tend = diagn.tendencies
 
    one_half = convert(NF,0.5)
    @. bernoulli_grid = one_half*(u_grid^2 + v_grid^2)      # = 1/2(u^2 + v^2) on grid
    spectral!(bernoulli,bernoulli_grid,S)                   # to spectral space
    bernoulli .+= geopot                                    # add geopotential Φ
    ∇²!(div_tend,bernoulli,S,add=true,flipsign=true)        # add -∇²(1/2(u^2 + v^2) + ϕ)
end

"""
    volume_flux_divergence!(diagn::DiagnosticVariablesLayer,
                            surface::SurfaceVariables,
                            model::ShallowWater)   

Computes the (negative) divergence of the volume fluxes `uh,vh` for the continuity equation, -∇⋅(uh,vh)."""
function volume_flux_divergence!(   diagn::DiagnosticVariablesLayer,
                                    surface::SurfaceVariables,
                                    model::ShallowWater)                        

    @unpack pres_grid, pres_tend = surface
    @unpack orography = model.boundaries.orography
    H₀ = model.constants.layer_thickness

    # compute dynamic layer thickness h on the grid
    # pres_grid is η, the interface displacement, update to
    # layer thickness h = η + H, H is the layer thickness at rest
    # H = H₀ - orography, H₀ is the layer thickness without mountains
    pres_grid .+= H₀ .- orography
    
    # now do -∇⋅(uh,vh) and store in pres_tend
    flux_divergence!(pres_tend,pres_grid,diagn,model,add=false,flipsign=true)
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

function SpeedyTransforms.gridded!( diagn::DiagnosticVariables,     # all diagnostic variables
                                    progn::PrognosticVariables,     # all prognostic variables
                                    lf::Int,                        # leapfrog index
                                    model::ModelSetup,
                                    )

    # all variables on layers
    for (progn_layer,diagn_layer) in zip(progn.layers,diagn.layers)
        gridded!(diagn_layer,progn_layer,lf,model)
    end

    # surface only for ShallowWaterModel or PrimitiveEquation
    S = model.spectral_transform
    model isa Barotropic || gridded!(diagn.surface.pres_grid,progn.pres.leapfrog[lf],S)

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
function SpeedyTransforms.gridded!( diagn::DiagnosticVariablesLayer,   
                                    progn::PrognosticVariablesLeapfrog,
                                    lf::Int,                            # leapfrog index
                                    model::Barotropic,
                                    )
    
    @unpack vor_grid, u_grid, v_grid = diagn.grid_variables
    @unpack U, V = diagn.dynamics_variables
    S = model.spectral_transform

    vor_lf = progn.leapfrog[lf].vor     # relative vorticity at leapfrog step lf
    gridded!(vor_grid,vor_lf,S)         # get vorticity on grid from spectral vor
    
    # get spectral U,V from spectral vorticity via stream function Ψ
    # U = u*coslat = -coslat*∂Ψ/∂lat
    # V = v*coslat = ∂Ψ/∂lon, radius omitted in both cases
    UV_from_vor!(U,V,vor_lf,S)

    # transform from U,V in spectral to u,v on grid (U,V = u,v*coslat)
    gridded!(u_grid,U,S,unscale_coslat=true)
    gridded!(v_grid,V,S,unscale_coslat=true)

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
function SpeedyTransforms.gridded!( diagn::DiagnosticVariablesLayer,
                                    progn::PrognosticVariablesLeapfrog,
                                    lf::Int,                            # leapfrog index
                                    model::ShallowWater,                # everything that's constant
                                    )
    
    @unpack vor_grid, div_grid, u_grid, v_grid = diagn.grid_variables
    @unpack U, V = diagn.dynamics_variables
    S = model.spectral_transform

    vor_lf = progn.leapfrog[lf].vor     # pick leapfrog index without memory allocation
    div_lf = progn.leapfrog[lf].div   

    # get spectral U,V from vorticity and divergence via stream function Ψ and vel potential ϕ
    # U = u*coslat = -coslat*∂Ψ/∂lat + ∂ϕ/dlon
    # V = v*coslat =  coslat*∂ϕ/∂lat + ∂Ψ/dlon
    UV_from_vordiv!(U,V,vor_lf,div_lf,S)

    gridded!(vor_grid,vor_lf,S)         # get vorticity on grid from spectral vor
    gridded!(div_grid,div_lf,S)         # get divergence on grid from spectral div

    # transform from U,V in spectral to u,v on grid (U,V = u,v*coslat)
    gridded!(u_grid,U,S,unscale_coslat=true)
    gridded!(v_grid,V,S,unscale_coslat=true)

    return nothing
end

function SpeedyTransforms.gridded!( diagn::DiagnosticVariablesLayer,
                                    progn::PrognosticVariablesLeapfrog,
                                    lf::Int,                            # leapfrog index
                                    model::PrimitiveEquation,           # everything that's constant
                                    )
    
    @unpack vor_grid, div_grid, u_grid, v_grid = diagn.grid_variables
    @unpack temp_grid, humid_grid = diagn.grid_variables
    @unpack U, V = diagn.dynamics_variables

    S = model.spectral_transform
    wet_core = model isa PrimitiveWetCore

    vor_lf = progn.leapfrog[lf].vor     # pick leapfrog index without memory allocation
    div_lf = progn.leapfrog[lf].div
    temp_lf = progn.leapfrog[lf].temp
    wet_core && (humid_lf = progn.leapfrog[lf].humid)

    # get spectral U,V from vorticity and divergence via stream function Ψ and vel potential ϕ
    # U = u*coslat = -coslat*∂Ψ/∂lat + ∂ϕ/dlon
    # V = v*coslat =  coslat*∂ϕ/∂lat + ∂Ψ/dlon
    UV_from_vordiv!(U,V,vor_lf,div_lf,S)

    gridded!(vor_grid,vor_lf,S)         # get vorticity on grid from spectral vor
    gridded!(div_grid,div_lf,S)         # get divergence on grid from spectral div
    gridded!(temp_grid,temp_lf,S)       # (absolute) temperature
    wet_core && gridded!(humid_grid,humid_lf,S) # specific humidity (wet core only)

    # include humidity effect into temp for everything stability-related
    virtual_temperature!(diagn,temp_lf,model)   # temp = virt temp for dry core

    # transform from U,V in spectral to u,v on grid (U,V = u,v*coslat)
    gridded!(u_grid,U,S,unscale_coslat=true)
    gridded!(v_grid,V,S,unscale_coslat=true)

    return nothing
end