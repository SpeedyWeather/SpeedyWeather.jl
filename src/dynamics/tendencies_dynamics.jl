function pressure_gradient!(diagn::DiagnosticVariables,
                            progn::PrognosticVariables,
                            lf::Integer,        # leapfrog index
                            S::SpectralTransform)
    
    pres = progn.pres.timesteps[lf]                      # log of surface pressure
    ∇lnp_x_spec = diagn.layers[1].dynamics_variables.a  # reuse work arrays for gradients
    ∇lnp_y_spec = diagn.layers[1].dynamics_variables.b  # in spectral space
    @unpack ∇lnp_x, ∇lnp_y = diagn.surface              # but store in grid space

    ∇!(∇lnp_x_spec,∇lnp_y_spec,pres,S)                  # CALCULATE ∇ln(pₛ)
    gridded!(∇lnp_x,∇lnp_x_spec,S,unscale_coslat=true)  # transform to grid: zonal gradient
    gridded!(∇lnp_y,∇lnp_y_spec,S,unscale_coslat=true)  # meridional gradient
end

function pressure_flux!(diagn::DiagnosticVariablesLayer,
                        surf::SurfaceVariables)

    @unpack ∇lnp_x, ∇lnp_y = surf   # zonal, meridional gradient of log surface pressure
    @unpack u_grid, v_grid = diagn.grid_variables
    @unpack uv∇lnp = diagn.dynamics_variables
    
    @. uv∇lnp = u_grid*∇lnp_x + v_grid*∇lnp_y           # the (u,v)⋅∇lnpₛ term
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
    vertical_integration!(Diag::DiagnosticVariables,G::Geometry)

Calculates the vertically averaged (weighted by the thickness of the σ level)
velocities (*coslat) and divergence. E.g.

    u_mean = ∑_k=1^nlev Δσ_k * u_k

u,v are averaged in grid-point space, divergence in spectral space.
"""
function vertical_integration!( diagn::DiagnosticVariables{NF},
                                progn::PrognosticVariables{NF},
                                lf::Int,    # leapfrog index for D̄_spec
                                G::Geometry{NF}) where NF
    
    @unpack σ_levels_thick, nlev = G
    @unpack ∇lnp_x, ∇lnp_y = diagn.surface  # zonal, meridional grad of log surface pressure
    
    ū = diagn.surface.u_mean_grid           # rename for convenience
    v̄ = diagn.surface.v_mean_grid
    D̄ = diagn.surface.div_mean_grid
    D̄_spec = diagn.surface.div_mean

    @boundscheck nlev == diagn.nlev || throw(BoundsError)

    fill!(ū,0)     # reset accumulators from previous vertical average
    fill!(v̄,0)
    fill!(D̄,0)
    fill!(D̄_spec,0)

    @inbounds for k in 1:nlev

        # arrays for layer-thickness weighted column averages
        Δσₖ = σ_levels_thick[k]
        u = diagn.layers[k].grid_variables.u_grid
        v = diagn.layers[k].grid_variables.v_grid
        D = diagn.layers[k].grid_variables.div_grid
        D_spec = progn.layers[k].timesteps[lf].div
        
        # for the Σ_r=1^k-1 Δσᵣ(Dᵣ +  u̲⋅∇lnpₛ) vertical integration
        # Simmons and Burridge, 1981 eq 3.12 split into div and u̲⋅∇lnpₛ
        D̄ᵣ = diagn.layers[k].dynamics_variables.div_sum_above
        ūv̄∇lnpᵣ = diagn.layers[k].dynamics_variables.uv∇lnp_sum_above

        # GRID-POINT SPACE: u,v,D with thickness weighting Δσₖ
        # before this k's u,v,D are added to ū,v̄,D̄ store in the
        # sum_above fields for a 1:k-1 integration
        # which is =0 for k=1 as ū,v̄,D̄ accumulators are 0-initialised
        @. D̄ᵣ = D̄
        @. ūv̄∇lnpᵣ = ū*∇lnp_x + v̄*∇lnp_y

        @. ū += u*Δσₖ  # now add the k-th element to the sum
        @. v̄ += v*Δσₖ
        @. D̄ += D*Δσₖ
        
        # SPECTRAL SPACE: divergence
        @. D̄_spec += D_spec*Δσₖ
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
    
    # vertical averages need to be computed first!
    ū = surf.u_mean_grid            # rename for convenience
    v̄ = surf.v_mean_grid
    D̄ = surf.div_mean               # spectral

    # in grid-point space the the (ū,v̄)⋅∇lnpₛ term (swap sign in spectral)
    @. pres_tend_grid = ū*∇lnp_x + v̄*∇lnp_y
    spectral!(pres_tend,pres_tend_grid,model.spectral_transform)
    
    # for semi-implicit D̄ is calc at time step i-1 in vertical_integration!
    @. pres_tend = -pres_tend - D̄   # the -D̄ term in spectral and swap sign
    
    pres_tend[1] = zero(NF)         # for mass conservation
    spectral_truncation!(pres_tend) # remove lmax+1 row, only vectors use it
    return nothing
end

function vertical_velocity!(diagn::DiagnosticVariablesLayer,
                            surf::SurfaceVariables,
                            model::PrimitiveEquation)

    @unpack k = diagn                           # vertical level
    Δσₖ = model.geometry.σ_levels_thick[k]      # σ level thickness at k
    σk_half = model.geometry.σ_levels_half[k+1] # σ at k+1/2
    σ̇ = diagn.dynamics_variables.σ_tend         # vertical mass flux M = pₛσ̇ at k+1/2
    
                                                # sum of Δσ-weighted div, uv∇lnp from 1:k-1
    @unpack div_sum_above, uv∇lnp, uv∇lnp_sum_above = diagn.dynamics_variables
    @unpack div_grid = diagn.grid_variables

    ūv̄∇lnp = surf.pres_tend_grid                # calc'd in surface_pressure_tendency! (excl -D̄)
    D̄ = surf.div_mean_grid                      # vertical avrgd div to be added to ūv̄∇lnp

    # mass flux σ̇ is zero at k=1/2 (not explicitly stored) and k=nlev+1/2 (stored in layer k)
    # set to zero for bottom layer then, and exit immediately
    k == model.geometry.nlev && (fill!(σ̇,0); return nothing)

    # Hoskins and Simmons, 1975 just before eq. (6)
    σ̇ .= σk_half*(D̄ .+ ūv̄∇lnp) .-
        (div_sum_above .+ Δσₖ*div_grid) .-      # for 1:k integral add level k σₖ-weighted div 
        (uv∇lnp_sum_above .+ Δσₖ*uv∇lnp)        # and level k σₖ-weighted uv∇lnp here
end

# MULTI LAYER VERSION
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

# SINGLE LAYER VERSION
function vertical_advection!(   layer::DiagnosticVariablesLayer,
                                diagn::DiagnosticVariables,
                                model::PrimitiveEquation)
            
    @unpack k = layer       # which layer are we on?
    
    wet_core = model isa PrimitiveWetCore
    @unpack σ_levels_thick, nlev = model.geometry
     
    # for k==1 "above" term is 0, for k==nlev "below" term is zero
    # avoid out-of-bounds indexing with k_above, k_below as follows
    k_above = max(1,k-1)        # just saturate, because M_1/2 = 0 (which zeros that term)
    k_below = min(k+1,nlev)     # just saturate, because M_nlev+1/2 = 0 (which zeros that term)
        
    # mass fluxes, M_1/2 = M_nlev+1/2 = 0, but k=1/2 isn't explicitly stored
    σ_tend_above = diagn.layers[k_above].dynamics_variables.σ_tend
    σ_tend_below = layer.dynamics_variables.σ_tend
    
    # layer thickness Δσ on level k
    Δσₖ = σ_levels_thick[k]
    
    u_tend = layer.tendencies.u_tend_grid                   # zonal velocity
    u_above = diagn.layers[k_above].grid_variables.u_grid
    u = layer.grid_variables.u_grid
    u_below = diagn.layers[k_below].grid_variables.u_grid

    _vertical_advection!(u_tend,σ_tend_above,σ_tend_below,u_above,u,u_below,Δσₖ)

    v_tend = layer.tendencies.v_tend_grid                   # meridional velocity
    v_above = diagn.layers[k_above].grid_variables.v_grid
    v = layer.grid_variables.v_grid
    v_below = diagn.layers[k_below].grid_variables.v_grid

    _vertical_advection!(v_tend,σ_tend_above,σ_tend_below,v_above,v,v_below,Δσₖ)

    T_tend = layer.tendencies.temp_tend_grid                   # temperature
    T_above = diagn.layers[k_above].grid_variables.temp_grid
    T = layer.grid_variables.temp_grid
    T_below = diagn.layers[k_below].grid_variables.temp_grid

    _vertical_advection!(T_tend,σ_tend_above,σ_tend_below,T_above,T,T_below,Δσₖ)

    if wet_core
        q_tend = layer.tendencies.humid_tend_grid               # humidity
        q_above = diagn.layers[k_above].grid_variables.humid_grid
        q = layer.grid_variables.humid_grid
        q_below = diagn.layers[k_below].grid_variables.humid_grid
    
        _vertical_advection!(q_tend,σ_tend_above,σ_tend_below,q_above,q,q_below,Δσₖ)
    end
end

# SINGLE THREADED VERSION uses the layer below to store intermediate result
function _vertical_advection!(  ξ_tend_below::Grid,     # tendency of quantity ξ at k+1
                                ξ_tend_k::Grid,         # tendency of quantity ξ at k
                                σ_tend::Grid,           # vertical velocity at k+1/2
                                ξ_below::Grid,          # quantity ξ at k+1
                                ξ::Grid,                # quantity ξ at k
                                Δσₖ::NF                 # layer thickness on σ levels
                                ) where {NF<:AbstractFloat,Grid<:AbstractGrid{NF}}
    Δσₖ2⁻¹ = -1/2Δσₖ                                    # precompute
    @. ξ_tend_below = σ_tend * (ξ_below - ξ)            # store without Δσ-scaling in layer below
    @. ξ_tend_k = Δσₖ2⁻¹ * (ξ_tend_k + ξ_tend_below)    # combine with layer above and scale
end

# MULTI THREADED VERSION only writes into layer k
function _vertical_advection!(  ξ_tend::Grid,           # tendency of quantity ξ at k
                                σ_tend_above::Grid,     # vertical velocity at k-1/2
                                σ_tend_below::Grid,     # vertical velocity at k+1/2
                                ξ_above::Grid,          # quantity ξ at k-1
                                ξ::Grid,                # quantity ξ at k
                                ξ_below::Grid,          # quantity ξ at k+1
                                Δσₖ::NF                 # layer thickness on σ levels
                                ) where {NF<:AbstractFloat,Grid<:AbstractGrid{NF}}
    Δσₖ2⁻¹ = -1/2Δσₖ                                    # precompute

    # += as the tendencies already contain the parameterizations
    @. ξ_tend += Δσₖ2⁻¹ * (σ_tend_below*(ξ_below - ξ) + σ_tend_above*(ξ - ξ_above))
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

    curl!(vor_tend,u_tend,v_tend,S)         # ∂ζ/∂t = ∇×(u_tend,v_tend)
    divergence!(div_tend,u_tend,v_tend,S)   # ∂D/∂t = ∇⋅(u_tend,v_tend)

    # only vectors make use of the lmax+1 row, set to zero for scalars
    spectral_truncation!(vor_tend)           
    spectral_truncation!(div_tend)
end

"""
Compute the temperature tendency
"""
function temperature_tendency!( diagn::DiagnosticVariablesLayer,
                                surf::SurfaceVariables,
                                model::PrimitiveEquation)

    @unpack temp_tend, temp_tend_grid = diagn.tendencies
    @unpack div_grid, temp_grid = diagn.grid_variables
    @unpack uv∇lnp, uv∇lnp_sum_above, div_sum_above = diagn.dynamics_variables
    @unpack κ = model.constants                 # thermodynamic kappa
    @unpack k = diagn                           # model level
    @unpack nlev = model.geometry
    
    Tᵥ = diagn.grid_variables.temp_virt_grid    # anomaly wrt to Tₖ
    Tₖs = model.geometry.temp_ref_profile       # reference temperatures
    Tₖ = Tₖs[k]                                 # reference temperature at k
    
    # coefficients from Simmons and Burridge 1981
    σ_lnp_A = model.geometry.σ_lnp_A[k]         # eq. 3.12, -1/Δσₖ*ln(σ_k+1/2/σ_k-1/2)
    σ_lnp_B = model.geometry.σ_lnp_B[k]         # eq. 3.12 -αₖ
    
    # semi-implicit: terms here are explicit+implicit evaluated at time step i
    # implicit_correction! then calculated the implicit terms from Vi-1 minus Vi
    # to move the implicit terms to i-1 which is cheaper then the alternative below

    # Adiabatic term following Simmons and Burridge 1981 but for σ coordinates 
    # += as tend already contains parameterizations + vertical advection
    @. temp_tend_grid += temp_grid*div_grid +       # +T'D term of hori advection
        κ*(Tᵥ+Tₖ)*(                                 # +κTᵥ*Dlnp/Dt, adiabatic term
            σ_lnp_A * (div_sum_above+uv∇lnp_sum_above) +    # eq. 3.12 1st term
            σ_lnp_B * (div_grid+uv∇lnp) +                   # eq. 3.12 2nd term
            uv∇lnp)                                         # eq. 3.13
    
    # SEMI-IMPLICIT ALTERNATIVE
    # evaluate only the explicit terms at time step i and the implicit terms
    # in implicit_correction! at i-1, however, this is more expensive then above

    # ūv̄∇lnp = surf.pres_tend_grid
    # # for explicit vertical advection get reference temperature differences
    # k_above = max(1,k-1)                        # layer index above
    # k_below = min(k+1,nlev)                     # layer index below
    # ΔT_above = Tₖ - Tₖs[k_above]                # temperature difference to layer above
    # ΔT_below = Tₖs[k_below] - Tₖ                # and to layer below

    # # for explicit vertical advection terms
    # σₖ = model.geometry.σ_levels_full[k]        # should be Σ_r=1^k Δσᵣ for model top at >0hPa
    # σₖ_above = model.geometry.σ_levels_full[k_above]
    # Δσₖ2⁻¹ = -1/2model.geometry.σ_levels_thick[k]

    # # Hoskins and Simmons 1975, Appendix 1 but the adiabatic term therein as above
    # @. temp_tend_grid += temp_grid*div_grid +
    #     Δσₖ2⁻¹*ΔT_below*(      σₖ*ūv̄∇lnp - (uv∇lnp_sum_above + σₖ*uv∇lnp)) +
    #     Δσₖ2⁻¹*ΔT_above*(σₖ_above*ūv̄∇lnp -  uv∇lnp_sum_above) + 
    #     κ*Tₖ*(σ_lnp_A*uv∇lnp_sum_above + σ_lnp_B*uv∇lnp) + 
    #     κ*Tᵥ*(σ_lnp_A * (div_sum_above+uv∇lnp_sum_above) + σ_lnp_B * (div_grid+uv∇lnp)) +
    #     κ*(Tᵥ+Tₖ)*uv∇lnp

    spectral!(temp_tend,temp_tend_grid,model.spectral_transform)

    # now add the -∇⋅((u,v)*T') term
    flux_divergence!(temp_tend,temp_grid,diagn,model,add=true,flipsign=true)
    
    # only vectors make use of the lmax+1 row, set to zero for scalars
    spectral_truncation!(temp_tend) 
end

function humidity_tendency!(diagn::DiagnosticVariablesLayer,
                            model::PrimitiveWetCore)

    @unpack humid_tend, humid_tend_grid = diagn.tendencies
    @unpack humid_grid = diagn.grid_variables

    horizontal_advection!(humid_tend,humid_tend_grid,humid_grid,diagn,model)

    # only vectors make use of the lmax+1 row, set to zero for scalars
    spectral_truncation!(humid_tend) 
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
    @unpack geopot = diagn.dynamics_variables
    bernoulli = diagn.dynamics_variables.a                  # reuse work arrays for bernoulli potential
    bernoulli_grid = diagn.dynamics_variables.a_grid
    @unpack div_tend = diagn.tendencies
 
    half = convert(NF,0.5)
    @. bernoulli_grid = half*(u_grid^2 + v_grid^2)          # = 1/2(u^2 + v^2) on grid
    spectral!(bernoulli,bernoulli_grid,S)                   # to spectral space
    bernoulli .+= geopot                                    # add geopotential Φ
    ∇²!(div_tend,bernoulli,S,add=true,flipsign=true)        # add -∇²(1/2(u^2 + v^2) + ϕ)
end

function linear_pressure_gradient!( diagn::DiagnosticVariablesLayer,
                                    progn::PrognosticVariables,
                                    model::PrimitiveEquation,
                                    lf::Int)                # leapfrog index to evaluate tendencies on
    
    @unpack R_dry = model.constants
    Tₖ = model.geometry.temp_ref_profile[diagn.k]           # reference temperature at layer k      
    pres = progn.pres.timesteps[lf]
    @unpack geopot = diagn.dynamics_variables

    # -R_dry*Tₖ*∇²lnpₛ, linear part of the ∇⋅RTᵥ∇lnpₛ pressure gradient term
    # Tₖ being the reference temperature profile, the anomaly term T' = Tᵥ - Tₖ is calculated
    # vordiv_tendencies! include as R_dry*Tₖ*lnpₛ into the geopotential on which the operator
    # -∇² is applied in bernoulli_potential!
    @. geopot += R_dry*Tₖ*pres
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
    @unpack seasonal_cycle, equinox, axial_tilt = M.parameters.planet
    A = M.parameters.interface_relax_amplitude

    s = 45/23.5     # heuristic conversion to Legendre polynomials
    θ = seasonal_cycle ? s*axial_tilt*sin(Dates.days(time - equinox)/365.25*2π) : 0
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
    model isa Barotropic || gridded!(diagn.surface.pres_grid,progn.pres.timesteps[lf],S)

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
                                    progn::PrognosticLayerTimesteps,
                                    lf::Int,                            # leapfrog index
                                    model::Barotropic)
    
    @unpack vor_grid, u_grid, v_grid = diagn.grid_variables
    U = diagn.dynamics_variables.a      # reuse work arrays for velocities in spectral
    V = diagn.dynamics_variables.b      # U = u*coslat, V=v*coslat
    S = model.spectral_transform

    vor_lf = progn.timesteps[lf].vor     # relative vorticity at leapfrog step lf
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
                                    progn::PrognosticLayerTimesteps,
                                    lf::Int,                            # leapfrog index
                                    model::ShallowWater,                # everything that's constant
                                    )
    
    @unpack vor_grid, div_grid, u_grid, v_grid = diagn.grid_variables
    U = diagn.dynamics_variables.a      # reuse work arrays for velocities spectral
    V = diagn.dynamics_variables.b      # U = u*coslat, V=v*coslat
    S = model.spectral_transform

    vor_lf = progn.timesteps[lf].vor     # pick leapfrog index without memory allocation
    div_lf = progn.timesteps[lf].div   

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
                                    progn::PrognosticLayerTimesteps,
                                    lf::Int,                            # leapfrog index
                                    model::PrimitiveEquation,           # everything that's constant
                                    )
    
    @unpack vor_grid, div_grid, u_grid, v_grid = diagn.grid_variables
    @unpack temp_grid, humid_grid = diagn.grid_variables
    U = diagn.dynamics_variables.a      # reuse work arrays for velocities spectral
    V = diagn.dynamics_variables.b      # U = u*coslat, V=v*coslat

    S = model.spectral_transform
    wet_core = model isa PrimitiveWetCore

    vor_lf = progn.timesteps[lf].vor     # pick leapfrog index without memory allocation
    div_lf = progn.timesteps[lf].div
    temp_lf = progn.timesteps[lf].temp
    wet_core && (humid_lf = progn.timesteps[lf].humid)

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