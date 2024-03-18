"""
$(TYPEDSIGNATURES)
Calculate all tendencies for the BarotropicModel."""
function dynamics_tendencies!(  diagn::DiagnosticVariablesLayer,
                                progn::PrognosticVariablesLayer,
                                time::DateTime,
                                model::Barotropic)
    forcing!(diagn, progn, model.forcing, time, model)  # = (Fᵤ, Fᵥ) forcing for u, v
    drag!(diagn, progn, model.drag, time, model)         # drag term for u, v
    vorticity_flux!(diagn, model)                       # = ∇×(v(ζ+f) + Fᵤ, -u(ζ+f) + Fᵥ)
end

"""
$(TYPEDSIGNATURES)
Calculate all tendencies for the ShallowWaterModel."""
function dynamics_tendencies!(  
    diagn::DiagnosticVariablesLayer,
    progn::PrognosticVariablesLayer,
    surface::SurfaceVariables,
    pres::LowerTriangularMatrix,    # spectral pressure/η for geopotential
    time::DateTime,                 # time to evaluate the tendencies at
    model::ShallowWater,
)
    (; forcing, drag, planet, atmosphere, orography) = model
    (; spectral_transform, geometry) = model

    # for compatibility with other ModelSetups pressure pres = interface displacement η here
    forcing!(diagn, progn, forcing, time, model)    # = (Fᵤ, Fᵥ, Fₙ) forcing for u, v, η
    drag!(diagn, progn, drag, time, model)          # drag term for momentum u, v
    
    # = ∇×(v(ζ+f) + Fᵤ, -u(ζ+f) + Fᵥ), tendency for vorticity
    # = ∇⋅(v(ζ+f) + Fᵤ, -u(ζ+f) + Fᵥ), tendency for divergence
    vorticity_flux!(diagn, model)
                                           
    geopotential!(diagn, pres, planet)              # geopotential Φ = gη in shallow water
    bernoulli_potential!(diagn, spectral_transform) # = -∇²(E+Φ), tendency for divergence
    
    # = -∇⋅(uh, vh), tendency for "pressure" η
    volume_flux_divergence!(diagn, surface, orography, atmosphere, geometry, spectral_transform)
end

"""
$(TYPEDSIGNATURES)
Calculate all tendencies for the PrimitiveEquation model (wet or dry)."""
function dynamics_tendencies!(  diagn::DiagnosticVariables,
                                progn::PrognosticVariables,
                                model::PrimitiveEquation,
                                lf::Int=2)          # leapfrog index for tendencies

    O = model.orography
    G = model.geometry
    S = model.spectral_transform
    GP = model.geopotential
    A = model.atmosphere
    I = model.implicit
    (; surface ) = diagn

    # for semi-implicit corrections (α >= 0.5) linear gravity-wave related tendencies are
    # evaluated at previous timestep i-1 (i.e. lf=1 leapfrog time step) 
    # nonlinear terms and parameterizations are always evaluated at lf
    lf_implicit = model.implicit.α == 0 ? lf : 1

    pressure_gradient!(diagn, progn, lf, S)            # calculate ∇ln(pₛ)

    @floop for (diagn_layer, progn_layer) in zip(diagn.layers, progn.layers)
        pressure_flux!(diagn_layer, surface)         # calculate (uₖ, vₖ)⋅∇ln(pₛ)

        # calculate Tᵥ = T + Tₖμq in spectral as a approxmation to Tᵥ = T(1+μq) used for geopotential
        linear_virtual_temperature!(diagn_layer, progn_layer, model, lf_implicit)
        temperature_anomaly!(diagn_layer, I)         # temperature relative to profile
    end

    geopotential!(diagn, GP, O)                       # from ∂Φ/∂ln(pₛ) = -RTᵥ, used in bernoulli_potential!
    vertical_integration!(diagn, progn, lf_implicit, G)# get ū, v̄, D̄ on grid; and and D̄ in spectral
    surface_pressure_tendency!(surface, S)           # ∂ln(pₛ)/∂t = -(ū, v̄)⋅∇ln(pₛ) - D̄

    @floop for layer in diagn.layers
        vertical_velocity!(layer, surface, G)         # calculate σ̇ for the vertical mass flux M = pₛσ̇
                                                    # add the RTₖlnpₛ term to geopotential
        linear_pressure_gradient!(layer, progn.surface, lf_implicit, A, I)
    end                                             # wait all because vertical_velocity! needs to
                                                    # finish before vertical_advection!
    @floop for layer in diagn.layers
        vertical_advection!(layer, diagn, model)      # use σ̇ for the vertical advection of u, v, T, q

        vordiv_tendencies!(layer, surface, model)     # vorticity advection, pressure gradient term
        temperature_tendency!(layer, model)          # hor. advection + adiabatic term
        humidity_tendency!(layer, model)             # horizontal advection of humidity (nothing for wetcore)
        bernoulli_potential!(layer, S)               # add -∇²(E+ϕ+RTₖlnpₛ) term to div tendency
    end
end

"""
$(TYPEDSIGNATURES)
Set the tendencies in `diagn` to zero."""
function zero_tendencies!(diagn::DiagnosticVariables{NF, Grid, Model}) where {NF, Grid, Model<:Barotropic}
    for layer in diagn.layers
        fill!(layer.tendencies.u_tend_grid, 0)
        fill!(layer.tendencies.v_tend_grid, 0)
        fill!(layer.tendencies.vor_tend, 0)
    end
end

"""
$(TYPEDSIGNATURES)
Set the tendencies in `diagn` to zero."""
function zero_tendencies!(diagn::DiagnosticVariables{NF, Grid, Model}) where {NF, Grid, Model<:ShallowWater}
    for layer in diagn.layers
        fill!(layer.tendencies.u_tend_grid, 0)
        fill!(layer.tendencies.v_tend_grid, 0)
        fill!(layer.tendencies.vor_tend, 0)
        fill!(layer.tendencies.div_tend, 0)
    end
    fill!(diagn.surface.pres_tend_grid, 0)
    fill!(diagn.surface.pres_tend, 0)
end

"""
$(TYPEDSIGNATURES)
Set the tendencies in `diagn` to zero."""
function zero_tendencies!(diagn::DiagnosticVariables{NF, Grid, Model}) where {NF, Grid, Model<:PrimitiveDry}
    for layer in diagn.layers
        fill!(layer.tendencies.u_tend_grid, 0)
        fill!(layer.tendencies.v_tend_grid, 0)
        fill!(layer.tendencies.vor_tend, 0)
        fill!(layer.tendencies.div_tend, 0)
        fill!(layer.tendencies.temp_tend_grid, 0)
    end
    fill!(diagn.surface.pres_tend_grid, 0)
    fill!(diagn.surface.pres_tend, 0)
end

"""
$(TYPEDSIGNATURES)
Set the tendencies in `diagn` to zero."""
function zero_tendencies!(diagn::DiagnosticVariables{NF, Grid, Model}) where {NF, Grid, Model<:PrimitiveWet}
    for layer in diagn.layers
        fill!(layer.tendencies.u_tend_grid, 0)
        fill!(layer.tendencies.v_tend_grid, 0)
        fill!(layer.tendencies.vor_tend, 0)
        fill!(layer.tendencies.div_tend, 0)
        fill!(layer.tendencies.temp_tend_grid, 0)
        fill!(layer.tendencies.humid_tend_grid, 0)
    end
    fill!(diagn.surface.pres_tend_grid, 0)
    fill!(diagn.surface.pres_tend, 0)
end

function pressure_gradient!(diagn::DiagnosticVariables,
                            progn::PrognosticVariables,
                            lf::Integer,                   # leapfrog index
                            S::SpectralTransform)
    
    (; pres) = progn.surface.timesteps[lf]                 # log of surface pressure
    ∇lnp_x_spec = diagn.layers[1].dynamics_variables.a     # reuse work arrays for gradients
    ∇lnp_y_spec = diagn.layers[1].dynamics_variables.b     # in spectral space
    (; ∇lnp_x, ∇lnp_y) = diagn.surface                     # but store in grid space

    ∇!(∇lnp_x_spec, ∇lnp_y_spec, pres, S)                  # CALCULATE ∇ln(pₛ)
    gridded!(∇lnp_x, ∇lnp_x_spec, S, unscale_coslat=true)  # transform to grid: zonal gradient
    gridded!(∇lnp_y, ∇lnp_y_spec, S, unscale_coslat=true)  # meridional gradient
end

function pressure_flux!(diagn::DiagnosticVariablesLayer,
                        surf::SurfaceVariables)

    (; ∇lnp_x, ∇lnp_y ) = surf   # zonal, meridional gradient of log surface pressure
    (; u_grid, v_grid ) = diagn.grid_variables
    (; uv∇lnp ) = diagn.dynamics_variables
    
    @. uv∇lnp = u_grid*∇lnp_x + v_grid*∇lnp_y               # the (u, v)⋅∇lnpₛ term
end

"""Convert absolute and virtual temperature to anomalies wrt to the reference profile"""
function temperature_anomaly!(
    diagn::DiagnosticVariablesLayer,
    I::ImplicitPrimitiveEquation,
    )
                    
    Tₖ = I.temp_profile[diagn.k]    # reference temperature on this layer
    (; temp_grid, temp_virt_grid ) = diagn.grid_variables

    @. temp_grid -= Tₖ              # absolute temperature -> anomaly
    @. temp_virt_grid -= Tₖ         # virtual temperature -> anomaly
end

"""
    vertical_integration!(Diag::DiagnosticVariables, G::Geometry)

Calculates the vertically averaged (weighted by the thickness of the σ level)
velocities (*coslat) and divergence. E.g.

    u_mean = ∑_k=1^nlev Δσ_k * u_k

u, v are averaged in grid-point space, divergence in spectral space.
"""
function vertical_integration!( diagn::DiagnosticVariables{NF},
                                progn::PrognosticVariables{NF},
                                lf::Int,    # leapfrog index for D̄_spec
                                G::Geometry{NF}) where NF
    
    (; σ_levels_thick, nlev ) = G
    (; ∇lnp_x, ∇lnp_y ) = diagn.surface  # zonal, meridional grad of log surface pressure
    
    ū = diagn.surface.u_mean_grid           # rename for convenience
    v̄ = diagn.surface.v_mean_grid
    D̄ = diagn.surface.div_mean_grid
    D̄_spec = diagn.surface.div_mean

    @boundscheck nlev == diagn.nlev || throw(BoundsError)

    fill!(ū, 0)     # reset accumulators from previous vertical average
    fill!(v̄, 0)
    fill!(D̄, 0)
    fill!(D̄_spec, 0)

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

        # GRID-POINT SPACE: u, v, D with thickness weighting Δσₖ
        # before this k's u, v, D are added to ū, v̄, D̄ store in the
        # sum_above fields for a 1:k-1 integration
        # which is =0 for k=1 as ū, v̄, D̄ accumulators are 0-initialised
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

with ū, v̄ being the vertically averaged velocities; px, py the gradients
of the logarithm of surface pressure ln(p_s) and D̄ the vertically averaged divergence.
1. Calculate ∇ln(p_s) in spectral space, convert to grid.
2. Multiply ū, v̄ with ∇ln(p_s) in grid-point space, convert to spectral.
3. D̄ is subtracted in spectral space.
4. Set tendency of the l=m=0 mode to 0 for better mass conservation."""
function surface_pressure_tendency!(
    surf::SurfaceVariables,
    S::SpectralTransform,
)
    (; pres_tend, pres_tend_grid, ∇lnp_x, ∇lnp_y ) = surf
    
    # vertical averages need to be computed first!
    ū = surf.u_mean_grid            # rename for convenience
    v̄ = surf.v_mean_grid
    D̄ = surf.div_mean               # spectral

    # in grid-point space the the (ū, v̄)⋅∇lnpₛ term (swap sign in spectral)
    @. pres_tend_grid = ū*∇lnp_x + v̄*∇lnp_y
    spectral!(pres_tend, pres_tend_grid, S)
    
    # for semi-implicit D̄ is calc at time step i-1 in vertical_integration!
    @. pres_tend = -pres_tend - D̄   # the -D̄ term in spectral and swap sign
    
    pres_tend[1] = 0                # for mass conservation
    return nothing
end

function vertical_velocity!(
    diagn::DiagnosticVariablesLayer,
    surf::SurfaceVariables,
    G::Geometry,
)
    (; k) = diagn                                # vertical level
    Δσₖ = G.σ_levels_thick[k]                   # σ level thickness at k
    σk_half = G.σ_levels_half[k+1]              # σ at k+1/2
    σ̇ = diagn.dynamics_variables.σ_tend         # vertical mass flux M = pₛσ̇ at k+1/2
    
                                                # sum of Δσ-weighted div, uv∇lnp from 1:k-1
    (; div_sum_above, uv∇lnp, uv∇lnp_sum_above) = diagn.dynamics_variables
    (; div_grid) = diagn.grid_variables

    ūv̄∇lnp = surf.pres_tend_grid                # calc'd in surface_pressure_tendency! (excl -D̄)
    D̄ = surf.div_mean_grid                      # vertical avrgd div to be added to ūv̄∇lnp

    # mass flux σ̇ is zero at k=1/2 (not explicitly stored) and k=nlev+1/2 (stored in layer k)
    # set to zero for bottom layer then, and exit immediately
    k == G.nlev && (fill!(σ̇, 0); return nothing)

    # Hoskins and Simmons, 1975 just before eq. (6)
    σ̇ .= σk_half*(D̄ .+ ūv̄∇lnp) .-
        (div_sum_above .+ Δσₖ*div_grid) .-      # for 1:k integral add level k σₖ-weighted div 
        (uv∇lnp_sum_above .+ Δσₖ*uv∇lnp)        # and level k σₖ-weighted uv∇lnp here
end

"""
$(TYPEDSIGNATURES)
Function barrier to unpack `model`."""
function vordiv_tendencies!(
    diagn::DiagnosticVariablesLayer,
    surf::SurfaceVariables,
    model::PrimitiveEquation,
)
    (; coriolis, atmosphere, geometry, spectral_transform) = model
    vordiv_tendencies!(diagn, surf, coriolis, atmosphere, geometry, spectral_transform)
end

"""$(TYPEDSIGNATURES)
Tendencies for vorticity and divergence. Excluding Bernoulli potential with geopotential
and linear pressure gradient inside the Laplace operator, which are added later in
spectral space.

    u_tend +=  v*(f+ζ) - RTᵥ'*∇lnp_x
    v_tend += -u*(f+ζ) - RTᵥ'*∇lnp_y

`+=` because the tendencies already contain the parameterizations and vertical advection.
`f` is coriolis, `ζ` relative vorticity, `R` the gas constant `Tᵥ'` the virtual temperature
anomaly, `∇lnp` the gradient of surface pressure and `_x` and `_y` its zonal/meridional
components. The tendencies are then curled/dived to get the tendencies for vorticity/divergence in
spectral space

    ∂ζ/∂t = ∇×(u_tend, v_tend)
    ∂D/∂t = ∇⋅(u_tend, v_tend) + ...

`+ ...` because there's more terms added later for divergence."""
function vordiv_tendencies!(
    diagn::DiagnosticVariablesLayer,
    surf::SurfaceVariables,
    coriolis::AbstractCoriolis,
    atmosphere::AbstractAtmosphere,
    geometry::AbstractGeometry,
    S::SpectralTransform,
)
    (; R_dry) = atmosphere                      # gas constant for dry air
    (; f) = coriolis                            # coriolis parameter
    (; coslat⁻¹) = geometry

    (; u_tend_grid, v_tend_grid) = diagn.tendencies  # already contains vertical advection
    u = diagn.grid_variables.u_grid             # velocity
    v = diagn.grid_variables.v_grid             # velocity
    vor = diagn.grid_variables.vor_grid         # relative vorticity
    ∇lnp_x = surf.∇lnp_x                        # zonal gradient of logarithm of surface pressure
    ∇lnp_y = surf.∇lnp_y                        # meridional gradient thereof
    Tᵥ = diagn.grid_variables.temp_virt_grid    # virtual temperature (anomaly!)

    # precompute ring indices and boundscheck
    rings = eachring(u_tend_grid, v_tend_grid, u, v, vor, ∇lnp_x, ∇lnp_y, Tᵥ)

    @inbounds for (j, ring) in enumerate(rings)
        coslat⁻¹j = coslat⁻¹[j]
        f_j = f[j]
        for ij in ring
            ω = vor[ij] + f_j                   # absolute vorticity
            RTᵥ = R_dry*Tᵥ[ij]                  # dry gas constant * virtual temperature anomaly
            u_tend_grid[ij] = (u_tend_grid[ij] + v[ij]*ω - RTᵥ*∇lnp_x[ij])*coslat⁻¹j
            v_tend_grid[ij] = (v_tend_grid[ij] - u[ij]*ω - RTᵥ*∇lnp_y[ij])*coslat⁻¹j
        end
    end

    # divergence and curl of that u, v_tend vector for vor, div tendencies
    (; vor_tend, div_tend ) = diagn.tendencies
    u_tend = diagn.dynamics_variables.a
    v_tend = diagn.dynamics_variables.b

    spectral!(u_tend, u_tend_grid, S)
    spectral!(v_tend, v_tend_grid, S)

    curl!(vor_tend, u_tend, v_tend, S)         # ∂ζ/∂t = ∇×(u_tend, v_tend)
    divergence!(div_tend, u_tend, v_tend, S)   # ∂D/∂t = ∇⋅(u_tend, v_tend)
    return nothing
end

# function barrier
function tendencies_physics_only!(
    diagn::DiagnosticVariablesLayer,
    model::PrimitiveEquation
)
    wet_core = model isa PrimitiveWet
    tendencies_physics_only!(diagn, model.geometry, model.spectral_transform, wet_core)
end

"""For dynamics=false, after calling parameterization_tendencies! call this function
to transform the physics tendencies from grid-point to spectral space including the
necessary coslat⁻¹ scaling."""
function tendencies_physics_only!(
    diagn::DiagnosticVariablesLayer,
    G::AbstractGeometry,
    S::SpectralTransform,
    wet_core::Bool = true
)
    (; coslat⁻¹) = G

    # already contain parameterizations
    (; u_tend_grid, v_tend_grid, temp_tend_grid, humid_tend_grid) = diagn.tendencies

    # precompute ring indices and boundscheck
    rings = eachring(u_tend_grid, v_tend_grid)

    @inbounds for (j, ring) in enumerate(rings)
        coslat⁻¹j = coslat⁻¹[j]
        for ij in ring
            u_tend_grid[ij] *= coslat⁻¹j
            v_tend_grid[ij] *= coslat⁻¹j
        end
    end

    # divergence and curl of that u, v_tend vector for vor, div tendencies
    (; vor_tend, div_tend, temp_tend, humid_tend) = diagn.tendencies
    u_tend = diagn.dynamics_variables.a
    v_tend = diagn.dynamics_variables.b

    spectral!(u_tend, u_tend_grid, S)
    spectral!(v_tend, v_tend_grid, S)
    spectral!(temp_tend, temp_tend_grid, S)
    wet_core && spectral!(humid_tend, humid_tend_grid, S)

    curl!(vor_tend, u_tend, v_tend, S)         # ∂ζ/∂t = ∇×(u_tend, v_tend)
    divergence!(div_tend, u_tend, v_tend, S)   # ∂D/∂t = ∇⋅(u_tend, v_tend)
    return nothing
end

# function barrier
function temperature_tendency!(
    diagn::DiagnosticVariablesLayer,
    model::PrimitiveEquation,
)
    (; adiabatic_conversion, atmosphere, geometry, spectral_transform, implicit) = model
    temperature_tendency!(diagn, adiabatic_conversion, atmosphere, geometry,
                            spectral_transform, implicit)
end

"""
$(TYPEDSIGNATURES)
Compute the temperature tendency

    ∂T/∂t += -∇⋅((u, v)*T') + T'D + κTᵥ*Dlnp/Dt

`+=` because the tendencies already contain parameterizations and vertical advection.
`T'` is the anomaly with respect to the reference/average temperature. Tᵥ is the virtual
temperature used in the adiabatic term κTᵥ*Dlnp/Dt."""
function temperature_tendency!(
    diagn::DiagnosticVariablesLayer,
    adiabatic_conversion::AbstractAdiabaticConversion,
    atmosphere::AbstractAtmosphere,
    G::Geometry,
    S::SpectralTransform,
    I::ImplicitPrimitiveEquation,
)
    (; temp_tend, temp_tend_grid)  = diagn.tendencies
    (; div_grid, temp_grid) = diagn.grid_variables
    (; uv∇lnp, uv∇lnp_sum_above, div_sum_above) = diagn.dynamics_variables
    
    (; κ) = atmosphere                          # thermodynamic kappa = R_Dry/heat_capacity
    Tᵥ = diagn.grid_variables.temp_virt_grid    # anomaly wrt to Tₖ
    Tₖ = I.temp_profile[diagn.k]                # average layer temperature from reference profile
    
    # coefficients from Simmons and Burridge 1981
    σ_lnp_A = adiabatic_conversion.σ_lnp_A[diagn.k]         # eq. 3.12, -1/Δσₖ*ln(σ_k+1/2/σ_k-1/2)
    σ_lnp_B = adiabatic_conversion.σ_lnp_B[diagn.k]         # eq. 3.12 -αₖ
    
    # semi-implicit: terms here are explicit+implicit evaluated at time step i
    # implicit_correction! then calculated the implicit terms from Vi-1 minus Vi
    # to move the implicit terms to i-1 which is cheaper then the alternative below

    # Adiabatic conversion term following Simmons and Burridge 1981 but for σ coordinates 
    # += as tend already contains parameterizations + vertical advection
    @. temp_tend_grid += temp_grid*div_grid +       # +T'D term of hori advection
        κ*(Tᵥ+Tₖ)*(                                 # +κTᵥ*Dlnp/Dt, adiabatic term
            σ_lnp_A * (div_sum_above+uv∇lnp_sum_above) +    # eq. 3.12 1st term
            σ_lnp_B * (div_grid+uv∇lnp) +                   # eq. 3.12 2nd term
            uv∇lnp)                                         # eq. 3.13

    spectral!(temp_tend, temp_tend_grid, S)

    # now add the -∇⋅((u, v)*T') term
    flux_divergence!(temp_tend, temp_grid, diagn, G, S, add=true, flipsign=true)
    
    return nothing
end

function humidity_tendency!(diagn::DiagnosticVariablesLayer,
                            model::PrimitiveWet)
    G = model.geometry
    S = model.spectral_transform

    (; humid_tend, humid_tend_grid ) = diagn.tendencies
    (; humid_grid ) = diagn.grid_variables

    # add horizontal advection to parameterization + vertical advection tendencies
    horizontal_advection!(humid_tend, humid_tend_grid, humid_grid, diagn, G, S, add=true)
end

# no humidity tendency for dry core
humidity_tendency!(::DiagnosticVariablesLayer, ::PrimitiveDry) = nothing

function horizontal_advection!( 
    A_tend::LowerTriangularMatrix{Complex{NF}}, # Ouput: tendency to write into
    A_tend_grid::AbstractGrid{NF},              # Input: tendency incl prev terms
    A_grid::AbstractGrid{NF},                   # Input: grid field to be advected
    diagn::DiagnosticVariablesLayer,
    G::Geometry,
    S::SpectralTransform;
    add::Bool=true                              # add/overwrite A_tend_grid?
) where NF

    (; div_grid) = diagn.grid_variables
    
    @inline kernel(a, b, c) = add ? a+b*c : b*c

    # +A*div term of the advection operator
    @inbounds for ij in eachgridpoint(A_tend_grid, A_grid, div_grid)
        # add as tend already contains parameterizations + vertical advection
        A_tend_grid[ij] = kernel(A_tend_grid[ij], A_grid[ij], div_grid[ij])
    end

    spectral!(A_tend, A_tend_grid, S)  # for +A*div in spectral space
    
    # now add the -∇⋅((u, v)*A) term
    flux_divergence!(A_tend, A_grid, diagn, G, S, add=true, flipsign=true)
end

"""
$(TYPEDSIGNATURES)
Computes ∇⋅((u, v)*A) with the option to add/overwrite `A_tend` and to
`flip_sign` of the flux divergence by doing so.

- `A_tend =  ∇⋅((u, v)*A)` for `add=false`, `flip_sign=false`
- `A_tend = -∇⋅((u, v)*A)` for `add=false`, `flip_sign=true`
- `A_tend += ∇⋅((u, v)*A)` for `add=true`, `flip_sign=false`
- `A_tend -= ∇⋅((u, v)*A)` for `add=true`, `flip_sign=true`
"""
function flux_divergence!(  A_tend::LowerTriangularMatrix{Complex{NF}}, # Ouput: tendency to write into
                            A_grid::AbstractGrid{NF},                   # Input: grid field to be advected
                            diagn::DiagnosticVariablesLayer,        
                            G::Geometry{NF},
                            S::SpectralTransform{NF};
                            add::Bool=true,                 # add result to A_tend or overwrite for false
                            flipsign::Bool=true) where NF   # compute -∇⋅((u, v)*A) (true) or ∇⋅((u, v)*A)? 

    (; u_grid, v_grid) = diagn.grid_variables
    (; coslat⁻¹) = G

    # reuse general work arrays a, b, a_grid, b_grid
    uA = diagn.dynamics_variables.a             # = u*A in spectral
    vA = diagn.dynamics_variables.b             # = v*A in spectral
    uA_grid = diagn.dynamics_variables.a_grid   # = u*A on grid
    vA_grid = diagn.dynamics_variables.b_grid   # = v*A on grid

    rings = eachring(uA_grid, vA_grid, u_grid, v_grid, A_grid)  # precompute ring indices

    @inbounds for (j, ring) in enumerate(rings)
        coslat⁻¹j = coslat⁻¹[j]
        for ij in ring
            Acoslat⁻¹j = A_grid[ij]*coslat⁻¹j
            uA_grid[ij] = u_grid[ij]*Acoslat⁻¹j
            vA_grid[ij] = v_grid[ij]*Acoslat⁻¹j
        end
    end

    spectral!(uA, uA_grid, S)
    spectral!(vA, vA_grid, S)

    divergence!(A_tend, uA, vA, S; add, flipsign)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Compute the vorticity advection as the curl/div of the vorticity fluxes

`∂ζ/∂t = ∇×(u_tend, v_tend)`
`∂D/∂t = ∇⋅(u_tend, v_tend)`

with

`u_tend = Fᵤ + v*(ζ+f)`
`v_tend = Fᵥ - u*(ζ+f)`

with `Fᵤ, Fᵥ` from `u_tend_grid`/`v_tend_grid` that are assumed to be alread
set in `forcing!`. Set `div=false` for the BarotropicModel which doesn't
require the divergence tendency."""
function vorticity_flux_curldiv!(   diagn::DiagnosticVariablesLayer,
                                    coriolis::AbstractCoriolis,
                                    geometry::Geometry,
                                    S::SpectralTransform;
                                    div::Bool=true,     # also calculate div of vor flux?
                                    add::Bool=false)    # accumulate in vor/div tendencies?
    
    (; f) = coriolis
    (; coslat⁻¹) = geometry

    (; u_tend_grid, v_tend_grid) = diagn.tendencies  # already contains forcing
    u = diagn.grid_variables.u_grid             # velocity
    v = diagn.grid_variables.v_grid             # velocity
    vor = diagn.grid_variables.vor_grid         # relative vorticity

    # precompute ring indices and boundscheck
    rings = eachring(u_tend_grid, v_tend_grid, u, v, vor)

    @inbounds for (j, ring) in enumerate(rings)
        coslat⁻¹j = coslat⁻¹[j]
        f_j = f[j]
        for ij in ring
            ω = vor[ij] + f_j                   # absolute vorticity
            u_tend_grid[ij] = (u_tend_grid[ij] + v[ij]*ω)*coslat⁻¹j
            v_tend_grid[ij] = (v_tend_grid[ij] - u[ij]*ω)*coslat⁻¹j
        end
    end

    # divergence and curl of that u, v_tend vector for vor, div tendencies
    (; vor_tend, div_tend ) = diagn.tendencies
    u_tend = diagn.dynamics_variables.a
    v_tend = diagn.dynamics_variables.b

    spectral!(u_tend, u_tend_grid, S)
    spectral!(v_tend, v_tend_grid, S)

    curl!(vor_tend, u_tend, v_tend, S; add)                 # ∂ζ/∂t = ∇×(u_tend, v_tend)
    div && divergence!(div_tend, u_tend, v_tend, S; add)    # ∂D/∂t = ∇⋅(u_tend, v_tend)
    return nothing       
end

"""
$(TYPEDSIGNATURES)
Vorticity flux tendency in the shallow water equations

`∂ζ/∂t = ∇×(u_tend, v_tend)`
`∂D/∂t = ∇⋅(u_tend, v_tend)`

with

`u_tend = Fᵤ + v*(ζ+f)`
`v_tend = Fᵥ - u*(ζ+f)`

with Fᵤ, Fᵥ the forcing from `forcing!` already in `u_tend_grid`/`v_tend_grid` and
vorticity ζ, coriolis f."""
function vorticity_flux!(diagn::DiagnosticVariablesLayer, model::ShallowWater)
    C = model.coriolis
    G = model.geometry
    S = model.spectral_transform
    vorticity_flux_curldiv!(diagn, C, G, S, div=true, add=true)
end

"""
$(TYPEDSIGNATURES)
Vorticity flux tendency in the barotropic vorticity equation

`∂ζ/∂t = ∇×(u_tend, v_tend)`

with

`u_tend = Fᵤ + v*(ζ+f)`
`v_tend = Fᵥ - u*(ζ+f)`

with Fᵤ, Fᵥ the forcing from `forcing!` already in `u_tend_grid`/`v_tend_grid` and
vorticity ζ, coriolis f."""
function vorticity_flux!(diagn::DiagnosticVariablesLayer, model::Barotropic)
    C = model.coriolis
    G = model.geometry
    S = model.spectral_transform
    vorticity_flux_curldiv!(diagn, C, G, S, div=false, add=true)
end

"""
$(TYPEDSIGNATURES)
Computes the Laplace operator ∇² of the Bernoulli potential `B` in spectral space.
  1. computes the kinetic energy KE = ½(u²+v²) on the grid
  2. transforms KE to spectral space
  3. adds geopotential for the Bernoulli potential in spectral space
  4. takes the Laplace operator.
    
This version is used for both ShallowWater and PrimitiveEquation, only the geopotential
calculation in geopotential! differs."""
function bernoulli_potential!(  diagn::DiagnosticVariablesLayer,     
                                S::SpectralTransform,
                                )
    
    (; u_grid, v_grid ) = diagn.grid_variables
    (; geopot ) = diagn.dynamics_variables
    bernoulli = diagn.dynamics_variables.a                  # reuse work arrays for Bernoulli potential
    bernoulli_grid = diagn.dynamics_variables.a_grid
    (; div_tend ) = diagn.tendencies
 
    half = convert(eltype(u_grid), 0.5)
    @. bernoulli_grid = half*(u_grid^2 + v_grid^2)          # = ½(u² + v²) on grid
    spectral!(bernoulli, bernoulli_grid, S)                 # to spectral space
    bernoulli .+= geopot                                    # add geopotential Φ
    ∇²!(div_tend, bernoulli, S, add=true, flipsign=true)    # add -∇²(½(u² + v²) + ϕ)
end

"""
$(TYPEDSIGNATURES)
Add the linear contribution of the pressure gradient to the geopotential.
The pressure gradient in the divergence equation takes the form

    -∇⋅(Rd * Tᵥ * ∇lnpₛ) = -∇⋅(Rd * Tᵥ' * ∇lnpₛ) - ∇²(Rd * Tₖ * lnpₛ)

So that the second term inside the Laplace operator can be added to the geopotential.
Rd is the gas constant, Tᵥ the virtual temperature and Tᵥ' its anomaly wrt to the
average or reference temperature Tₖ, lnpₛ is the logarithm of surface pressure."""
function linear_pressure_gradient!( 
    diagn::DiagnosticVariablesLayer,
    surface::PrognosticSurfaceTimesteps,
    lf::Int,                # leapfrog index to evaluate tendencies on
    atmosphere::AbstractAtmosphere,
    I::ImplicitPrimitiveEquation,
)                          
    (; R_dry) = atmosphere                   # dry gas constant 
    Tₖ = I.temp_profile[diagn.k]             # reference profile at layer k      
    (; pres) = surface.timesteps[lf]         # logarithm of surface pressure
    (; geopot) = diagn.dynamics_variables

    # -R_dry*Tₖ*∇²lnpₛ, linear part of the ∇⋅RTᵥ∇lnpₛ pressure gradient term
    # Tₖ being the reference temperature profile, the anomaly term T' = Tᵥ - Tₖ is calculated
    # vordiv_tendencies! include as R_dry*Tₖ*lnpₛ into the geopotential on which the operator
    # -∇² is applied in bernoulli_potential!
    @. geopot += R_dry*Tₖ*pres
end

"""
$(TYPEDSIGNATURES)
Computes the (negative) divergence of the volume fluxes `uh, vh` for the continuity equation, -∇⋅(uh, vh)."""
function volume_flux_divergence!(   diagn::DiagnosticVariablesLayer,
                                    surface::SurfaceVariables,
                                    orog::AbstractOrography,
                                    atmosphere::AbstractAtmosphere,
                                    G::AbstractGeometry,
                                    S::SpectralTransform)                        

    (; pres_grid, pres_tend ) = surface
    (; orography ) = orog
    H = atmosphere.layer_thickness

    # compute dynamic layer thickness h on the grid
    # pres_grid is η, the interface displacement, update to
    # layer thickness h = η + H - Hb
    # H is the layer thickness at rest without mountains
    # Hb the orography
    pres_grid .+= H .- orography
    
    # now do -∇⋅(uh, vh) and store in pres_tend
    flux_divergence!(pres_tend, pres_grid, diagn, G, S, add=true, flipsign=true)
end

"""
$(TYPEDSIGNATURES)
Propagate the spectral state of `progn` to `diagn` using time step/leapfrog index `lf`.
Function barrier that calls gridded! for the respective `model`."""
function SpeedyTransforms.gridded!( 
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Int,
    model::ModelSetup)

    # all variables on layers
    for (progn_layer, diagn_layer) in zip(progn.layers, diagn.layers)
        progn_layer_lf = progn_layer.timesteps[lf]
        gridded!(diagn_layer, progn_layer_lf, model)
    end

    # surface only for ShallowWaterModel or PrimitiveEquation
    S = model.spectral_transform
    (; pres_grid) = diagn.surface
    (; pres) = progn.surface.timesteps[lf]
    model isa Barotropic || gridded!(pres_grid, pres, S)

    return nothing
end

"""
$(TYPEDSIGNATURES)
Propagate the spectral state of the prognostic variables `progn` to the
diagnostic variables in `diagn` for the barotropic vorticity model.
Updates grid vorticity, spectral stream function and spectral and grid velocities u, v."""
function SpeedyTransforms.gridded!( diagn::DiagnosticVariablesLayer,   
                                    progn::PrognosticVariablesLayer,
                                    model::Barotropic)
    
    (; vor_grid, u_grid, v_grid ) = diagn.grid_variables
    (; vor ) = progn                    # relative vorticity
    U = diagn.dynamics_variables.a      # reuse work arrays for velocities in spectral
    V = diagn.dynamics_variables.b      # U = u*coslat, V=v*coslat
    S = model.spectral_transform

    gridded!(vor_grid, vor, S)          # get vorticity on grid from spectral vor
    
    # get spectral U, V from spectral vorticity via stream function Ψ
    # U = u*coslat = -coslat*∂Ψ/∂lat
    # V = v*coslat = ∂Ψ/∂lon, radius omitted in both cases
    UV_from_vor!(U, V, vor, S)

    # transform from U, V in spectral to u, v on grid (U, V = u, v*coslat)
    gridded!(u_grid, U, S, unscale_coslat=true)
    gridded!(v_grid, V, S, unscale_coslat=true)
 
    return nothing
end

"""
$(TYPEDSIGNATURES)
Propagate the spectral state of the prognostic variables `progn` to the
diagnostic variables in `diagn` for the shallow water model. Updates grid vorticity,
grid divergence, grid interface displacement (`pres_grid`) and the velocities
u, v."""
function SpeedyTransforms.gridded!( diagn::DiagnosticVariablesLayer,
                                    progn::PrognosticVariablesLayer,
                                    model::ShallowWater,                # everything that's constant
                                    )
    
    (; vor_grid, div_grid, u_grid, v_grid ) = diagn.grid_variables
    (; vor, div) = progn
    U = diagn.dynamics_variables.a      # reuse work arrays for velocities spectral
    V = diagn.dynamics_variables.b      # U = u*coslat, V=v*coslat
    S = model.spectral_transform

    # get spectral U, V from vorticity and divergence via stream function Ψ and vel potential ϕ
    # U = u*coslat = -coslat*∂Ψ/∂lat + ∂ϕ/dlon
    # V = v*coslat =  coslat*∂ϕ/∂lat + ∂Ψ/dlon
    UV_from_vordiv!(U, V, vor, div, S)

    gridded!(vor_grid, vor, S)            # get vorticity on grid from spectral vor
    gridded!(div_grid, div, S)            # get divergence on grid from spectral div

    # transform from U, V in spectral to u, v on grid (U, V = u, v*coslat)
    gridded!(u_grid, U, S, unscale_coslat=true)
    gridded!(v_grid, V, S, unscale_coslat=true)

    return nothing
end

"""
$(TYPEDSIGNATURES)
Propagate the spectral state of the prognostic variables `progn` to the
diagnostic variables in `diagn` for primitive equation models. Updates grid vorticity,
grid divergence, grid temperature, pressure (`pres_grid`) and the velocities
u, v."""
function SpeedyTransforms.gridded!( diagn::DiagnosticVariablesLayer,
                                    progn::PrognosticVariablesLayer,
                                    model::PrimitiveEquation)           # everything that's constant
    
    (; vor_grid, div_grid, u_grid, v_grid ) = diagn.grid_variables
    (; temp_grid, humid_grid ) = diagn.grid_variables
    (; temp_grid_prev, u_grid_prev, v_grid_prev) = diagn.grid_variables
    (; vor, div, temp, humid) = progn
    U = diagn.dynamics_variables.a      # reuse work arrays for velocities spectral
    V = diagn.dynamics_variables.b      # U = u*coslat, V=v*coslat

    # retain previous time step for vertical advection
    @. temp_grid_prev = temp_grid
    @. u_grid_prev = u_grid
    @. v_grid_prev = v_grid

    S = model.spectral_transform
    wet_core = model isa PrimitiveWet

    # get spectral U, V from vorticity and divergence via stream function Ψ and vel potential ϕ
    # U = u*coslat = -coslat*∂Ψ/∂lat + ∂ϕ/dlon
    # V = v*coslat =  coslat*∂ϕ/∂lat + ∂Ψ/dlon
    UV_from_vordiv!(U, V, vor, div, S)

    gridded!(vor_grid, vor, S)                # get vorticity on grid from spectral vor
    gridded!(div_grid, div, S)                # get divergence on grid from spectral div
    gridded!(temp_grid, temp, S)              # (absolute) temperature
    
    if wet_core                             # specific humidity (wet core only)
        gridded!(humid_grid, humid, S)        
        hole_filling!(humid_grid, model.hole_filling, model)  # remove negative humidity
    end

    # include humidity effect into temp for everything stability-related
    temperature_average!(diagn, temp, S)
    virtual_temperature!(diagn, temp, model)  # temp = virt temp for dry core

    # transform from U, V in spectral to u, v on grid (U, V = u, v*coslat)
    gridded!(u_grid, U, S, unscale_coslat=true)
    gridded!(v_grid, V, S, unscale_coslat=true)

    return nothing
end

""" 
$(TYPEDSIGNATURES)
Calculates the average temperature of a layer from the l=m=0 harmonic
and stores the result in `diagn.temp_average`"""
function temperature_average!(
    diagn::DiagnosticVariablesLayer,
    temp::LowerTriangularMatrix,
    S::SpectralTransform,
)
    
    # average from l=m=0 harmonic divided by norm of the sphere
    diagn.temp_average[] = real(temp[1])/S.norm_sphere
end 