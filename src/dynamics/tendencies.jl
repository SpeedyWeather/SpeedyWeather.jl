"""$(TYPEDSIGNATURES)
Calculate all tendencies for the BarotropicModel."""
function dynamics_tendencies!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Integer,                    # leapfrog index to evaluate tendencies at
    model::Barotropic,
)
    forcing!(diagn, progn, lf, model)   # = (Fᵤ, Fᵥ) forcing for u, v
    drag!(diagn, progn, lf, model)      # drag term for u, v
    vorticity_flux!(diagn, model)       # = ∇×(v(ζ+f) + Fᵤ, -u(ζ+f) + Fᵥ)
    tracer_advection!(diagn, model)
end

"""
$(TYPEDSIGNATURES)
Calculate all tendencies for the ShallowWaterModel."""
function dynamics_tendencies!(  
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Integer,                    # leapfrog index to evaluate tendencies at
    model::ShallowWater,
)
    (; planet, atmosphere, orography) = model
    (; spectral_transform, geometry) = model

    # for compatibility with other AbstractModels pressure pres = interface displacement η here
    forcing!(diagn, progn, lf, model)   # = (Fᵤ, Fᵥ, Fₙ) forcing for u, v, η
    drag!(diagn, progn, lf, model)      # drag term for u, v

    # = ∇×(v(ζ+f) + Fᵤ, -u(ζ+f) + Fᵥ), tendency for vorticity
    # = ∇⋅(v(ζ+f) + Fᵤ, -u(ζ+f) + Fᵥ), tendency for divergence
    vorticity_flux!(diagn, model)

    pres_lf = get_step(progn.pres, lf)              # interface displacement η at leapfrog step lf
    geopotential!(diagn, pres_lf, planet)           # geopotential Φ = gη in shallow water
    bernoulli_potential!(diagn, spectral_transform) # = -∇²(E+Φ), tendency for divergence
    
    # = -∇⋅(uh, vh), tendency for "pressure" η
    volume_flux_divergence!(diagn, orography, atmosphere, geometry, spectral_transform)
    
    # advect all tracers
    tracer_advection!(diagn, model)
end

"""$(TYPEDSIGNATURES)
Calculate all tendencies for the PrimitiveEquation model (wet or dry)."""
function dynamics_tendencies!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Integer,                # leapfrog index for tendencies
    model::PrimitiveEquation,
)

    (; orography, geometry, spectral_transform, geopotential, atmosphere, implicit) = model

    # for semi-implicit corrections (α >= 0.5) linear gravity-wave related tendencies are
    # evaluated at previous timestep i-1 (i.e. lf=1 leapfrog time step) 
    # nonlinear terms and parameterizations are always evaluated at lf
    lf_implicit = implicit.α == 0 ? lf : 1

    # calculate ∇ln(pₛ), then (u_k, v_k)⋅∇ln(p_s)
    pressure_gradient_flux!(diagn, progn, lf, spectral_transform)

    # calculate Tᵥ = T + Tₖμq in spectral as a approxmation to Tᵥ = T(1+μq) used for geopotential
    linear_virtual_temperature!(diagn, progn, lf_implicit, model)

    # temperature relative to profile
    temperature_anomaly!(diagn, implicit)

    # from ∂Φ/∂ln(pₛ) = -RTᵥ for bernoulli_potential!
    geopotential!(diagn, geopotential, orography)

    # get ū, v̄, D̄ on grid; D̄ in spectral
    vertical_integration!(diagn, progn, lf_implicit, geometry)

    # ∂ln(pₛ)/∂t = -(ū, v̄)⋅∇ln(pₛ) - D̄
    surface_pressure_tendency!(diagn, spectral_transform)

    # calculate vertical velocity σ̇ in sigma coordinates for the vertical mass flux M = p_s*σ̇
    vertical_velocity!(diagn, geometry)       
    
    # add the RTₖlnpₛ term to geopotential
    linear_pressure_gradient!(diagn, progn, lf_implicit, atmosphere, implicit)

    # use σ̇ for the vertical advection of u, v, T, q
    vertical_advection!(diagn, model)

    # vorticity advection, pressure gradient term
    vordiv_tendencies!(diagn, model)

    # hor. advection + adiabatic term
    temperature_tendency!(diagn, model)

    # horizontal advection of humidity (nothing for wetcore)
    humidity_tendency!(diagn, model)

    # add -∇²(E+ϕ+RTₖlnpₛ) term to div tendency
    bernoulli_potential!(diagn, spectral_transform)

    # advect all tracers
    tracer_advection!(diagn, model)

    return nothing
end

"""$(TYPEDSIGNATURES)
Compute the gradient ∇lnp_s of the logarithm of surface pressure, followed by
its flux, (u,v) * ∇lnp_s."""
function pressure_gradient_flux!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Integer,                   # leapfrog index
    S::SpectralTransform,
)

    (; scratch_memory) = diagn.dynamics

    # PRESSURE GRADIENT
    pres = get_step(progn.pres, lf)         # log of surface pressure at leapfrog step lf
    ∇lnp_x_spec = diagn.dynamics.a_2D       # reuse 2D work arrays for gradients
    ∇lnp_y_spec = diagn.dynamics.b_2D       # in spectral space
    (; ∇lnp_x, ∇lnp_y) = diagn.dynamics     # but store in grid space

    ∇!(∇lnp_x_spec, ∇lnp_y_spec, pres, S)                   # CALCULATE ∇ln(pₛ)
    transform!(∇lnp_x, ∇lnp_x_spec, scratch_memory, S, unscale_coslat=true) # transform to grid: zonal gradient
    transform!(∇lnp_y, ∇lnp_y_spec, scratch_memory, S, unscale_coslat=true) # meridional gradient

    (; u_grid, v_grid ) = diagn.grid
    (; uv∇lnp ) = diagn.dynamics

    # PRESSURE GRADIENT FLUX
    @inbounds for k in eachlayer(u_grid, v_grid, uv∇lnp)
        for ij in eachgridpoint(u_grid, v_grid, uv∇lnp)
            # the (u, v)⋅∇lnp_s term
            uv∇lnp[ij, k] = u_grid[ij, k]*∇lnp_x[ij] + v_grid[ij, k]*∇lnp_y[ij]
        end
    end
end

"""$(TYPEDSIGNATURES)
Convert absolute and virtual temperature to anomalies wrt to the reference profile"""
function temperature_anomaly!(
    diagn::DiagnosticVariables,
    implicit::ImplicitPrimitiveEquation,
)                    
    (; temp_profile) = implicit     # reference temperature profile
    (; temp_grid, temp_virt_grid ) = diagn.grid

    @inbounds for k in eachlayer(temp_grid, temp_virt_grid)
        Tₖ = temp_profile[k]
        for ij in eachgridpoint(temp_grid, temp_virt_grid)
            temp_grid[ij, k]      -= Tₖ  # absolute temperature -> anomaly
            temp_virt_grid[ij, k] -= Tₖ  # virtual temperature -> anomaly
        end
    end
end

"""$(TYPEDSIGNATURES)
Calculates the vertically averaged (weighted by the thickness of the σ level)
velocities (*coslat) and divergence. E.g.

    u_mean = ∑_k=1^nlayers Δσ_k * u_k

u, v are averaged in grid-point space, divergence in spectral space.
"""
function vertical_integration!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Integer,                    # leapfrog index for D̄_spec
    geometry::Geometry,
)
    (; σ_levels_thick, nlayers ) = geometry
    (; ∇lnp_x, ∇lnp_y ) = diagn.dynamics    # zonal, meridional grad of log surface pressure
    (; u_grid, v_grid, div_grid) = diagn.grid
    (; u_mean_grid, v_mean_grid, div_mean_grid, div_mean) = diagn.dynamics
    (; div_sum_above, uv∇lnp_sum_above) = diagn.dynamics
    div = get_step(progn.div, lf)

    @boundscheck nlayers == diagn.nlayers || throw(BoundsError)

    fill!(u_mean_grid, 0)           # reset accumulators from previous vertical average
    fill!(v_mean_grid, 0)
    fill!(div_mean_grid, 0)
    fill!(div_mean, 0)

    @inbounds for k in 1:nlayers    # integrate from top to bottom

        # arrays for layer-thickness weighted column averages
        Δσₖ = σ_levels_thick[k]

        # GRID-POINT SPACE: u, v, D with thickness weighting Δσₖ
        # before this k's u, v, D are added to ū, v̄, D̄ store in the
        # sum_above fields for a 1:k-1 integration
        # which is =0 for k=1 as ū, v̄, D̄ accumulators are 0-initialised
        for ij in eachgridpoint(u_mean_grid, v_mean_grid, div_mean_grid)
            # for the Σ_r=1^k-1 Δσᵣ(Dᵣ +  u̲⋅∇lnpₛ) vertical integration
            # Simmons and Burridge, 1981 eq 3.12 split into div and u̲⋅∇lnpₛ
            div_sum_above[ij, k] = div_mean_grid[ij]
            uv∇lnp_sum_above[ij, k] = u_mean_grid[ij]*∇lnp_x[ij] + v_mean_grid[ij]*∇lnp_y[ij]

            u_mean_grid[ij] += u_grid[ij, k]*Δσₖ  # now add the k-th element to the sum
            v_mean_grid[ij] += v_grid[ij, k]*Δσₖ
            div_mean_grid[ij] += div_grid[ij, k]*Δσₖ
        end
        
        # SPECTRAL SPACE: divergence
        for lm in eachharmonic(div, div_mean)
            div_mean[lm] += div[lm, k]*Δσₖ
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

with ū, v̄ being the vertically averaged velocities; px, py the gradients
of the logarithm of surface pressure ln(p_s) and D̄ the vertically averaged divergence.
1. Calculate ∇ln(p_s) in spectral space, convert to grid.
2. Multiply ū, v̄ with ∇ln(p_s) in grid-point space, convert to spectral.
3. D̄ is subtracted in spectral space.
4. Set tendency of the l=m=0 mode to 0 for better mass conservation."""
function surface_pressure_tendency!(
    diagn::DiagnosticVariables,
    S::SpectralTransform,
)
    (; pres_tend, pres_tend_grid) = diagn.tendencies
    (; ∇lnp_x, ∇lnp_y, u_mean_grid, v_mean_grid, div_mean, scratch_memory) = diagn.dynamics
    
    # in grid-point space the the (ū, v̄)⋅∇lnpₛ term (swap sign in spectral)
    @. pres_tend_grid = u_mean_grid*∇lnp_x + v_mean_grid*∇lnp_y
    transform!(pres_tend, pres_tend_grid, scratch_memory, S)
    
    # for semi-implicit div_mean is calc at time step i-1 in vertical_integration!
    @. pres_tend = -pres_tend - div_mean    # add the -div_mean term in spectral, swap sign
    
    pres_tend[1] = 0                # for mass conservation
    return nothing
end

function vertical_velocity!(
    diagn::DiagnosticVariables,
    geometry::Geometry,
)
    (; σ_levels_thick, σ_levels_half, nlayers) = geometry
    (; σ_tend) = diagn.dynamics                  
    
    # sum of Δσ-weighted div, uv∇lnp from 1:k-1
    (; div_sum_above, uv∇lnp, uv∇lnp_sum_above) = diagn.dynamics
    (; div_mean_grid) = diagn.dynamics          # vertical avrgd div to be added to ūv̄∇lnp
    (; σ_tend) = diagn.dynamics                 # vertical mass flux M = pₛσ̇ at k+1/2
    (; div_grid) = diagn.grid
    ūv̄∇lnp = diagn.tendencies.pres_tend_grid    # calc'd in surface_pressure_tendency! (excl -D̄)

    grids_match(σ_tend, div_sum_above, div_grid, uv∇lnp_sum_above, uv∇lnp) ||
        throw(DimensionMismatch(σ_tend, div_sum_above, div_grid, uv∇lnp_sum_above, uv∇lnp))

    @inbounds for k in 1:nlayers-1
        Δσₖ = σ_levels_thick[k]
        σₖ_half = σ_levels_half[k+1]

        for ij in eachgridpoint(σ_tend)
            # Hoskins and Simmons, 1975 just before eq. (6)
            σ_tend[ij, k] = σₖ_half*(div_mean_grid[ij] + ūv̄∇lnp[ij]) -
                (div_sum_above[ij, k] + Δσₖ*div_grid[ij, k]) -
                (uv∇lnp_sum_above[ij, k] + Δσₖ*uv∇lnp[ij, k])
        end
    end

    # mass flux σ̇ is zero at k=1/2 (not explicitly stored) and k=nlayers+1/2 (stored in layer k)
    # set to zero for bottom layer then
    σ_tend[:, nlayers] .= 0
    return nothing
end

"""
$(TYPEDSIGNATURES)
Function barrier to unpack `model`."""
function vordiv_tendencies!(
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)
    (; coriolis, atmosphere, geometry, spectral_transform) = model
    vordiv_tendencies!(diagn, coriolis, atmosphere, geometry, spectral_transform)
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
    diagn::DiagnosticVariables,
    coriolis::AbstractCoriolis,
    atmosphere::AbstractAtmosphere,
    geometry::AbstractGeometry,
    S::SpectralTransform,
)
    (; R_dry) = atmosphere                      # gas constant for dry air
    (; f) = coriolis                            # coriolis parameter
    (; coslat⁻¹) = geometry

     # tendencies already contain parameterizations + advection, therefore accumulate
    (; u_tend_grid, v_tend_grid) = diagn.tendencies
    (; u_grid, v_grid, vor_grid, temp_virt_grid) = diagn.grid   # velocities, vorticity
    (; ∇lnp_x, ∇lnp_y, scratch_memory) = diagn.dynamics         # zonal/meridional gradient of logarithm of surface pressure

    # precompute ring indices and boundscheck
    rings = eachring(u_tend_grid, v_tend_grid, u_grid, v_grid, vor_grid, temp_virt_grid)

    @inbounds for k in eachlayer(u_tend_grid, v_tend_grid)
        for (j, ring) in enumerate(rings)
            coslat⁻¹j = coslat⁻¹[j]
            f_j = f[j]
            for ij in ring
                ω = vor_grid[ij, k] + f_j               # absolute vorticity
                RTᵥ = R_dry*temp_virt_grid[ij, k]       # dry gas constant * virtual temperature anomaly(!)
                u_tend_grid[ij, k] = (u_tend_grid[ij, k] + v_grid[ij, k]*ω - RTᵥ*∇lnp_x[ij])*coslat⁻¹j
                v_tend_grid[ij, k] = (v_tend_grid[ij, k] - u_grid[ij, k]*ω - RTᵥ*∇lnp_y[ij])*coslat⁻¹j
            end
        end
    end

    # divergence and curl of that u, v_tend vector for vor, div tendencies
    (; vor_tend, div_tend ) = diagn.tendencies
    u_tend = diagn.dynamics.a
    v_tend = diagn.dynamics.b

    transform!(u_tend, u_tend_grid, scratch_memory, S)
    transform!(v_tend, v_tend_grid, scratch_memory, S)

    curl!(vor_tend, u_tend, v_tend, S)         # ∂ζ/∂t = ∇×(u_tend, v_tend)
    divergence!(div_tend, u_tend, v_tend, S)   # ∂D/∂t = ∇⋅(u_tend, v_tend)
    return nothing
end

"""$(TYPEDSIGNATURES)
For dynamics=false, after calling parameterization_tendencies! call this function
to transform the physics tendencies from grid-point to spectral space including the
necessary coslat⁻¹ scaling."""
function physics_tendencies_only!(
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)
    (; scratch_memory) = diagn.dynamics
    (; coslat⁻¹) = model.geometry
    S = model.spectral_transform

    # already contain parameterizations
    (; u_tend_grid, v_tend_grid, temp_tend_grid, humid_tend_grid) = diagn.tendencies
    RingGrids._scale_lat!(u_tend_grid, coslat⁻¹)
    RingGrids._scale_lat!(v_tend_grid, coslat⁻¹)

    # divergence and curl of that u, v_tend vector for vor, div tendencies
    (; vor_tend, div_tend, temp_tend, humid_tend) = diagn.tendencies
    u_tend = diagn.dynamics.a
    v_tend = diagn.dynamics.b

    transform!(u_tend, u_tend_grid, scratch_memory, S)
    transform!(v_tend, v_tend_grid, scratch_memory, S)
    transform!(temp_tend, temp_tend_grid, scratch_memory, S)
    model isa PrimitiveWet && transform!(humid_tend, humid_tend_grid, scratch_memory, S)

    curl!(vor_tend, u_tend, v_tend, S)         # ∂ζ/∂t = ∇×(u_tend, v_tend)
    divergence!(div_tend, u_tend, v_tend, S)   # ∂D/∂t = ∇⋅(u_tend, v_tend)
    return nothing
end

# function barrier
function temperature_tendency!(
    diagn::DiagnosticVariables,
    model::PrimitiveEquation,
)
    (; adiabatic_conversion, atmosphere, implicit, geometry, spectral_transform) = model
    temperature_tendency!(diagn, adiabatic_conversion, atmosphere, implicit,
                            geometry, spectral_transform)
end

"""
$(TYPEDSIGNATURES)
Compute the temperature tendency

    ∂T/∂t += -∇⋅((u, v)*T') + T'D + κTᵥ*Dlnp/Dt

`+=` because the tendencies already contain parameterizations and vertical advection.
`T'` is the anomaly with respect to the reference/average temperature. Tᵥ is the virtual
temperature used in the adiabatic term κTᵥ*Dlnp/Dt."""
function temperature_tendency!(
    diagn::DiagnosticVariables,
    adiabatic_conversion::AbstractAdiabaticConversion,
    atmosphere::AbstractAtmosphere,
    implicit::ImplicitPrimitiveEquation,
    G::Geometry,
    S::SpectralTransform,
)
    (; temp_tend, temp_tend_grid)  = diagn.tendencies
    (; div_grid, temp_grid) = diagn.grid
    (; uv∇lnp, uv∇lnp_sum_above, div_sum_above, scratch_memory) = diagn.dynamics
    (; κ) = atmosphere              # thermodynamic kappa = R_Dry/heat_capacity
    Tᵥ = diagn.grid.temp_virt_grid  # anomaly wrt to Tₖ
    (; temp_profile) = implicit

    # semi-implicit: terms here are explicit+implicit evaluated at time step i
    # implicit_correction! then calculated the implicit terms from Vi-1 minus Vi
    # to move the implicit terms to i-1 which is cheaper then the alternative below

    @inbounds for k in eachlayer(temp_tend_grid, temp_grid, div_grid, Tᵥ, div_sum_above, uv∇lnp_sum_above)
        Tₖ = temp_profile[k]    # average layer temperature from reference profile
        
        # coefficients from Simmons and Burridge 1981
        σ_lnp_A = adiabatic_conversion.σ_lnp_A[k]   # eq. 3.12, -1/Δσₖ*ln(σ_k+1/2/σ_k-1/2)
        σ_lnp_B = adiabatic_conversion.σ_lnp_B[k]   # eq. 3.12 -αₖ
        
        # Adiabatic conversion term following Simmons and Burridge 1981 but for σ coordinates 
        # += as tend already contains parameterizations + vertical advection
        for ij in eachgridpoint(temp_tend_grid, temp_grid, div_grid, Tᵥ)
            temp_tend_grid[ij, k] +=
                temp_grid[ij, k] * div_grid[ij, k] +            # +T'D term of hori advection
                κ * (Tᵥ[ij, k] + Tₖ)*(                          # +κTᵥ*Dlnp/Dt, adiabatic term
                σ_lnp_A * (div_sum_above[ij, k] + uv∇lnp_sum_above[ij, k]) +    # eq. 3.12 1st term
                σ_lnp_B * (div_grid[ij, k] + uv∇lnp[ij, k]) +                   # eq. 3.12 2nd term
                uv∇lnp[ij, k])                                                  # eq. 3.13
        end
    end

    transform!(temp_tend, temp_tend_grid, scratch_memory, S)

    # now add the -∇⋅((u, v)*T') term
    flux_divergence!(temp_tend, temp_grid, diagn, G, S, add=true, flipsign=true)
    
    return nothing
end

function humidity_tendency!(diagn::DiagnosticVariables,
                            model::PrimitiveWet)
    G = model.geometry
    S = model.spectral_transform

    (; humid_tend, humid_tend_grid ) = diagn.tendencies
    (; humid_grid ) = diagn.grid

    # add horizontal advection to parameterization + vertical advection tendencies
    horizontal_advection!(humid_tend, humid_tend_grid, humid_grid, diagn, G, S, add=true)
end

# no humidity tendency for dry core
humidity_tendency!(::DiagnosticVariables, ::PrimitiveDry) = nothing

function tracer_advection!(
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    G = model.geometry
    S = model.spectral_transform

    for (name, tracer) in model.tracers
        tracer_tend = diagn.tendencies.tracers_tend[name]
        tracer_tend_grid = diagn.tendencies.tracers_tend_grid[name]
        tracer_grid = diagn.grid.tracers_grid[name]
        
        # add horizontal advection to parameterization + vertical advection + forcing/drag tendencies
        tracer.active && horizontal_advection!(tracer_tend, tracer_tend_grid, tracer_grid, diagn, G, S, add=true)
    end
end

function horizontal_advection!( 
    A_tend::LowerTriangularArray,       # Ouput: tendency to write into
    A_tend_grid::AbstractField,         # Input: tendency incl prev terms
    A_grid::AbstractField,              # Input: grid field to be advected
    diagn::DiagnosticVariables,
    G::Geometry,
    S::SpectralTransform;
    add::Bool=true,                     # add/overwrite A_tend_grid?
)

    (; div_grid) = diagn.grid
    (; scratch_memory) = diagn.dynamics
    
    kernel = add ? (a,b,c) -> a+b*c : (a,b,c) -> b*c

    for k in eachlayer(A_tend_grid, A_grid, div_grid)
        # +A*div term of the advection operator
        @inbounds for ij in eachgridpoint(A_tend_grid, A_grid, div_grid)
            # add as tend already contains parameterizations + vertical advection
            A_tend_grid[ij, k] = kernel(A_tend_grid[ij, k], A_grid[ij, k], div_grid[ij, k])
        end
    end

    transform!(A_tend, A_tend_grid, scratch_memory, S)  # for +A*div in spectral space
    
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
function flux_divergence!(
    A_tend::LowerTriangularArray,   # Ouput: tendency to write into
    A_grid::AbstractField,          # Input: grid field to be advected
    diagn::DiagnosticVariables,     # for u_grid, v_grid
    G::Geometry,
    S::SpectralTransform;
    add::Bool=true,                 # add result to A_tend or overwrite for false
    flipsign::Bool=true,            # compute -∇⋅((u, v)*A) (true) or ∇⋅((u, v)*A)? 
)
    (; u_grid, v_grid) = diagn.grid
    (; scratch_memory) = diagn.dynamics
    (; coslat⁻¹) = G

    # reuse general work arrays a, b, a_grid, b_grid
    uA = diagn.dynamics.a           # = u*A in spectral
    vA = diagn.dynamics.b           # = v*A in spectral
    uA_grid = diagn.dynamics.a_grid # = u*A on grid
    vA_grid = diagn.dynamics.b_grid # = v*A on grid

    # precomputed ring indices and check grids_match
    rings = eachring(A_grid, u_grid, v_grid)
    @inbounds for k in eachlayer(u_grid, v_grid)
        for (j, ring) in enumerate(rings)
            coslat⁻¹j = coslat⁻¹[j]
            for ij in ring
                Acoslat⁻¹j = A_grid[ij, k]*coslat⁻¹j
                uA_grid[ij, k] = u_grid[ij, k]*Acoslat⁻¹j
                vA_grid[ij, k] = v_grid[ij, k]*Acoslat⁻¹j
            end
        end
    end

    transform!(uA, uA_grid, scratch_memory, S)
    transform!(vA, vA_grid, scratch_memory, S)

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
function vorticity_flux_curldiv!(   
                                    diagn::DiagnosticVariables,
                                    coriolis::AbstractCoriolis,
                                    geometry::Geometry,
                                    S::SpectralTransform;
                                    div::Bool=true,     # also calculate div of vor flux?
                                    add::Bool=false)    # accumulate in vor/div tendencies?
    
    (; f) = coriolis
    (; coslat⁻¹) = geometry

    (; u_tend_grid, v_tend_grid) = diagn.tendencies     # already contains forcing
    u = diagn.grid.u_grid                               # velocity
    v = diagn.grid.v_grid                               # velocity
    vor = diagn.grid.vor_grid                           # relative vorticity
    (; whichring) = u.grid                              # precomputed ring indices
    # Launch the kernel for vorticity flux calculation
    arch = S.architecture

    launch!(arch, RingGridWorkOrder, size(u), _vorticity_flux_kernel!,
            u_tend_grid, v_tend_grid, u, v, vor, f, coslat⁻¹, whichring)
    synchronize(arch)
    
    # divergence and curl of that u, v_tend vector for vor, div tendencies
    (; vor_tend, div_tend ) = diagn.tendencies
    u_tend = diagn.dynamics.a
    v_tend = diagn.dynamics.b

    transform!(u_tend, u_tend_grid, S)
    transform!(v_tend, v_tend_grid, S)

    curl!(vor_tend, u_tend, v_tend, S; add)                 # ∂ζ/∂t = ∇×(u_tend, v_tend)
    div && divergence!(div_tend, u_tend, v_tend, S; add)    # ∂D/∂t = ∇⋅(u_tend, v_tend)
    return nothing       
end

@kernel inbounds=true function _vorticity_flux_kernel!(
    u_tend_grid, v_tend_grid, u, v, vor, @Const(f), @Const(coslat⁻¹), @Const(whichring)
)
    # Get indices
    ij, k = @index(Global, NTuple)
    j = whichring[ij]
    
    # Get the coriolis parameter and cosine latitude factor for this latitude
    f_j = f[j]
    coslat⁻¹j = coslat⁻¹[j]
    
    # Calculate absolute vorticity
    ω = vor[ij, k] + f_j
    
    # Update tendencies
    u_tend_grid[ij, k] = (u_tend_grid[ij, k] + v[ij, k] * ω) * coslat⁻¹j
    v_tend_grid[ij, k] = (v_tend_grid[ij, k] - u[ij, k] * ω) * coslat⁻¹j
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
function vorticity_flux!(diagn::DiagnosticVariables, model::ShallowWater)
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
function vorticity_flux!(diagn::DiagnosticVariables, model::Barotropic)
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
function bernoulli_potential!(
    diagn::DiagnosticVariables,   
    S::SpectralTransform,
)   
    (; u_grid, v_grid ) = diagn.grid
    (; geopot, scratch_memory ) = diagn.dynamics
    bernoulli = diagn.dynamics.a                            # reuse work arrays a, a_grid
    bernoulli_grid = diagn.dynamics.a_grid
    (; div_tend ) = diagn.tendencies
 
    half = convert(eltype(bernoulli_grid), 0.5)
    @. bernoulli_grid = half*(u_grid^2 + v_grid^2)          # = ½(u² + v²) on grid
    transform!(bernoulli, bernoulli_grid, scratch_memory, S)                # to spectral space
    bernoulli .+= geopot                                    # add geopotential Φ
    ∇²!(div_tend, bernoulli, S, add=true, flipsign=true)    # add -∇²(½(u² + v²) + ϕ)
    return nothing
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
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Int,                # leapfrog index to evaluate tendencies on
    atmosphere::AbstractAtmosphere,
    implicit::ImplicitPrimitiveEquation,
)                          
    (; R_dry) = atmosphere                  # dry gas constant 
    (; temp_profile) = implicit             # reference profile at layer k
    pres = get_step(progn.pres, lf)         # logarithm of surface pressure at leapfrog index lf
    (; geopot) = diagn.dynamics

    # -R_dry*Tₖ*∇²lnpₛ, linear part of the ∇⋅RTᵥ∇lnpₛ pressure gradient term
    # Tₖ being the reference temperature profile, the anomaly term T' = Tᵥ - Tₖ is calculated
    # vordiv_tendencies! include as R_dry*Tₖ*lnpₛ into the geopotential on which the operator
    # -∇² is applied in bernoulli_potential!
    @inbounds for k in eachmatrix(geopot)
        R_dryTₖ = R_dry*temp_profile[k]
        for lm in eachharmonic(pres)
            geopot[lm, k] += R_dryTₖ*pres[lm]
        end
    end
end

"""
$(TYPEDSIGNATURES)
Computes the (negative) divergence of the volume fluxes `uh, vh` for the continuity equation, -∇⋅(uh, vh)."""
function volume_flux_divergence!(
    diagn::DiagnosticVariables,
    orog::AbstractOrography,
    atmosphere::AbstractAtmosphere,
    G::AbstractGeometry,
    S::SpectralTransform
)                        

    (; pres_grid) = diagn.grid
    (; pres_tend ) = diagn.tendencies
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
Propagate the spectral state of the prognostic variables `progn` to the
diagnostic variables in `diagn` for the barotropic vorticity model.
Updates grid vorticity, spectral stream function and spectral and grid velocities u, v."""
function SpeedyTransforms.transform!(
    diagn::DiagnosticVariables,   
    progn::PrognosticVariables,
    lf::Integer,
    model::Barotropic;
    kwargs...
)    
    (; vor_grid, u_grid, v_grid ) = diagn.grid
    (; scratch_memory) = diagn.dynamics  

    vor = get_step(progn.vor, lf)   # relative vorticity at leapfrog step lf
    U = diagn.dynamics.a            # reuse work arrays for velocities in spectral
    V = diagn.dynamics.b            # reuse work arrays for velocities in spectral
                                    # U = u*coslat, V=v*coslat
    S = model.spectral_transform

    transform!(vor_grid, vor, scratch_memory, S)    # get vorticity on grid from spectral vor
    
    # get spectral U, V from spectral vorticity via stream function Ψ
    # U = u*coslat = -coslat*∂Ψ/∂lat
    # V = v*coslat = ∂Ψ/∂lon, radius omitted in both cases
    UV_from_vor!(U, V, vor, S)

    # transform from U, V in spectral to u, v on grid (U, V = u, v*coslat)
    transform!(u_grid, U, scratch_memory, S, unscale_coslat=true)
    transform!(v_grid, V, scratch_memory, S, unscale_coslat=true)
 
    for (name, tracer) in model.tracers
        tracer_var = get_step(progn.tracers[name], lf)  # tracer at leapfrog step lf
        tracer.active && transform!(diagn.grid.tracers_grid[name], tracer_var, scratch_memory, S)
    end

    # transform random pattern for random process unless NoRandomProcess
    transform!(diagn, progn, lf, model.random_process, S)

    return nothing
end

"""
$(TYPEDSIGNATURES)
Propagate the spectral state of the prognostic variables `progn` to the
diagnostic variables in `diagn` for the shallow water model. Updates grid vorticity,
grid divergence, grid interface displacement (`pres_grid`) and the velocities
u, v."""
function SpeedyTransforms.transform!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Integer,
    model::ShallowWater;
    kwargs...
)
    (; vor_grid, div_grid, pres_grid, u_grid, v_grid ) = diagn.grid
    vor =  get_step(progn.vor, lf)  # relative vorticity at leapfrog step lf
    div =  get_step(progn.div, lf)  # divergence at leapfrog step lf
    pres = get_step(progn.pres, lf) # interface displacement η at leapfrog step lf

    (; scratch_memory) = diagn.dynamics # reuse work arrays for velocities spectral
                                              # U = u*coslat, V=v*coslat
    U = diagn.dynamics.a 
    V = diagn.dynamics.b                                

    S = model.spectral_transform
    
    transform!(vor_grid,  vor, scratch_memory, S)  # get vorticity on grid from spectral vor
    transform!(div_grid,  div, scratch_memory, S)  # get divergence on grid from spectral div
    transform!(pres_grid, pres, scratch_memory, S)  # get η on grid from spectral η
    
    # get spectral U, V from vorticity and divergence via stream function Ψ and vel potential ϕ
    # U = u*coslat = -coslat*∂Ψ/∂lat + ∂ϕ/dlon
    # V = v*coslat =  coslat*∂ϕ/∂lat + ∂Ψ/dlon
    UV_from_vordiv!(U, V, vor, div, S)

    # transform from U, V in spectral to u, v on grid (U, V = u, v*coslat)
    transform!(u_grid, U, scratch_memory, S, unscale_coslat=true)
    transform!(v_grid, V, scratch_memory, S, unscale_coslat=false)

    for (name, tracer) in model.tracers
        tracer_var = get_step(progn.tracers[name], lf)  # tracer at leapfrog step lf
        tracer.active && transform!(diagn.grid.tracers_grid[name], tracer_var, scratch_memory, S)
    end

    # transform random pattern for random process unless NoRandomProcess
    transform!(diagn, progn, lf, model.random_process, S)

    return nothing
end

"""
$(TYPEDSIGNATURES)
Propagate the spectral state of the prognostic variables `progn` to the
diagnostic variables in `diagn` for primitive equation models. Updates grid vorticity,
grid divergence, grid temperature, pressure (`pres_grid`) and the velocities
u, v."""
function SpeedyTransforms.transform!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Integer,
    model::PrimitiveEquation;
    initialize::Bool = false,
)   
    (; vor_grid, div_grid, pres_grid, u_grid, v_grid, temp_grid, humid_grid, 
    temp_grid_prev, humid_grid_prev, u_grid_prev, v_grid_prev) = diagn.grid

    vor   = get_step(progn.vor, lf)     # relative vorticity at leapfrog step lf
    div   = get_step(progn.div, lf)     # divergence at leapfrog step lf
    temp  = get_step(progn.temp, lf)    # temperature at leapfrog step lf
    humid = get_step(progn.humid, lf)   # humidity at leapfrog step lf
    pres  = get_step(progn.pres, lf)    # logarithm of surface pressure at leapfrog step lf

    (; scratch_memory) = diagn.dynamics

    U = diagn.dynamics.a            # reuse work arrays
    V = diagn.dynamics.b            # U = u*coslat, V=v*coslat
    S = model.spectral_transform

    # retain previous time step for vertical advection and some parameterizations
    if initialize == false              # only store prev after initial step
        @. temp_grid_prev = temp_grid   # this is temperature anomaly wrt to implicit reference profile!
        @. humid_grid_prev = humid_grid
        @. u_grid_prev = u_grid
        @. v_grid_prev = v_grid

        for (name, tracer) in model.tracers
            if tracer.active
                diagn.grid.tracers_grid_prev[name] .= diagn.grid.tracers_grid[name]
            end
        end
    end

    transform!(vor_grid,  vor,  scratch_memory, S)  # get vorticity on grid from spectral vor
    transform!(div_grid,  div,  scratch_memory, S)  # get divergence on grid from spectral div
    transform!(temp_grid, temp, scratch_memory, S)  # -- temperature --
    transform!(pres_grid, pres, scratch_memory, S)  # -- pressure --
    
    if model isa PrimitiveWet
        transform!(humid_grid, humid, scratch_memory, S)        
        hole_filling!(humid_grid, model.hole_filling, model)  # remove negative humidity
    end

    # get spectral U, V from vorticity and divergence via stream function Ψ and vel potential ϕ
    # U = u*coslat = -coslat*∂Ψ/∂lat + ∂ϕ/dlon
    # V = v*coslat =  coslat*∂ϕ/∂lat + ∂Ψ/dlon
    UV_from_vordiv!(U, V, vor, div, S)
    
    # transform from U, V in spectral to u, v on grid (U, V = u, v*coslat)
    transform!(u_grid, U, scratch_memory, S, unscale_coslat=true)
    transform!(v_grid, V, scratch_memory, S, unscale_coslat=true)
    
    # include humidity effect into temp for everything stability-related
    temperature_average!(diagn, temp, S)
    virtual_temperature!(diagn, model)      # temp = virt temp for dry core

    if initialize   # at initial step store prev <- current
        # subtract the reference temperature profile Tₖ as temp_grid is too after every time step
        # the Tₖ is added for the physics parameterizations again
        # technically Tₖ is from model.implicit (which is held constant throughout) integration
        # but given it's the initial step here using the instantaneous diagn.temp_average is the same
        @inbounds for k in eachlayer(temp_grid_prev, temp_grid)
            Tₖ = diagn.temp_average[k]
            for ij in eachgridpoint(temp_grid_prev, temp_grid)
                temp_grid_prev[ij, k] = temp_grid[ij, k] - Tₖ
            end
        end

        @. humid_grid_prev = humid_grid
        @. u_grid_prev = u_grid
        @. v_grid_prev = v_grid

        for (name, tracer) in model.tracers
            if tracer.active
                diagn.grid.tracers_grid_prev[name] .= diagn.grid.tracers_grid[name]
            end
        end
    end

    for (name, tracer) in model.tracers
        tracer_var = get_step(progn.tracers[name], lf)  # tracer at leapfrog step lf
        tracer.active && transform!(diagn.grid.tracers_grid[name], tracer_var, scratch_memory, S)
    end

    # transform random pattern for random process unless NoRandomProcess
    transform!(diagn, progn, lf, model.random_process, S)

    return nothing
end

""" 
$(TYPEDSIGNATURES)
Calculates the average temperature of a layer from the l=m=0 harmonic
and stores the result in `diagn.temp_average`"""
function temperature_average!(
    diagn::DiagnosticVariables,
    temp::LowerTriangularArray,
    S::SpectralTransform,
)
    @inbounds for k in eachmatrix(temp)
        # average from l=m=0 harmonic divided by norm of the sphere
        diagn.temp_average[k] = real(temp[1, k])/S.norm_sphere
    end
end 