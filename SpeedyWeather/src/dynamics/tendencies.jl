"""$(TYPEDSIGNATURES)
Calculate tendencies in grid space for the Barotropic model."""
function grid_tendencies!(vars::Variables, model::Barotropic)
    vorticity_flux_grid_tendencies!(vars, model.coriolis, model.geometry)
    return nothing
end

"""$(TYPEDSIGNATURES)
Calculate tendencies in grid space for the ShallowWater model."""
function grid_tendencies!(vars::Variables, model::ShallowWater)
    vorticity_flux_grid_tendencies!(vars, model.coriolis, model.geometry)
    bernoulli_grid_potential!(vars, model)   # kinetic_energy_grid = ½(u²+v²)+Φ
    return nothing
end

"""$(TYPEDSIGNATURES)
Calculate tendencies in grid space for the PrimitiveEquation model."""
function grid_tendencies!(vars::Variables, model::PrimitiveEquation)
    vordiv_grid_tendencies!(vars, model)             # u/v_tend_grid + kinetic_energy_grid (fused: ½(u²+v²))
    temperature_grid_tendency!(vars, model)          # temp_tend_grid + uT_anomaly_grid, vT_anomaly_grid (fused)
    humidity_grid_tendency!(vars, model)             # humid_tend_grid + uq_grid, vq_grid (fused; no-op for dry)
    surface_pressure_grid_tendency!(vars, model)     # pres_tend_grid += (ū,v̄)·∇lnpₛ
    return nothing
end

"""$(TYPEDSIGNATURES)
Reads transformed spectral tendencies and accumulates the final spectral tendencies for the BarotropicModel."""
function spectral_tendencies!(vars::Variables, model::Barotropic)
    vorticity_flux_spectral_tendencies!(vars, model.spectral_transform; div = false, add = true)
    return nothing
end

"""$(TYPEDSIGNATURES)
Reads transformed spectral tendencies and accumulates the final spectral tendencies for the ShallowWaterModel."""
function spectral_tendencies!(vars::Variables, model::ShallowWater)
    vorticity_flux_spectral_tendencies!(vars, model.spectral_transform; div = true, add = true)
    bernoulli_spectral_potential!(vars, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Reads transformed spectral tendencies and accumulates the final spectral tendencies for the PrimitiveEquationModel."""
function spectral_tendencies!(vars::Variables, model::PrimitiveEquation)
    vordiv_spectral_tendencies!(vars, model)
    temperature_spectral_tendency!(vars, model)
    humidity_spectral_tendency!(vars, model)         # no-op for PrimitiveDry
    bernoulli_spectral_potential!(vars, model)
    surface_pressure_spectral_tendency!(vars, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Calculate all tendencies for the BarotropicModel."""
function dynamics_tendencies!(
        vars::Variables,
        lf::Integer,                    # leapfrog index to evaluate tendencies at
        model::Barotropic,
    )
    forcing!(vars, lf, model)           # = (Fᵤ, Fᵥ) forcing for u, v
    drag!(vars, lf, model)              # drag term for u, v

    # compute tendencies in grid space: u, v
    grid_tendencies!(vars, model)

    # batched transform of grid tendencies to spectral space
    transform!(parent(vars.fused.spectral_tendencies),
               parent(vars.fused.grid_tendencies),
               vars.scratch.transform_memory, model.spectral_transform)

    # accumulates into final spectral tendencies: vorticity
    spectral_tendencies!(vars, model)

    tracer_advection!(vars, model)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Calculate all tendencies for the ShallowWaterModel."""
function dynamics_tendencies!(
        vars::Variables,
        lf::Integer,                        # leapfrog index to evaluate tendencies at
        model::ShallowWater,
    )
    (; planet, atmosphere, orography) = model
    (; spectral_transform, geometry) = model

    # for compatibility with other AbstractModels pressure pres = interface displacement η here
    forcing!(vars, lf, model)   # = (Fᵤ, Fᵥ, Fₙ) forcing for u, v, η
    drag!(vars, lf, model)      # drag term for u, v

    geopotential!(vars, planet)             # geopotential Φ = gη in shallow water

    # compute tendencies in grid space: u, v, kinetic_energy
    grid_tendencies!(vars, model)
    
    # batched transform of grid tendencies to spectral space
    transform!(parent(vars.fused.spectral_tendencies),
               parent(vars.fused.grid_tendencies),
               vars.scratch.transform_memory, spectral_transform)

    # accumulates into final spectral tendencies: vorticity, divergence
    spectral_tendencies!(vars, model)

    # = -∇⋅(uh, vh), tendency for "pressure" η
    volume_flux_divergence!(vars, orography, atmosphere, geometry, spectral_transform)

    # advect all tracers
    tracer_advection!(vars, model)

    return nothing
end

"""$(TYPEDSIGNATURES)
Calculate all tendencies for the PrimitiveEquation model (wet or dry)."""
function dynamics_tendencies!(
        vars::Variables,
        lf::Integer,                # leapfrog index for tendencies
        model::PrimitiveEquation,
    )
    forcing!(vars, lf, model)
    drag!(vars, lf, model)

    (; orography, geometry, spectral_transform, geopotential, atmosphere, implicit) = model

    # for semi-implicit corrections (α >= 0.5) linear gravity-wave related tendencies are
    # evaluated at previous timestep i-1 (i.e. lf=1 leapfrog time step)
    # nonlinear terms and parameterizations are always evaluated at lf
    lf_implicit = implicit.α == 0 ? lf : 1

    # calculate ∇ln(pₛ), then (u_k, v_k)⋅∇ln(p_s)
    pressure_gradient_flux!(vars, lf, spectral_transform)

    # calculate Tᵥ = T + Tₖμq in spectral as a approxmation to Tᵥ = T(1+μq) used for geopotential
    linear_virtual_temperature!(vars, lf_implicit, model)

    # temperature relative to profile
    # TODO: broadcast with LTA doesn't work here becasue of a broadcast conflict (temp profile and temp_grid are different dimensions and array types)
    vars.grid.temperature.data .-= implicit.temp_profile'

    # from ∂Φ/∂ln(pₛ) = -RTᵥ for bernoulli_potential!
    geopotential!(vars, geopotential, orography)

    # get ū, v̄, D̄ on grid; D̄ in spectral
    vertical_integration!(vars, lf_implicit, geometry)

    # calculate vertical velocity σ̇ in sigma coordinates for the vertical mass flux M = pₛ * σ̇
    vertical_velocity!(vars, geometry)

    # add the RTₖlnpₛ term to geopotential
    linear_pressure_gradient!(vars, lf_implicit, atmosphere, implicit)

    # use σ̇ for the vertical advection of u, v, T, q
    vertical_advection!(vars, model)


    # compute tendencies in grid space: u, v, temperature, pressure, u·T'·coslat⁻¹, v·T'·coslat⁻¹, kinetic energy, (wet model: humidity, u·q·coslat⁻¹, v·q·coslat⁻¹)
    grid_tendencies!(vars, model)

    
    # batched transform of grid tendencies to spectral space
    transform!(parent(vars.fused.spectral_tendencies),
               parent(vars.fused.grid_tendencies),
               vars.scratch.transform_memory, spectral_transform)

    # accumulates into final spectral tendencies: vorticity, divergence, temperature, divergence, pressure (, humidity)
    spectral_tendencies!(vars, model)

    # advect all tracers
    tracer_advection!(vars, model)

    # back to absolute temperature
    vars.grid.temperature.data .+= implicit.temp_profile'

    return nothing
end

"""$(TYPEDSIGNATURES)
Compute the gradient ∇ln(pₛ) of the logarithm of surface pressure,
followed by its flux, (u,v) * ∇ln(pₛ)."""
function pressure_gradient_flux!(
        vars::Variables,
        lf::Integer,                   # leapfrog index
        S::AbstractSpectralTransform,
    )
    progn = vars.prognostic
    scratch_memory = vars.scratch.transform_memory

    # PRESSURE GRADIENT
    pres = get_step(progn.pressure, lf)             # log of surface pressure at leapfrog step lf
    dpres_dx_spec = vars.dynamics.dpres_dx_spec     # view of slot 1 in :dpres_grad_spec parent
    dpres_dy_spec = vars.dynamics.dpres_dy_spec     # view of slot 2 in :dpres_grad_spec parent
    (; dpres_dx, dpres_dy) = vars.dynamics          # views of slot 1 / 2 in :dpres_grad parent

    ∇!(dpres_dx_spec, dpres_dy_spec, pres, S)       # CALCULATE ∇ln(pₛ)

    # One batched spectral→grid transform for both gradients 
    transform!(parent(vars.fused.dpres_grad),
               parent(vars.fused.dpres_grad_spec),
               scratch_memory, S, unscale_coslat = true)

    (; u, v) = vars.grid
    uv∇lnp = vars.dynamics.pres_flux

    # PRESSURE GRADIENT FLUX
    uv∇lnp .= u .* dpres_dx .+ v .* dpres_dy

    return nothing
end

"""$(TYPEDSIGNATURES)
Calculates the vertically averaged (weighted by the thickness of the σ level)
velocities (`*coslat`) and divergence. E.g.

    u_mean = ∑_k=1^nlayers Δσ_k * u_k

u, v are averaged in grid-point space, divergence in spectral space.
"""
@inline vertical_integration!(
    vars::Variables,
    lf::Integer,                    # leapfrog index for D̄_spec
    geometry::Geometry,
) = vertical_integration!(geometry.spectral_grid.architecture, vars, lf, geometry)

# For the vertical integration and vertical average, the kernel version is unreasonably slow
# on CPU, that's why we have two seperate versions for this function
function vertical_integration!(
        ::CPU,
        vars::Variables,
        lf::Integer,                        # leapfrog index for D̄_spec
        geometry::Geometry,
    )
    (; σ_levels_thick, nlayers) = geometry
    (; dpres_dx, dpres_dy) = vars.dynamics      # zonal, meridional grad of log surface pressure
    (; u, v) = vars.grid
    div_grid = vars.grid.divergence
    (; u_mean_grid, v_mean_grid, div_mean_grid, div_mean) = vars.dynamics
    (; div_sum_above, pres_flux_sum_above) = vars.dynamics
    div = get_step(vars.prognostic.divergence, lf)

    fill!(u_mean_grid, 0)                   # reset accumulators from previous vertical average
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
            pres_flux_sum_above[ij, k] = u_mean_grid[ij] * dpres_dx[ij] + v_mean_grid[ij] * dpres_dy[ij]

            u_mean_grid[ij] += u[ij, k] * Δσₖ  # now add the k-th element to the sum
            v_mean_grid[ij] += v[ij, k] * Δσₖ
            div_mean_grid[ij] += div_grid[ij, k] * Δσₖ
        end

        # SPECTRAL SPACE: divergence
        for lm in eachharmonic(div, div_mean)
            div_mean[lm] += div[lm, k] * Δσₖ
        end
    end
    return nothing
end

function vertical_integration!(
        ::GPU,
        vars::Variables,
        lf::Integer,                    # leapfrog index for D̄_spec
        geometry::Geometry,
    )

    (; σ_levels_thick, nlayers) = geometry
    (; dpres_dx, dpres_dy) = vars.dynamics    # zonal, meridional grad of log surface pressure
    (; u, v) = vars.grid
    div_grid = vars.grid.divergence
    (; u_mean_grid, v_mean_grid, div_mean_grid, div_mean) = vars.dynamics
    (; div_sum_above, pres_flux_sum_above) = vars.dynamics
    div = get_step(vars.prognostic.divergence, lf)

    fill!(u_mean_grid, 0)           # reset accumulators from previous vertical average
    fill!(v_mean_grid, 0)
    fill!(div_mean_grid, 0)
    fill!(div_mean, 0)

    # GRID-POINT SPACE: u, v, D with thickness weighting Δσₖ
    arch = architecture(u_mean_grid)
    launch!(
        arch, RingGridWorkOrder, (size(u_mean_grid, 1),), _vertical_integration_kernel!,
        u_mean_grid, v_mean_grid, div_mean_grid, div_sum_above, pres_flux_sum_above,
        u, v, div_grid, dpres_dx, dpres_dy, σ_levels_thick, nlayers
    )

    # SPECTRAL SPACE: divergence (computed with kernel)
    launch!(
        arch, SpectralWorkOrder, (size(div_mean, 1),), _vertical_integration_spectral_kernel!,
        div_mean, div, σ_levels_thick, nlayers
    )

    return nothing
end

@kernel inbounds = true function _vertical_integration_kernel!(
        u_mean_grid,            # Output: vertically averaged zonal velocity
        v_mean_grid,            # Output: vertically averaged meridional velocity
        div_mean_grid,          # Output: vertically averaged divergence
        div_sum_above,          # Output: sum of div from layers above
        pres_flux_sum_above,    # Output: sum of uv∇lnp from layers above
        u_grid,                 # Input: zonal velocity
        v_grid,                 # Input: meridional velocity
        div_grid,               # Input: divergence
        dpres_dx,                 # Input: zonal gradient of log surface pressure
        dpres_dy,                 # Input: meridional gradient of log surface pressure
        σ_levels_thick,         # Input: layer thicknesses
        nlayers,                # Input: number of layers
    )
    ij = @index(Global, Linear)  # global index: grid point ij

    # Initialize accumulators for this grid point
    u_mean = zero(eltype(u_mean_grid))
    v_mean = zero(eltype(v_mean_grid))
    div_mean = zero(eltype(div_mean_grid))

    # Loop over layers, integrating from top to bottom
    for k in 1:nlayers
        Δσₖ = σ_levels_thick[k]

        # Store sum from layers 1:k-1 before adding k-th layer
        # for the Σ_r=1^k-1 Δσᵣ(Dᵣ + u̲⋅∇lnpₛ) vertical integration
        # Simmons and Burridge, 1981 eq 3.12 split into div and u̲⋅∇lnpₛ
        div_sum_above[ij, k] = div_mean
        pres_flux_sum_above[ij, k] = u_mean * dpres_dx[ij] + v_mean * dpres_dy[ij]

        # Add k-th layer contribution to the running sum
        u_mean += u_grid[ij, k] * Δσₖ
        v_mean += v_grid[ij, k] * Δσₖ
        div_mean += div_grid[ij, k] * Δσₖ
    end

    # Store final accumulated values
    u_mean_grid[ij] = u_mean
    v_mean_grid[ij] = v_mean
    div_mean_grid[ij] = div_mean
end

@kernel inbounds = true function _vertical_integration_spectral_kernel!(
        div_mean,               # Output: vertically averaged divergence (spectral)
        div,                    # Input: divergence (spectral)
        σ_levels_thick, # Input: layer thicknesses
        nlayers,       # Input: number of layers
    )
    lm = @index(Global, Linear)  # global index: harmonic lm

    # Initialize accumulator for this harmonic
    div_sum = zero(eltype(div_mean))

    # Loop over layers, accumulating weighted divergence
    for k in 1:nlayers
        Δσₖ = σ_levels_thick[k]
        div_sum += div[lm, k] * Δσₖ
    end

    # Store final accumulated value
    div_mean[lm] = div_sum
end

"""
$(TYPEDSIGNATURES)

TODO: OLD VERSION / SEQUENTIAL VERSION: MIGHT BE DELETED

Computes the tendency of the logarithm of surface pressure as

    -(ū*px + v̄*py) - D̄

with ū, v̄ being the vertically averaged velocities; px, py the gradients
of the logarithm of surface pressure ln(pₛ) and D̄ the vertically averaged divergence.
1. Calculate ∇ln(pₛ) in spectral space, convert to grid.
2. Multiply ū, v̄ with ∇ln(pₛ) in grid-point space, convert to spectral.
3. D̄ is subtracted in spectral space.
4. Set tendency of the l=m=0 mode to 0 for better mass conservation."""
function surface_pressure_tendency!(
        vars::Variables,
        S::AbstractSpectralTransform,
    )
    surface_pressure_grid_tendency!(vars)

    pres_tend = vars.tendencies.pressure
    pres_tend_grid = vars.tendencies.grid.pressure
    scratch_memory = vars.scratch.transform_memory
    transform!(pres_tend, pres_tend_grid, scratch_memory, S)

    surface_pressure_spectral_tendency!(vars)
    return nothing
end

"""$(TYPEDSIGNATURES)

The tendency of the logarithm of surface pressure is computed as

    -(ū*px + v̄*py) - D̄

with ū, v̄ being the vertically averaged velocities; px, py the gradients
of the logarithm of surface pressure ln(pₛ) and D̄ the vertically averaged divergence.

Here, computes the grid tendency contribuation of the logarithm of surface pressure by:
* ∇ln(pₛ)/px,py is previously computed in grid space in `pressure_gradient_flux!``, 
* Multiply ū, v̄ with ∇ln(pₛ) in grid-point space."""
function surface_pressure_grid_tendency!(vars::Variables)
    pres_tend_grid = vars.tendencies.grid.pressure
    (; dpres_dx, dpres_dy, u_mean_grid, v_mean_grid) = vars.dynamics
    @. pres_tend_grid += u_mean_grid * dpres_dx + v_mean_grid * dpres_dy
    return nothing
end

surface_pressure_grid_tendency!(vars::Variables, ::PrimitiveEquation) =
    surface_pressure_grid_tendency!(vars)

"""$(TYPEDSIGNATURES)

The tendency of the logarithm of surface pressure is computed as

    -(ū*px + v̄*py) - D̄

with ū, v̄ being the vertically averaged velocities; px, py the gradients
of the logarithm of surface pressure ln(pₛ) and D̄ the vertically averaged divergence.

Here, we accumlates the spectral tendency of the logarithm of surface pressure: 
* D̄ is subtracted in spectral space.
* Set tendency of the l=m=0 mode to 0 for better mass conservation."""
function surface_pressure_spectral_tendency!(vars::Variables)
    pres_tend = vars.tendencies.pressure
    div_mean = vars.dynamics.div_mean
    @. pres_tend = -pres_tend - div_mean
    pres_tend.data[1:1] .= 0
    return nothing
end

surface_pressure_spectral_tendency!(vars::Variables, ::PrimitiveEquation) =
    surface_pressure_spectral_tendency!(vars)


"""$(TYPEDSIGNATURES)
Compute vertical velocity."""
function vertical_velocity!(
        vars::Variables,
        geometry::Geometry,
    )
    (; σ_levels_thick, σ_levels_half, nlayers) = geometry

    # sum of Δσ-weighted div, uv∇lnp from 1:k-1
    (; div_sum_above, pres_flux, pres_flux_sum_above) = vars.dynamics
    (; div_mean_grid) = vars.dynamics           # vertical avrgd div to be added to ūv̄∇lnp
    div_grid = vars.grid.divergence

    # vertical velocity in sigma coordinates, positive down
    (; w) = vars.dynamics                       # = vertical mass flux M = pₛσ̇ at k+1/2

    # to calculate u_mean_grid*dpres_dx + v_mean_grid*dpres_dy again
    (; dpres_dx, dpres_dy, u_mean_grid, v_mean_grid) = vars.dynamics
    ūv̄∇lnp = vars.scratch.grid.a_2D             # use scratch memory
    @. ūv̄∇lnp = u_mean_grid * dpres_dx + v_mean_grid * dpres_dy

    grids_match(w, div_sum_above, div_grid, pres_flux_sum_above, pres_flux) ||
        throw(DimensionMismatch(w, div_sum_above, div_grid, pres_flux_sum_above, pres_flux))

    # Hoskins and Simmons, 1975 just before eq. (6)
    Δσₖ = view(σ_levels_thick, 1:(nlayers - 1))'
    σₖ_half = view(σ_levels_half, 2:nlayers)'
    # TODO: broadcast issue here, that's why the .data are neeeded
    w.data[:, 1:(nlayers - 1)] .= σₖ_half .* (div_mean_grid.data .+ ūv̄∇lnp.data) .-
        (div_sum_above.data[:, 1:(nlayers - 1)] .+ Δσₖ .* div_grid.data[:, 1:(nlayers - 1)]) .-
        (pres_flux_sum_above.data[:, 1:(nlayers - 1)] .+ Δσₖ .* pres_flux.data[:, 1:(nlayers - 1)])

    # mass flux σ̇ is zero at k=1/2 (not explicitly stored) and k=nlayers+1/2 (stored in layer k)
    # set to zero for bottom layer then
    w.data[:, nlayers] .= 0
    return nothing
end

"""
$(TYPEDSIGNATURES)
Add the linear contribution of the pressure gradient to the geopotential.
The pressure gradient in the divergence equation takes the form

    -∇⋅(Rd * Tᵥ * ∇lnpₛ) = -∇⋅(Rd * Tᵥ' * ∇lnpₛ) - ∇²(Rd * Tₖ * lnpₛ)

So that the second term inside the Laplace operator can be added to the geopotential.
Rd is the gas constant, Tᵥ the virtual temperature, Tᵥ' its anomaly wrt to the
average or reference temperature Tₖ, and ln(pₛ) is the logarithm of surface pressure."""
function linear_pressure_gradient!(
        vars::Variables,
        lf::Int,                # leapfrog index to evaluate tendencies on
        atmosphere::AbstractAtmosphere,
        implicit::AbstractImplicit,
    )
    (; R_dry) = atmosphere                      # dry gas constant
    Tₖ = implicit.temp_profile                  # reference profile at layer k
    lnpₛ = get_step(vars.prognostic.pressure, lf)   # logarithm of surface pressure at leapfrog index lf
    Φ = vars.dynamics.geopotential

    # -R_dry*Tₖ*∇²lnpₛ, linear part of the ∇⋅RTᵥ∇lnpₛ pressure gradient term
    # Tₖ being the reference temperature profile, the anomaly term T' = Tᵥ - Tₖ is calculated
    # vordiv_tendencies! include as R_dry*Tₖ*lnpₛ into the geopotential on which the operator
    # -∇² is applied in bernoulli_potential!
    # TODO: Broadcast issue with LTA, conflicting broadcast styles
    Φ.data .+= R_dry .* Tₖ' .* lnpₛ.data

    return nothing
end

"""
$(TYPEDSIGNATURES)
Function barrier to unpack `model`."""
function vordiv_tendencies!(
        vars::Variables,
        model::PrimitiveEquation,
    )
    (; coriolis, atmosphere, geometry, implicit, spectral_transform) = model
    return vordiv_tendencies!(vars, coriolis, atmosphere, geometry, implicit, spectral_transform)
end

"""$(TYPEDSIGNATURES)

TODO: OLD VERSION / SEQUENTIAL VERSION: MIGHT BE DELETED

Tendencies for vorticity and divergence. Excluding Bernoulli potential with geopotential
and linear pressure gradient inside the Laplace operator, which are added later in
spectral space. 

    u_tend +=  v*(f+ζ) - RTᵥ'*∇lnpₛ_x
    v_tend += -u*(f+ζ) - RTᵥ'*∇lnpₛ_y

`+=` because the tendencies already contain the parameterizations and vertical advection.
`f` is coriolis, `ζ` relative vorticity, `R` the gas constant `Tᵥ'` the virtual temperature
anomaly, `∇lnpₛ` the gradient of surface pressure and `_x` and `_y` its zonal/meridional
components. The tendencies are then curled/dived to get the tendencies for vorticity/divergence in
spectral space

    ∂ζ/∂t = ∇×(u_tend, v_tend)
    ∂D/∂t = ∇⋅(u_tend, v_tend) + ...

`+ ...` because there's more terms added later for divergence."""
function vordiv_tendencies!(
        vars::Variables,
        coriolis::AbstractCoriolis,
        atmosphere::AbstractAtmosphere,
        geometry::AbstractGeometry,
        implicit::AbstractImplicit,
        S::AbstractSpectralTransform,
    )
    vordiv_grid_tendencies!(vars, coriolis, atmosphere, geometry, implicit)

    u_tend_grid = vars.tendencies.grid.u
    v_tend_grid = vars.tendencies.grid.v
    u_tend = vars.dynamics.u_tendency
    v_tend = vars.dynamics.v_tendency
    scratch_memory = vars.scratch.transform_memory
    transform!(u_tend, u_tend_grid, scratch_memory, S)
    transform!(v_tend, v_tend_grid, scratch_memory, S)

    vordiv_spectral_tendencies!(vars, S)
    return nothing
end

# TODO: Might rename this just to u, v spectral tendencies? Because that's what it does, but 
# then it's also nice to always have the symmetric with *_grid_tendencies! _spectral_tendencies!

"""$(TYPEDSIGNATURES)

Tendencies for vorticity and divergence, here just the gridded u and v tendencies are computed.
Additionally computes the kinetic energy on the fly as well taht is saved for later computations.  

Launches `_vordiv_tendencies_kernel!` to add the vorticity flux and pressure gradient terms to `u_tend_grid, v_tend_grid` (which already contain forcing,
drag, and vertical advection contributions); Excludes Bernoulli potential with geopotential
and linear pressure gradient inside the Laplace operator, which are added later in spectral
space:

    u_tend_grid += v·(f + ζ) - R·Tᵥ'·∂lnpₛ/∂x
    v_tend_grid += -u·(f + ζ) - R·Tᵥ'·∂lnpₛ/∂y

`+=` because the tendencies already contain the parameterizations and vertical advection.
`f` is coriolis, `ζ` relative vorticity, `R` the gas constant, `Tᵥ'` the virtual temperature
anomaly, `∇lnpₛ` the gradient of surface pressure and `_x`, `_y` its zonal/meridional components.

The kinetic energy is computed on fly according to `1/2(u^2 + v^2)` and added to the `kinetic_energy_grid`.
"""
function vordiv_grid_tendencies!(
        vars::Variables,
        coriolis::AbstractCoriolis,
        atmosphere::AbstractAtmosphere,
        geometry::AbstractGeometry,
        implicit::AbstractImplicit,
    )
    (; f) = coriolis
    Tₖ = implicit.temp_profile
    (; coslat⁻¹) = geometry

    # tendencies already contain parameterizations + advection, therefore accumulate
    u_tend_grid = vars.tendencies.grid.u
    v_tend_grid = vars.tendencies.grid.v
    (; u, v) = vars.grid
    vor = vars.grid.vorticity
    temp = vars.grid.temperature
    vars.scratch.grid.a .= 0
    humid = haskey(vars.grid, :humidity) ? vars.grid.humidity : vars.scratch.grid.a
    (; dpres_dx, dpres_dy) = vars.dynamics
    kinetic_energy_grid = vars.dynamics.grid.kinetic_energy

    (; whichring) = u_tend_grid.grid
    arch = architecture(u_tend_grid)
    launch!(
        arch, RingGridWorkOrder, size(u_tend_grid), _vordiv_tendencies_bernoulli_kernel!,
        u_tend_grid, v_tend_grid, kinetic_energy_grid, u, v, vor, temp, humid,
        dpres_dx, dpres_dy, Tₖ, f, coslat⁻¹, whichring, atmosphere,
    )
    return nothing
end

vordiv_grid_tendencies!(vars::Variables, model::PrimitiveEquation) =
    vordiv_grid_tendencies!(vars, model.coriolis, model.atmosphere, model.geometry, model.implicit)

"""$(TYPEDSIGNATURES)

Tendencies for vorticity and divergence. Given the gridded u and v tendencies, the tendencies are  
curled/dived to get the tendencies for vorticity/divergence in spectral space

    ∂ζ/∂t = ∇×(u_tend, v_tend)
    ∂D/∂t = ∇⋅(u_tend, v_tend) + ...

`+ ...` because there's more terms added later for divergence."""
function vordiv_spectral_tendencies!(vars::Variables, S::AbstractSpectralTransform)
    vor_tend = vars.tendencies.vorticity
    div_tend = vars.tendencies.divergence
    u_tend = vars.dynamics.u_tendency
    v_tend = vars.dynamics.v_tendency
    curl!(vor_tend, u_tend, v_tend, S, add = true)
    divergence!(div_tend, u_tend, v_tend, S, add = true)
    return nothing
end

vordiv_spectral_tendencies!(vars::Variables, model::PrimitiveEquation) =
    vordiv_spectral_tendencies!(vars, model.spectral_transform)

@kernel inbounds = true function _vordiv_tendencies_bernoulli_kernel!(
        u_tend_grid,            # Input/Output: zonal wind tendency
        v_tend_grid,            # Input/Output: meridional wind tendency
        kinetic_energy_grid,    # Output: kinetic energy ½(u²+v²)
        u_grid,                 # Input: zonal velocity
        v_grid,                 # Input: meridional velocity
        vor_grid,               # Input: relative vorticity
        temp_grid,              # Input: temperature anomaly
        humid_grid,             # Input: humidity
        dpres_dx,               # Input: zonal gradient of log surface pressure
        dpres_dy,               # Input: meridional gradient of log surface pressure
        Tₖ,                     # Input: reference temperature profile
        f,                      # Input: coriolis parameter
        coslat⁻¹,               # Input: 1/cos(latitude) for scaling
        whichring,              # Input: mapping from grid point to latitude ring
        atmosphere,             # Input: atmosphere for R_dry and μ_virt_temp
    )
    ij, k = @index(Global, NTuple)
    j = whichring[ij]           # latitude ring index for this grid point
    coslat⁻¹j = coslat⁻¹[j]     # get coslat⁻¹ for this latitude
    f_j = f[j]                  # coriolis parameter for this latitude

    u_ij = u_grid[ij, k]        # shared load
    v_ij = v_grid[ij, k]        # shared load
    ω = vor_grid[ij, k] + f_j   # absolute vorticity

    # compute virtual temperature on the fly, temp_grid is anomaly
    (; R_dry) = atmosphere
    Tᵥ = virtual_temperature(temp_grid[ij, k] + Tₖ[k], humid_grid[ij, k], atmosphere)
    RTᵥ = R_dry * (Tᵥ - Tₖ[k])    # dry gas constant * virtual temperature anomaly
    u_tend_grid[ij, k] = (u_tend_grid[ij, k] + v_ij * ω - RTᵥ * dpres_dx[ij]) * coslat⁻¹j
    v_tend_grid[ij, k] = (v_tend_grid[ij, k] - u_ij * ω - RTᵥ * dpres_dy[ij]) * coslat⁻¹j

    # kinetic energy ½(u²+v²); KE-only here, geopotential is added in spectral space for PrimEq
    half = convert(eltype(kinetic_energy_grid), 0.5)
    kinetic_energy_grid[ij, k] = half * (u_ij * u_ij + v_ij * v_ij)
end

"""$(TYPEDSIGNATURES)
For `dynamics=false`, after calling `parameterization_tendencies!` call this function
to transform the physics tendencies from grid-point to spectral space including the
necessary coslat⁻¹ scaling."""
function parameterization_tendencies_only!(
        vars::Variables,
        model::PrimitiveEquation,
    )
    scratch_memory = vars.scratch.transform_memory
    (; coslat⁻¹) = model.geometry
    S = model.spectral_transform

    # already contain parameterizations
    u_tend_grid = vars.tendencies.grid.u
    v_tend_grid = vars.tendencies.grid.v
    temp_tend_grid = vars.tendencies.grid.temperature
    RingGrids._scale_lat!(u_tend_grid, coslat⁻¹)
    RingGrids._scale_lat!(v_tend_grid, coslat⁻¹)

    # divergence and curl of that u, v_tend vector for vor, div tendencies
    vor_tend = vars.tendencies.vorticity
    div_tend = vars.tendencies.divergence
    temp_tend = vars.tendencies.temperature
    u_tend = vars.dynamics.u_tendency        # spectral u-tendency 
    v_tend = vars.dynamics.v_tendency        # spectral v-tendency

    transform!(u_tend, u_tend_grid, scratch_memory, S)
    transform!(v_tend, v_tend_grid, scratch_memory, S)
    transform!(temp_tend, temp_tend_grid, scratch_memory, S)

    # humidity only for models that have humidity
    if haskey(vars.tendencies, :humidity)
        humid_tend = vars.tendencies.humidity
        humid_tend_grid = vars.tendencies.grid.humidity
        transform!(humid_tend, humid_tend_grid, scratch_memory, S)
    end

    curl!(vor_tend, u_tend, v_tend, S)         # ∂ζ/∂t = ∇×(u_tend, v_tend)
    divergence!(div_tend, u_tend, v_tend, S)   # ∂D/∂t = ∇⋅(u_tend, v_tend)
    return nothing
end

# function barrier
function temperature_tendency!(
        vars::Variables,
        model::PrimitiveEquation,
    )
    (; adiabatic_conversion, atmosphere, implicit, geometry, spectral_transform) = model
    return temperature_tendency!(
        vars, adiabatic_conversion, atmosphere, implicit,
        geometry, spectral_transform
    )
end

"""$(TYPEDSIGNATURES)

TODO: OLD VERSION / SEQUENTIAL VERSION: MIGHT BE DELETED

Compute the temperature tendency.

    ∂T/∂t += -∇⋅((u, v)*T') + T'D + κTᵥ*Dlnp/Dt

`+=` because the tendencies already contain parameterizations and vertical advection.
`T'` is the anomaly with respect to the reference/average temperature. Tᵥ is the virtual
temperature used in the adiabatic term κTᵥ*Dlnp/Dt."""
function temperature_tendency!(
        vars::Variables,
        adiabatic_conversion::AbstractAdiabaticConversion,
        atmosphere::AbstractAtmosphere,
        implicit::AbstractImplicit,
        G::Geometry,
        S::AbstractSpectralTransform,
    )
    temperature_grid_tendency!(vars, adiabatic_conversion, atmosphere, implicit, G)

    temp_tend = vars.tendencies.temperature
    temp_tend_grid = vars.tendencies.grid.temperature
    uT_grid = vars.dynamics.grid.uT_anomaly
    vT_grid = vars.dynamics.grid.vT_anomaly
    uT_spec = vars.dynamics.uT_anomaly
    vT_spec = vars.dynamics.vT_anomaly
    scratch_memory = vars.scratch.transform_memory
    transform!(temp_tend, temp_tend_grid, scratch_memory, S)
    transform!(uT_spec, uT_grid, scratch_memory, S)
    transform!(vT_spec, vT_grid, scratch_memory, S)

    temperature_spectral_tendency!(vars, S)
    return nothing
end

"""$(TYPEDSIGNATURES)

Compute the temperature tendency.

    ∂T/∂t += -∇⋅((u, v)*T') + T'D + κTᵥ*Dlnp/Dt

`+=` because the tendencies already contain parameterizations and vertical advection.
`T'` is the anomaly with respect to the reference/average temperature. Tᵥ is the virtual
temperature used in the adiabatic term κTᵥ*Dlnp/Dt.

Here, the gridded tendency contributions are computed in a single fused kernel that:
* adds the adiabatic + T'D terms to the temperature tendency
* writes `(uT_anomaly_grid, vT_anomaly_grid) = coslat⁻¹·(u·T', v·T')` for the flux divergence."""
function temperature_grid_tendency!(
        vars::Variables,
        adiabatic_conversion::AbstractAdiabaticConversion,
        atmosphere::AbstractAtmosphere,
        implicit::AbstractImplicit,
        G::Geometry,
    )
    temp_tend_grid = vars.tendencies.grid.temperature
    div_grid = vars.grid.divergence
    temp = vars.grid.temperature
    (; u, v) = vars.grid

    vars.scratch.grid.a .= 0
    humid = haskey(vars.grid, :humidity) ? vars.grid.humidity : vars.scratch.grid.a

    (; pres_flux, pres_flux_sum_above, div_sum_above) = vars.dynamics
    uT_anomaly_grid = vars.dynamics.grid.uT_anomaly
    vT_anomaly_grid = vars.dynamics.grid.vT_anomaly
    (; temp_profile) = implicit
    (; coslat⁻¹) = G
    (; whichring) = temp_tend_grid.grid

    # semi-implicit: terms here are explicit+implicit evaluated at time step i
    # implicit_correction! then calculated the implicit terms from Vi-1 minus Vi
    # to move the implicit terms to i-1 which is cheaper then the alternative below

    # Fused: adiabatic + T'D term on temp_tend_grid, plus (uT, vT) = coslat⁻¹·(u·T', v·T')
    arch = architecture(temp_tend_grid)
    launch!(
        arch, RingGridWorkOrder, size(temp_tend_grid), _temperature_tendency_flux_kernel!,
        temp_tend_grid, uT_anomaly_grid, vT_anomaly_grid,
        temp, u, v, div_grid, humid, div_sum_above, pres_flux_sum_above, pres_flux,
        temp_profile, adiabatic_conversion.σ_lnp_A, adiabatic_conversion.σ_lnp_B,
        coslat⁻¹, whichring, atmosphere
    )
    return nothing
end

temperature_grid_tendency!(vars::Variables, model::PrimitiveEquation) =
    temperature_grid_tendency!(vars, model.adiabatic_conversion, model.atmosphere, model.implicit, model.geometry)

"""$(TYPEDSIGNATURES)

Compute the temperature tendency.

    ∂T/∂t += -∇⋅((u, v)*T') + T'D + κTᵥ*Dlnp/Dt

`+=` because the tendencies already contain parameterizations and vertical advection.
`T'` is the anomaly with respect to the reference/average temperature. Tᵥ is the virtual
temperature used in the adiabatic term κTᵥ*Dlnp/Dt.

Here, the previosuly computed gridded tendency contribution are accumulated and the spectral
tendency computed with divergence is computed via the `flux_spectral_divergence!` function."""
function temperature_spectral_tendency!(vars::Variables, S::AbstractSpectralTransform)
    temp_tend = vars.tendencies.temperature
    uT_spec = vars.dynamics.uT_anomaly
    vT_spec = vars.dynamics.vT_anomaly
    flux_spectral_divergence!(temp_tend, uT_spec, vT_spec, S; add = true, flipsign = true)
    return nothing
end

temperature_spectral_tendency!(vars::Variables, model::PrimitiveEquation) =
    temperature_spectral_tendency!(vars, model.spectral_transform)

@kernel inbounds = true function _temperature_tendency_flux_kernel!(
        temp_tend_grid,             # Input/Output: temperature tendency
        uT_anomaly_grid,            # Output: coslat⁻¹·u·T' for flux divergence
        vT_anomaly_grid,            # Output: coslat⁻¹·v·T' for flux divergence
        temp_grid,                  # Input: temperature anomaly
        u_grid,                     # Input: zonal velocity
        v_grid,                     # Input: meridional velocity
        div_grid,                   # Input: divergence
        humid_grid,                 # Input: humidity
        div_sum_above,              # Input: sum of div from layers above
        pres_flux_sum_above,        # Input: sum of pres_flux from layers above
        pres_flux,                  # Input: (u,v)⋅∇lnp term
        temp_profile,               # Input: reference temperature profile
        σ_lnp_A,                    # Input: adiabatic conversion coefficient A
        σ_lnp_B,                    # Input: adiabatic conversion coefficient B
        coslat⁻¹,                   # Input: 1/cos(latitude) for scaling
        whichring,                  # Input: mapping from grid point to latitude ring
        atmosphere,                 # Input: atmosphere for κ and μ_virt_temp
    )

    ij, k = @index(Global, NTuple)
    j = whichring[ij]
    coslat⁻¹j = coslat⁻¹[j]
    Tₖ = temp_profile[k]    # average layer temperature from reference profile

    # coefficients from Simmons and Burridge 1981
    σ_lnp_A_k = σ_lnp_A[k]   # eq. 3.12, -1/Δσₖ*ln(σ_k+1/2/σ_k-1/2)
    σ_lnp_B_k = σ_lnp_B[k]   # eq. 3.12 -αₖ

    T = temp_grid[ij, k]    # shared load: anomaly T'

    # Adiabatic conversion term following Simmons and Burridge 1981 but for σ coordinates
    # += as tend already contains parameterizations + vertical advection
    Tᵥ = virtual_temperature(T + Tₖ, humid_grid[ij, k], atmosphere)
    (; κ) = atmosphere
    temp_tend_grid[ij, k] +=
        T * div_grid[ij, k] +                               # +T'D term of hori advection
        κ * Tᵥ * (                                          # +κTᵥ*Dlnp/Dt, adiabatic term
        σ_lnp_A_k * (div_sum_above[ij, k] + pres_flux_sum_above[ij, k]) +  # eq. 3.12 1st term
            σ_lnp_B_k * (div_grid[ij, k] + pres_flux[ij, k]) +             # eq. 3.12 2nd term
            pres_flux[ij, k]
    )                                                        # eq. 3.13

    # (uT_anomaly, vT_anomaly) = coslat⁻¹·(u·T', v·T') — shares T load with adiabatic block
    Tcoslat⁻¹j = T * coslat⁻¹j
    uT_anomaly_grid[ij, k] = u_grid[ij, k] * Tcoslat⁻¹j
    vT_anomaly_grid[ij, k] = v_grid[ij, k] * Tcoslat⁻¹j
end

#TODO: OLD VERSION / SEQUENTIAL VERSION: MIGHT BE DELETED
function humidity_tendency!(
        vars::Variables,
        model::PrimitiveWet
    )
    G = model.geometry
    S = model.spectral_transform

    humidity_grid_tendency!(vars, model)

    humid_tend = vars.tendencies.humidity
    humid_tend_grid = vars.tendencies.grid.humidity
    uq_grid = vars.dynamics.grid.uq
    vq_grid = vars.dynamics.grid.vq
    uq_spec = vars.dynamics.uq
    vq_spec = vars.dynamics.vq
    scratch_memory = vars.scratch.transform_memory
    transform!(humid_tend, humid_tend_grid, scratch_memory, S)
    transform!(uq_spec, uq_grid, scratch_memory, S)
    transform!(vq_spec, vq_grid, scratch_memory, S)

    humidity_spectral_tendency!(vars, model)
    return nothing
end

# no humidity tendency for dry core
humidity_tendency!(::Variables, ::PrimitiveDry) = nothing

"""$(TYPEDSIGNATURES)

Computes the gridded contribution to the humidity tendency `humid_tend_grid`. Adds the
`+q·div` advection term to `humid_tend_grid` and writes the `(uq, vq)` flux intermediates
to the grid-side named slots — fused into a single kernel that shares the `humid_grid` load."""
function humidity_grid_tendency!(vars::Variables, model::PrimitiveWet)
    G = model.geometry
    humid_tend_grid = vars.tendencies.grid.humidity
    humid_grid = vars.grid.humidity
    div_grid = vars.grid.divergence
    (; u, v) = vars.grid
    uq_grid = vars.dynamics.grid.uq
    vq_grid = vars.dynamics.grid.vq
    (; coslat⁻¹) = G
    (; whichring) = humid_tend_grid.grid

    arch = architecture(humid_tend_grid)
    launch!(
        arch, RingGridWorkOrder, size(humid_tend_grid), _humidity_advection_flux_kernel!,
        humid_tend_grid, uq_grid, vq_grid,
        humid_grid, u, v, div_grid, coslat⁻¹, whichring
    )
    return nothing
end
humidity_grid_tendency!(::Variables, ::PrimitiveDry) = nothing

@kernel inbounds = true function _humidity_advection_flux_kernel!(
        humid_tend_grid,            # Input/Output: humidity tendency, accumulates +q·D
        uq_grid,                    # Output: coslat⁻¹·u·q
        vq_grid,                    # Output: coslat⁻¹·v·q
        humid_grid,                 # Input: humidity
        u_grid,                     # Input: zonal velocity
        v_grid,                     # Input: meridional velocity
        div_grid,                   # Input: divergence
        coslat⁻¹,                   # Input: 1/cos(latitude) for scaling
        whichring,                  # Input: mapping from grid point to latitude ring
    )
    ij, k = @index(Global, NTuple)
    j = whichring[ij]
    coslat⁻¹j = coslat⁻¹[j]

    q = humid_grid[ij, k]   # shared load

    # += q·D term of the advection operator (tend already contains parameterizations + vertical advection)
    humid_tend_grid[ij, k] = muladd(q, div_grid[ij, k], humid_tend_grid[ij, k])

    # (uq, vq) = coslat⁻¹·(u·q, v·q)
    qcoslat⁻¹j = q * coslat⁻¹j
    uq_grid[ij, k] = u_grid[ij, k] * qcoslat⁻¹j
    vq_grid[ij, k] = v_grid[ij, k] * qcoslat⁻¹j
end

"""$(TYPEDSIGNATURES)

Computes the spectral humidity tendency via the `horizontal_spectral_advection!`
Adds `-∇⋅(uq, vq)` to the previously computed gridded and transformed tendency."""
function humidity_spectral_tendency!(vars::Variables, model::PrimitiveWet)
    S = model.spectral_transform
    humid_tend = vars.tendencies.humidity
    uq_spec = vars.dynamics.uq
    vq_spec = vars.dynamics.vq
    horizontal_spectral_advection!(humid_tend, uq_spec, vq_spec, S)
    return nothing
end
humidity_spectral_tendency!(::Variables, ::PrimitiveDry) = nothing

function tracer_advection!(
        vars::Variables,
        model::AbstractModel,
    )
    G = model.geometry
    S = model.spectral_transform

    for (name, tracer) in model.tracers
        name_grid = Symbol(name, "_grid")
        tracer_tend = vars.tendencies.tracers[name]
        tracer_tend_grid = vars.tendencies.tracers[name_grid]
        tracer_grid = vars.grid.tracers[name]

        # add horizontal advection to parameterization + vertical advection + forcing/drag tendencies
        tracer.active && horizontal_advection!(tracer_tend, tracer_tend_grid, tracer_grid, vars, G, S, add = true)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)

#TODO: OLD VERSION / SEQUENTIAL VERSION: MIGHT BE DELETED

Compute the horizontal advection (combined grid + transform + spec reduction).
"""
function horizontal_advection!(
        A_tend::LowerTriangularArray,       # Output: tendency to write into
        A_tend_grid::AbstractField,         # Input: tendency incl prev terms
        A_grid::AbstractField,              # Input: grid field to be advected
        vars::Variables,
        G::Geometry,
        S::AbstractSpectralTransform;
        add::Bool = true,                   # add/overwrite A_tend_grid?
        # Forwarded to `flux_divergence!`; defaults to the unfused scratch.{a,b} for unfused callers.
        uA = vars.scratch.a,
        vA = vars.scratch.b,
        uA_grid = vars.scratch.grid.a,
        vA_grid = vars.scratch.grid.b,
    )
    horizontal_grid_advection!(A_tend_grid, A_grid, vars, G; add, uA_grid, vA_grid)

    scratch_memory = vars.scratch.transform_memory
    transform!(A_tend, A_tend_grid, scratch_memory, S)  # for +A*div in spectral space
    transform!(uA, uA_grid, scratch_memory, S)
    transform!(vA, vA_grid, scratch_memory, S)

    horizontal_spectral_advection!(A_tend, uA, vA, S)
    return nothing
end

"""$(TYPEDSIGNATURES)

Compute the gridded contribution to the horizontal advection. 

Writes `+A*div` to `A_tend_grid` and `(u*A, v*A)` to `(uA_grid, vA_grid)`
"""
function horizontal_grid_advection!(
        A_tend_grid::AbstractField,         # Input/Output: A_tend on grid, accumulates +A*div term
        A_grid::AbstractField,              # Input: grid field to be advected
        vars::Variables,
        G::Geometry;
        add::Bool = true,                   # use muladd (true) or overwrite (false) for the +A*div term
        uA_grid = vars.dynamics.grid.uT_anomaly,   # caller picks the correct named slot
        vA_grid = vars.dynamics.grid.vT_anomaly,
    )

    # barotropic model doesn't have divergence, the +A*div term is then zero
    if haskey(vars.grid, :divergence)
        div_grid = vars.grid.divergence

        kernel_func = add ? muladd : @inline (a, b, c) -> a * b

        # Launch kernel to compute +A*div term of the advection operator
        arch = architecture(A_tend_grid)
        launch!(
            arch, RingGridWorkOrder, size(A_tend_grid), _horizontal_advection_kernel!,
            A_tend_grid, A_grid, div_grid, kernel_func
        )
    end

    # write u*A and v*A on grid
    flux_grid_divergence!(uA_grid, vA_grid, A_grid, vars, G)
    return nothing
end

"""$(TYPEDSIGNATURES)

Compute the spectral tendencies due to horizontal advection. 

Computes `A_tend += -∇⋅(uA, vA)`. `A_tend` is assumed to already
hold the spectral form of `A_tend_grid` (= forcing + parameterizations 
+ the `+A*div` term written by `horizontal_grid_advection!`).
"""
function horizontal_spectral_advection!(
        A_tend::LowerTriangularArray,
        uA::LowerTriangularArray,
        vA::LowerTriangularArray,
        S::AbstractSpectralTransform,
    )
    # A_tend += -∇⋅(uA, vA)
    flux_spectral_divergence!(A_tend, uA, vA, S; add = true, flipsign = true)
    return nothing
end

@kernel inbounds = true function _horizontal_advection_kernel!(
        A_tend_grid,    # Input/Output: tendency grid
        A_grid,         # Input: field to be advected
        div_grid,       # Input: divergence field
        kernel_func,    # Input: kernel function (muladd or multiply)
    )
    ij, k = @index(Global, NTuple)  # global index: grid point ij, layer k
    # +A*div term of the advection operator
    # add as tend already contains parameterizations + vertical advection
    A_tend_grid[ij, k] = kernel_func(A_grid[ij, k], div_grid[ij, k], A_tend_grid[ij, k])
end

"""$(TYPEDSIGNATURES)
Computes `∇⋅((u, v)*A)` with the option to add/overwrite `A_tend` and to
`flip_sign` of the flux divergence by doing so.

- `A_tend =  ∇⋅((u, v)*A)` for `add=false`, `flip_sign=false`
- `A_tend = -∇⋅((u, v)*A)` for `add=false`, `flip_sign=true`
- `A_tend += ∇⋅((u, v)*A)` for `add=true`, `flip_sign=false`
- `A_tend -= ∇⋅((u, v)*A)` for `add=true`, `flip_sign=true`
"""
function flux_divergence!(
        A_tend::LowerTriangularArray,   # Output: tendency to write into
        A_grid::AbstractField,          # Input: grid field to be advected
        vars::Variables,                # for u,v on grid and scratch memory
        G::Geometry,
        S::AbstractSpectralTransform;
        add::Bool = true,               # add result to A_tend or overwrite for false
        flipsign::Bool = true,          # compute -∇⋅((u, v)*A) (true) or ∇⋅((u, v)*A)?
        # Named slots for the spec/grid intermediates u*A and v*A. Default to the unfused
        # :a/:b scratches (used by tracer_advection!, volume_flux_divergence! for η);
        # named callers (temperature_tendency!, humidity_tendency!) pass slots from the
        # :spectral_tendencies/:grid_tendencies fuse parents (uT_anomaly/vT_anomaly, uq/vq).
        # When the slots are fused, dycore callers prefer the grid/spectral split below.
        uA = vars.scratch.a,            # = u*A in spectral
        vA = vars.scratch.b,            # = v*A in spectral
        uA_grid = vars.scratch.grid.a,  # = u*A on grid
        vA_grid = vars.scratch.grid.b,  # = v*A on grid
    )
    flux_grid_divergence!(uA_grid, vA_grid, A_grid, vars, G)

    scratch_memory = vars.scratch.transform_memory
    transform!(uA, uA_grid, scratch_memory, S)
    transform!(vA, vA_grid, scratch_memory, S)

    flux_spectral_divergence!(A_tend, uA, vA, S; add, flipsign)
    return nothing
end

"""$(TYPEDSIGNATURES)
Gridded half of `flux_divergence!`: writes `u*A` and `v*A` to `(uA_grid, vA_grid)` on the grid,
with the standard coslat⁻¹ scaling."""
function flux_grid_divergence!(
        uA_grid::AbstractField,         # Output: u*A on grid (named slot in :grid_tendencies)
        vA_grid::AbstractField,         # Output: v*A on grid
        A_grid::AbstractField,          # Input: grid field to be advected
        vars::Variables,
        G::Geometry,
    )
    (; u, v) = vars.grid
    (; coslat⁻¹) = G
    (; whichring) = A_grid.grid
    arch = architecture(A_grid)
    launch!(
        arch, RingGridWorkOrder, size(A_grid), _flux_divergence_kernel!,
        uA_grid, vA_grid, A_grid, u, v, coslat⁻¹, whichring
    )
    return nothing
end

#TODO: just call `divergence!` directly instead or keep this to have the *_spectral_ *_grid_ symmetry? 
"""$(TYPEDSIGNATURES)
Purely spectral half of `flux_divergence!`: Just computes the actual divergence"""
function flux_spectral_divergence!(
        A_tend::LowerTriangularArray,
        uA::LowerTriangularArray,
        vA::LowerTriangularArray,
        S::AbstractSpectralTransform;
        add::Bool = true,
        flipsign::Bool = true,
    )
    divergence!(A_tend, uA, vA, S; add, flipsign)
    return nothing
end

@kernel inbounds = true function _flux_divergence_kernel!(
        uA_grid,                # Output: u*A on grid
        vA_grid,                # Output: v*A on grid
        A_grid,                 # Input: field to be advected
        u_grid,                 # Input: zonal velocity
        v_grid,                 # Input: meridional velocity
        coslat⁻¹,       # Input: 1/cos(latitude) for scaling
        whichring,      # Input: mapping from grid point to latitude ring
    )
    I = @index(Global, Cartesian)

    j = whichring[I[1]]               # latitude ring index for this grid point
    coslat⁻¹j = coslat⁻¹[j]         # get coslat⁻¹ for this latitude
    Acoslat⁻¹j = A_grid[I] * coslat⁻¹j
    uA_grid[I] = u_grid[I] * Acoslat⁻¹j
    vA_grid[I] = v_grid[I] * Acoslat⁻¹j
end

"""
$(TYPEDSIGNATURES)

#TODO: OLD VERSION / SEQUENTIAL VERSION: MIGHT BE DELETED

Compute the vorticity advection as the curl/div of the vorticity fluxes

    ∂ζ/∂t = ∇×(u_tend, v_tend)
    ∂D/∂t = ∇⋅(u_tend, v_tend)

with

    u_tend = Fᵤ + v*(ζ+f)
    v_tend = Fᵥ - u*(ζ+f)

with `Fᵤ, Fᵥ` from `u_tend_grid`/`v_tend_grid` that are assumed to be alread
set in `forcing!`. Set `div=false` for the BarotropicModel which doesn't
require the divergence tendency."""
function vorticity_flux_curldiv!(
        vars::Variables,
        coriolis::AbstractCoriolis,
        geometry::Geometry,
        S::AbstractSpectralTransform;
        div::Bool = true,     # also calculate div of vor flux?
        add::Bool = false
    )    # accumulate in vor/div tendencies?
    vorticity_flux_grid_tendencies!(vars, coriolis, geometry)

    u_tend_grid = vars.tendencies.grid.u
    v_tend_grid = vars.tendencies.grid.v
    u_tend = vars.dynamics.u_tendency
    v_tend = vars.dynamics.v_tendency
    scratch_memory = vars.scratch.transform_memory
    transform!(u_tend, u_tend_grid, scratch_memory, S)
    transform!(v_tend, v_tend_grid, scratch_memory, S)

    vorticity_flux_spectral_tendencies!(vars, S; div, add)
    return nothing
end

"""$(TYPEDSIGNATURES)

Compute the vorticity advection as the curl/div of the vorticity fluxes

    ∂ζ/∂t = ∇×(u_tend, v_tend)
    ∂D/∂t = ∇⋅(u_tend, v_tend)

with

    u_tend = Fᵤ + v*(ζ+f)
    v_tend = Fᵥ - u*(ζ+f)

with `Fᵤ, Fᵥ` from `u_tend_grid`/`v_tend_grid` that are assumed to be alread
set in `forcing!`. Set `div=false` for the BarotropicModel which doesn't
require the divergence tendency.

Here, we only compute the gridded contriubation `u_tend_grid``, `v_tend_grid`` on top of 
the forcing already accumulated therein."""
function vorticity_flux_grid_tendencies!(
        vars::Variables,
        coriolis::AbstractCoriolis,
        geometry::Geometry,
    )
    (; f) = coriolis
    (; coslat⁻¹) = geometry
    u_tend_grid = vars.tendencies.grid.u          # already contains forcing
    v_tend_grid = vars.tendencies.grid.v          # already contains forcing
    (; u, v) = vars.grid
    vor = vars.grid.vorticity
    (; whichring) = u.grid

    arch = architecture(u_tend_grid)
    launch!(
        arch, RingGridWorkOrder, size(u), _vorticity_flux_kernel!,
        u_tend_grid, v_tend_grid, u, v, vor, f, coslat⁻¹, whichring
    )
    return nothing
end

"""$(TYPEDSIGNATURES)

Compute the vorticity advection as the curl/div of the vorticity fluxes

    ∂ζ/∂t = ∇×(u_tend, v_tend)
    ∂D/∂t = ∇⋅(u_tend, v_tend)

with

    u_tend = Fᵤ + v*(ζ+f)
    v_tend = Fᵥ - u*(ζ+f)

with `Fᵤ, Fᵥ` from `u_tend_grid`/`v_tend_grid` that are assumed to be alread
set in `forcing!`. Set `div=false` for the BarotropicModel which doesn't
require the divergence tendency.

Here, only final spectral contributions `∂ζ/∂t`, `∂D/∂tare computed given previously 
`u_tend`, `v_tend`."""
function vorticity_flux_spectral_tendencies!(
        vars::Variables,
        S::AbstractSpectralTransform;
        div::Bool = true,
        add::Bool = false,
    )
    vor_tend = vars.tendencies.vorticity
    u_tend = vars.dynamics.u_tendency
    v_tend = vars.dynamics.v_tendency

    curl!(vor_tend, u_tend, v_tend, S; add)                   # ∂ζ/∂t = ∇×(u_tend, v_tend)

    if div                                                    # not needed/available in Barotropic
        div_tend = vars.tendencies.divergence
        divergence!(div_tend, u_tend, v_tend, S; add)         # ∂D/∂t = ∇⋅(u_tend, v_tend)
    end
    return nothing
end

@kernel inbounds = true function _vorticity_flux_kernel!(
        u_tend_grid, v_tend_grid, u, v, vor, f, coslat⁻¹, whichring
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

    ∂ζ/∂t = ∇×(u_tend, v_tend)
    ∂D/∂t = ∇⋅(u_tend, v_tend)

with

    u_tend = Fᵤ + v*(ζ+f)
    v_tend = Fᵥ - u*(ζ+f)

with Fᵤ, Fᵥ the forcing from `forcing!` already in `u_tend_grid`/`v_tend_grid` and
vorticity ζ, coriolis f."""
vorticity_flux!(vars::Variables, model::ShallowWater) =
    vorticity_flux_curldiv!(vars, model.coriolis, model.geometry, model.spectral_transform, div = true, add = true)

"""
$(TYPEDSIGNATURES)
Vorticity flux tendency in the barotropic vorticity equation

    ∂ζ/∂t = ∇×(u_tend, v_tend)

with

    u_tend = Fᵤ + v*(ζ+f)
    v_tend = Fᵥ - u*(ζ+f)

with Fᵤ, Fᵥ the forcing from `forcing!` already in `u_tend_grid`/`v_tend_grid` and
vorticity ζ, coriolis f."""
vorticity_flux!(vars::Variables, model::Barotropic) =
    vorticity_flux_curldiv!(vars, model.coriolis, model.geometry, model.spectral_transform, div = false, add = true)

#TODO: OLD VERSION / SEQUENTIAL VERSION: MIGHT BE DELETED

function bernoulli_potential!(vars::Variables, model::ShallowWater)
    S = model.spectral_transform
    bernoulli_grid_potential!(vars, model)

    bernoulli = vars.dynamics.kinetic_energy
    bernoulli_grid = vars.dynamics.grid.kinetic_energy
    scratch_memory = vars.scratch.transform_memory
    transform!(bernoulli, bernoulli_grid, scratch_memory, S)

    bernoulli_spectral_potential!(vars, model)
    return nothing
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
        vars::Variables,
        S::AbstractSpectralTransform,
    )
    bernoulli_grid_potential!(vars, S)

    bernoulli = vars.dynamics.kinetic_energy
    bernoulli_grid = vars.dynamics.grid.kinetic_energy
    scratch_memory = vars.scratch.transform_memory

    # TODO
    # Tₖ*lnpₛ on grid, use broadcasting as T is 3D but surface pressure is 2D
    # Add the linear contribution of the pressure gradient to the geopotential.
    # The pressure gradient in the divergence equation takes the form
    #     -∇⋅(Rd * Tᵥ * ∇lnpₛ) = -∇⋅(Rd * Tᵥ' * ∇lnpₛ) - ∇²(Rd * Tₖ * lnpₛ)
    # So that the second term inside the Laplace operator can be added to the geopotential.
    # Rd is the gas constant, Tᵥ the virtual temperature and Tᵥ' its anomaly wrt to the
    # average or reference temperature Tₖ, pₛ is the surface pressure.
    # broadcast 1D Tₖ (1 value per layer) over 2D pₛ (one value per grid point) to 3D
    # Tₖ = model.implicit.temp_profile                # average layer temperature (1D)
    # pₛ = diagn.grid.pres_grid_prev                  # 2D not prev is in Pa
    # RdTlnpₛ .= R_dry * Tₖ' .* log.(pₛ)

    transform!(bernoulli, bernoulli_grid, scratch_memory, S)

    bernoulli_spectral_potential!(vars, S)
    return nothing
end

"""$(TYPEDSIGNATURES)
Gridded contribution of `bernoulli_potential!` for ShallowWater: writes `kinetic_energy_grid = ½(u²+v²)+Φ`
(geopotential Φ is included because in ShallowWater geopotential lives only on the grid)."""
function bernoulli_grid_potential!(vars::Variables, ::ShallowWater)
    (; u, v) = vars.grid
    Φ = vars.grid.geopotential
    bernoulli_grid = vars.dynamics.grid.kinetic_energy
    half = convert(eltype(bernoulli_grid), 0.5)
    @. bernoulli_grid = half * (u^2 + v^2) + Φ
    return nothing
end

"""$(TYPEDSIGNATURES)
Gridded contribution of `bernoulli_potential!` for PrimitiveEquation: writes `kinetic_energy_grid = ½(u²+v²)`."""
function bernoulli_grid_potential!(vars::Variables, ::Union{PrimitiveEquation, AbstractSpectralTransform})
    (; u, v) = vars.grid
    bernoulli_grid = vars.dynamics.grid.kinetic_energy
    half = convert(eltype(bernoulli_grid), 0.5)
    @. bernoulli_grid = half * (u^2 + v^2)
    return nothing
end

"""$(TYPEDSIGNATURES)
Spectral half of `bernoulli_potential!`. `kinetic_energy` is assumed to already hold the spec
transform of `kinetic_energy_grid` (from the mega-batched transform). Adds `-∇²(KE)` into
`div_tend`. For SW, geopotential is already absorbed into `kinetic_energy_grid` on the grid
side; for PrimitiveEquation, the spectral geopotential is added here first."""
function bernoulli_spectral_potential!(vars::Variables, model::ShallowWater)
    _bernoulli_spectral_potential!(vars, model.spectral_transform)
    return nothing
end

function bernoulli_spectral_potential!(vars::Variables, model::PrimitiveEquation)
    bernoulli_spectral_potential!(vars, model.spectral_transform)
    return nothing
end

function bernoulli_spectral_potential!(vars::Variables, S::AbstractSpectralTransform)
    # PrimitiveEquation path: add spectral geopotential to KE before the Laplacian.
    bernoulli = vars.dynamics.kinetic_energy
    geopot = vars.dynamics.geopotential
    bernoulli .+= geopot
    _bernoulli_spectral_potential!(vars, S)
    return nothing
end

function _bernoulli_spectral_potential!(vars::Variables, S::AbstractSpectralTransform)
    bernoulli = vars.dynamics.kinetic_energy
    div_tend = vars.tendencies.divergence
    ∇²!(div_tend, bernoulli, S, add = true, flipsign = true)
    return nothing
end



"""$(TYPEDSIGNATURES)
Computes the (negative) divergence of the volume fluxes `uh, vh` for the continuity equation, -∇⋅(uh, vh)."""
function volume_flux_divergence!(
        vars::Variables,
        orog::AbstractOrography,
        atmosphere::AbstractAtmosphere,
        G::AbstractGeometry,
        S::AbstractSpectralTransform
    )

    (; η) = vars.grid
    η_tend = vars.tendencies.η
    (; orography) = orog
    H = atmosphere.layer_thickness

    # compute dynamic layer thickness h on the grid
    # η is the interface displacement, update to layer thickness h = η + H - Hb
    # H is the layer thickness at rest without mountains, Hb the orography
    η .+= H .- orography

    # now do -∇⋅(uh, vh) and store in η_tend
    flux_divergence!(η_tend, η, vars, G, S, add = true, flipsign = true)
    return nothing
end


"""
$(TYPEDSIGNATURES)
Calculates the average temperature of a layer from the l=m=0 harmonic
and stores the result in `diagn.temp_average`"""
function temperature_average!(
        vars::Variables,
        temp::LowerTriangularArray,
        S::AbstractSpectralTransform,
    )
    # average from l=m=0 harmonic divided by norm of the sphere
    @. vars.grid.temp_average = real(temp[1, :]) / S.norm_sphere
    return nothing
end

function reset_tendencies!(vars::Variables; value = 0)
    _reset_tendencies!(vars.tendencies, value)
    return vars
end

# recursively fill all arrays in a NamedTuple, unpacking nested NamedTuples
# this avoids Union-typed iteration which Enzyme cannot differentiate
@inline _reset_tendencies!(nt::NamedTuple, value) = _reset_tendencies_inner!(values(nt), value)
@inline _reset_tendencies_inner!(::Tuple{}, _) = nothing
@inline function _reset_tendencies_inner!(t::Tuple, value)
    _reset_tendency!(first(t), value)
    _reset_tendencies_inner!(Base.tail(t), value)
    return nothing
end

# dispatch on element type: nested NamedTuple vs array
@inline _reset_tendency!(nt::NamedTuple, value) = _reset_tendencies_inner!(values(nt), value)
@inline _reset_tendency!(a::AbstractArray, value) = fill!(a, value)
