"""$(TYPEDSIGNATURES)
Calculate tendencies in grid space for the Barotropic model."""
function grid_tendencies!(vars::Variables, model::Barotropic)
    vorticity_flux_grid_tendencies!(vars, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Calculate tendencies in grid space for the ShallowWater model."""
function grid_tendencies!(vars::Variables, model::ShallowWater)
    vorticity_flux_grid_tendencies!(vars, model)
    bernoulli_grid_potential!(vars, model, model.time_stepping)   # kinetic_energy_grid = ½(u²+v²)+Φ
    return nothing
end

"""$(TYPEDSIGNATURES)
Calculate tendencies in grid space for the PrimitiveEquation model."""
function grid_tendencies!(vars::Variables, model::PrimitiveEquation)
    vordiv_grid_tendencies!(vars, model)             # u_tend_grid, v_tend_grid
    temperature_grid_tendency!(vars, model)          # temp_tend_grid + uT_anomaly_grid, vT_anomaly_grid
    humidity_grid_tendency!(vars, model)             # humid_tend_grid + uq_grid, vq_grid (no-op for dry)
    bernoulli_grid_potential!(vars, model, model.time_stepping)           # kinetic_energy_grid = ½(u²+v²)
    surface_pressure_grid_tendency!(vars, model)     # pres_tend_grid += (ū,v̄)·∇lnpₛ
    return nothing
end

"""$(TYPEDSIGNATURES)
Reads transformed spectral tendencies and accumulates the final spectral tendencies for the BarotropicModel."""
function spectral_tendencies!(vars::Variables, model::Barotropic)
    vorticity_flux_spectral_tendencies!(vars, model.spectral_transform, model.time_stepping; div = false, add = true)
    return nothing
end

"""$(TYPEDSIGNATURES)
Reads transformed spectral tendencies and accumulates the final spectral tendencies for the ShallowWaterModel."""
function spectral_tendencies!(vars::Variables, model::ShallowWater)
    vorticity_flux_spectral_tendencies!(vars, model.spectral_transform, model.time_stepping; div = true, add = true)
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
        model::Barotropic,
    )
    (; time_stepping) = model

    forcing!(vars, model)               # = (Fᵤ, Fᵥ) forcing for u, v
    drag!(vars, model)                  # drag term for u, v
    scale_tendencies!(vars, model)      # dynamical core uses a scaled time step, Δt/radius

    # compute tendencies in grid space: u, v
    grid_tendencies!(vars, model)

    # batched transform of grid tendencies to spectral space
    transform!(get_tendency_step(parent(vars.fused.spectral_tendencies), time_stepping, DynamicalCore()),
               get_tendency_step(parent(vars.fused.grid_tendencies), time_stepping, DynamicalCore()),
               vars.scratch.transform_memory, model.spectral_transform)

    # accumulates into final spectral tendencies: vorticity
    spectral_tendencies!(vars, model)

    tracer_advection!(vars, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Calculate all tendencies for the ShallowWaterModel."""
function dynamics_tendencies!(
        vars::Variables,
        model::ShallowWater,
    )
    (; spectral_transform, time_stepping) = model

    forcing!(vars, model)               # = (Fᵤ, Fᵥ, Fₙ) forcing for u, v, η
    drag!(vars, model)                  # drag term for u, v
    scale_tendencies!(vars, model)      # dynamical core uses a scaled time step, Δt/radius

    geopotential!(vars, model)          # geopotential Φ = gη in shallow water

    # compute tendencies in grid space: u, v, kinetic_energy
    grid_tendencies!(vars, model)
    
    # batched transform of grid tendencies to spectral space
    transform!(get_tendency_step(parent(vars.fused.spectral_tendencies), time_stepping, DynamicalCore()),
               get_tendency_step(parent(vars.fused.grid_tendencies), time_stepping, DynamicalCore()),
               vars.scratch.transform_memory, spectral_transform)

    # accumulates into final spectral tendencies: vorticity, divergence
    spectral_tendencies!(vars, model)

    # = -∇⋅(uh, vh), tendency for interface displacement η
    volume_flux_divergence!(vars, model)

    # advect all tracers
    tracer_advection!(vars, model)

    return nothing
end

"""$(TYPEDSIGNATURES)
Calculate all tendencies for the PrimitiveEquation model (wet or dry)."""
function dynamics_tendencies!(
        vars::Variables,
        model::PrimitiveEquation,
    )
    forcing!(vars, model)
    drag!(vars, model)
    scale_tendencies!(vars, model)

    (; orography, geometry, spectral_transform, geopotential, atmosphere, implicit, time_stepping) = model

    # calculate ∇ln(pₛ), then (u_k, v_k)⋅∇ln(p_s)
    pressure_gradient_flux!(vars, spectral_transform, time_stepping)

    # calculate Tᵥ = T + Tₖμq in spectral as a approxmation to Tᵥ = T(1+μq) used for geopotential
    linear_virtual_temperature!(vars, model)

    # temperature relative to profile
    # TODO: broadcast with LTA doesn't work here becasue of a broadcast conflict (temp profile and temp_grid are different dimensions and array types)
    T = get_prognostic_step(vars.grid.temperature, time_stepping, DynamicalCore())
    T.data .-= implicit.temp_profile'

    # from ∂Φ/∂ln(pₛ) = -RTᵥ for bernoulli_potential!
    geopotential!(vars, geopotential, orography)

    # get ū, v̄, D̄ on grid; D̄ in spectral
    vertical_integration!(vars, geometry, time_stepping)

    # calculate vertical velocity σ̇ in sigma coordinates for the vertical mass flux M = pₛ * σ̇
    vertical_velocity!(vars, geometry, time_stepping)

    # add the RTₖlnpₛ term to geopotential
    linear_pressure_gradient!(vars, atmosphere, implicit, time_stepping)

    # use σ̇ for the vertical advection of u, v, T, q
    vertical_advection!(vars, model)


    # compute tendencies in grid space: u, v, temperature, pressure, u·T'·coslat⁻¹, v·T'·coslat⁻¹, kinetic energy, (wet model: humidity, u·q·coslat⁻¹, v·q·coslat⁻¹)
    grid_tendencies!(vars, model)

    
    # batched transform of grid tendencies to spectral space
    transform!(get_tendency_step(parent(vars.fused.spectral_tendencies), time_stepping, DynamicalCore()),
               get_tendency_step(parent(vars.fused.grid_tendencies), time_stepping, DynamicalCore()),
               vars.scratch.transform_memory, spectral_transform)

    # accumulates into final spectral tendencies: vorticity, divergence, temperature, divergence, pressure (, humidity)
    spectral_tendencies!(vars, model)

    # advect all tracers
    tracer_advection!(vars, model)

    # back to absolute temperature
    T.data .+= implicit.temp_profile'

    return nothing
end

"""$(TYPEDSIGNATURES)
Compute the gradient ∇ln(pₛ) of the logarithm of surface pressure,
followed by its flux, (u,v) * ∇ln(pₛ)."""
function pressure_gradient_flux!(
        vars::Variables,
        S::AbstractSpectralTransform,
        time_stepping::AbstractTimeStepper,
    )
    progn = vars.prognostic
    scratch_memory = vars.scratch.transform_memory

    # PRESSURE GRADIENT
    pres = get_prognostic_step(progn.pressure, time_stepping, DynamicalCore())
    dpres_dx_spec = vars.dynamics.dpres_dx_spec     # view of slot 1 in :dpres_grad_spec parent
    dpres_dy_spec = vars.dynamics.dpres_dy_spec     # view of slot 2 in :dpres_grad_spec parent
    (; dpres_dx, dpres_dy) = vars.dynamics          # views of slot 1 / 2 in :dpres_grad parent

    ∇!(dpres_dx_spec, dpres_dy_spec, pres, S)       # CALCULATE ∇ln(pₛ)

    # One batched spectral→grid transform for both gradients 
    transform!(parent(vars.fused.dpres_grad),
               parent(vars.fused.dpres_grad_spec),
               scratch_memory, S, unscale_coslat = true)

    u = get_prognostic_step(vars.grid.u, time_stepping, DynamicalCore())
    v = get_prognostic_step(vars.grid.v, time_stepping, DynamicalCore())
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
    geometry::Geometry,
    time_stepping::AbstractTimeStepper,
) = vertical_integration!(geometry.spectral_grid.architecture, vars, geometry, time_stepping)

# For the vertical integration and vertical average, the kernel version is unreasonably slow
# on CPU, that's why we have two seperate versions for this function
function vertical_integration!(
        ::CPU,
        vars::Variables,
        geometry::Geometry,
        time_stepping::AbstractTimeStepper,
    )
    # TODO: σ_levels_thick is used here as the sigma-coordinate layer thickness Δσₖ for the
    # Δσ-weighted vertical integrals (Simmons & Burridge 1981, eq. 3.12). Generalising to
    # hybrid coordinates requires replacing Δσₖ with the actual pressure thickness divided
    # by surface pressure, which varies in space and time.
    (; σ_levels_thick, nlayers) = geometry
    (; dpres_dx, dpres_dy) = vars.dynamics      # zonal, meridional grad of log surface pressure
    u = get_prognostic_step(vars.grid.u, time_stepping, DynamicalCore())
    v = get_prognostic_step(vars.grid.v, time_stepping, DynamicalCore())
    div_grid = get_prognostic_step(vars.grid.divergence, time_stepping, DynamicalCore())
    (; u_mean_grid, v_mean_grid, div_mean_grid, div_mean) = vars.dynamics
    (; div_sum_above, pres_flux_sum_above) = vars.dynamics
    div = get_prognostic_step(vars.prognostic.divergence, time_stepping, LinearDynamicalCore())

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
        geometry::Geometry,
        time_stepping::AbstractTimeStepper,
    )

    # TODO: same as CPU version — Δσₖ weights baked into sigma-coordinate continuity equation.
    (; σ_levels_thick, nlayers) = geometry
    (; dpres_dx, dpres_dy) = vars.dynamics    # zonal, meridional grad of log surface pressure
    u = get_prognostic_step(vars.grid.u, time_stepping, DynamicalCore())
    v = get_prognostic_step(vars.grid.v, time_stepping, DynamicalCore())
    div_grid = get_prognostic_step(vars.grid.divergence, time_stepping, DynamicalCore())
    (; u_mean_grid, v_mean_grid, div_mean_grid, div_mean) = vars.dynamics
    (; div_sum_above, pres_flux_sum_above) = vars.dynamics
    div = get_prognostic_step(vars.prognostic.divergence, time_stepping, LinearDynamicalCore())

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
        σ_levels_thick,         # Input: layer thicknesses
        nlayers,                # Input: number of layers
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

"""$(TYPEDSIGNATURES)

The tendency of the logarithm of surface pressure is computed as

    -(ū*px + v̄*py) - D̄

with ū, v̄ being the vertically averaged velocities; px, py the gradients
of the logarithm of surface pressure ln(pₛ) and D̄ the vertically averaged divergence.

Here, computes the grid tendency contribuation of the logarithm of surface pressure by:
* ∇ln(pₛ)/px,py is previously computed in grid space in `pressure_gradient_flux!``, 
* Multiply ū, v̄ with ∇ln(pₛ) in grid-point space."""
function surface_pressure_grid_tendency!(vars::Variables, time_stepping::AbstractTimeStepper)
    pres_tend_grid = get_tendency_step(vars.tendencies.grid.pressure, time_stepping, DynamicalCore())
    (; dpres_dx, dpres_dy, u_mean_grid, v_mean_grid) = vars.dynamics
    @. pres_tend_grid += u_mean_grid * dpres_dx + v_mean_grid * dpres_dy
    return nothing
end

surface_pressure_grid_tendency!(vars::Variables, model::PrimitiveEquation) =
    surface_pressure_grid_tendency!(vars, model.time_stepping)

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
        time_stepping::AbstractTimeStepper,
    )
    # TODO: σ_levels_thick and σ_levels_half are used here to compute the sigma-coordinate
    # vertical velocity σ̇ (Hoskins & Simmons 1975, eq. before 6). Generalising to hybrid
    # coordinates requires a reformulation of the vertical velocity equation.
    (; σ_levels_thick, σ_levels_half, nlayers) = geometry

    # sum of Δσ-weighted div, uv∇lnp from 1:k-1
    (; div_sum_above, pres_flux, pres_flux_sum_above) = vars.dynamics
    (; div_mean_grid) = vars.dynamics           # vertical avrgd div to be added to ūv̄∇lnp
    div_grid = get_prognostic_step(vars.grid.divergence, time_stepping, DynamicalCore())

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
    # @views so the RHS `[:, 1:nlayers-1]` slices are views, not materialized copies,
    # letting the whole dotted expression fuse into one allocation-free broadcast
    @views w.data[:, 1:(nlayers - 1)] .= σₖ_half .* (div_mean_grid.data .+ ūv̄∇lnp.data) .-
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
        atmosphere::AbstractAtmosphere,
        implicit::AbstractImplicit,
        time_stepping::AbstractTimeStepper,
    )
    (; R_dry) = atmosphere                      # dry gas constant
    Tₖ = implicit.temp_profile                  # reference profile at layer k

    # for Leapfrog this term is evaluated at the previous time step and the
    # implicit corrections will move it to the current as done for all linear gravity-wave related terms
    lnpₛ = get_prognostic_step(vars.prognostic.pressure, time_stepping, LinearDynamicalCore())
    Φ = vars.dynamics.spectral_geopotential

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
    (; coriolis, atmosphere, geometry, implicit, spectral_transform, time_stepping) = model
    return vordiv_tendencies!(vars, coriolis, atmosphere, geometry, implicit, spectral_transform, time_stepping)
end

# TODO: Might rename this just to u, v spectral tendencies? Because that's what it does, but 
# then it's also nice to always have the symmetric with *_grid_tendencies! _spectral_tendencies!
"""$(TYPEDSIGNATURES)

Tendencies for vorticity and divergence, here just the gridded u and v tendencies are computed. 

Launches `_vordiv_tendencies_kernel!` to add the vorticity flux and pressure gradient terms to `u_tend_grid, v_tend_grid` (which already contain forcing,
drag, and vertical advection contributions); Excludes Bernoulli potential with geopotential
and linear pressure gradient inside the Laplace operator, which are added later in spectral
space:

    u_tend_grid += v·(f + ζ) - R·Tᵥ'·∂lnpₛ/∂x
    v_tend_grid += -u·(f + ζ) - R·Tᵥ'·∂lnpₛ/∂y

`+=` because the tendencies already contain the parameterizations and vertical advection.
`f` is coriolis, `ζ` relative vorticity, `R` the gas constant, `Tᵥ'` the virtual temperature
anomaly, `∇lnpₛ` the gradient of surface pressure and `_x`, `_y` its zonal/meridional components."""
function vordiv_grid_tendencies!(
        vars::Variables,
        coriolis::AbstractCoriolis,
        atmosphere::AbstractAtmosphere,
        geometry::AbstractGeometry,
        implicit::AbstractImplicit,
        time_stepping::AbstractTimeStepper,
    )
    (; f) = coriolis
    Tₖ = implicit.temp_profile
    scale = vars.prognostic.scale[]             # scale Coriolis on the fly as is vorticity
    (; coslat⁻¹) = geometry

    # tendencies already contain parameterizations + advection, therefore accumulate
    u_tend_grid = get_tendency_step(vars.tendencies.grid.u, time_stepping, DynamicalCore())
    v_tend_grid = get_tendency_step(vars.tendencies.grid.v, time_stepping, DynamicalCore())
    u = get_prognostic_step(vars.grid.u, time_stepping, DynamicalCore())
    v = get_prognostic_step(vars.grid.v, time_stepping, DynamicalCore())
    vor = get_prognostic_step(vars.grid.vorticity, time_stepping, DynamicalCore())
    temp = get_prognostic_step(vars.grid.temperature, time_stepping, DynamicalCore())

    # if no humidity, pass a zero scratch array to the kernel
    humid = haskey(vars.grid, :humidity) ?
        get_prognostic_step(vars.grid.humidity, time_stepping, DynamicalCore()) :
        fill!(vars.scratch.grid.a, 0)

    (; dpres_dx, dpres_dy) = vars.dynamics              # zonal/meridional gradient of logarithm of surface pressure

    (; whichring) = u_tend_grid.grid
    arch = architecture(u_tend_grid)
    launch!(
        arch, RingGridWorkOrder, size(u_tend_grid), _vordiv_tendencies_kernel!,
        u_tend_grid, v_tend_grid, u, v, vor, temp, humid,
        dpres_dx, dpres_dy, Tₖ, f, scale, coslat⁻¹, whichring, atmosphere,
    )
    return nothing
end

vordiv_grid_tendencies!(vars::Variables, model::PrimitiveEquation) =
    vordiv_grid_tendencies!(vars, model.coriolis, model.atmosphere, model.geometry, model.implicit, model.time_stepping)

"""$(TYPEDSIGNATURES)

Tendencies for vorticity and divergence. Given the gridded u and v tendencies, the tendencies are  
curled/dived to get the tendencies for vorticity/divergence in spectral space

    ∂ζ/∂t = ∇×(u_tend, v_tend)
    ∂D/∂t = ∇⋅(u_tend, v_tend) + ...

`+ ...` because there's more terms added later for divergence."""
function vordiv_spectral_tendencies!(vars::Variables, time_stepping::AbstractTimeStepper, S::AbstractSpectralTransform)
    vor_tend = get_tendency_step(vars.tendencies.vorticity, time_stepping, DynamicalCore())
    div_tend = get_tendency_step(vars.tendencies.divergence, time_stepping, DynamicalCore()) 
    u_tend = get_tendency_step(vars.dynamics.u_tendency, time_stepping, DynamicalCore())
    v_tend = get_tendency_step(vars.dynamics.v_tendency, time_stepping, DynamicalCore())
    curl!(vor_tend, u_tend, v_tend, S, add = true)
    divergence!(div_tend, u_tend, v_tend, S, add = true)
    return nothing
end

vordiv_spectral_tendencies!(vars::Variables, model::PrimitiveEquation) =
    vordiv_spectral_tendencies!(vars, model.time_stepping, model.spectral_transform)

@kernel inbounds = true function _vordiv_tendencies_kernel!(
        u_tend_grid,            # Input/Output: zonal wind tendency
        v_tend_grid,            # Input/Output: meridional wind tendency
        u_grid,                 # Input: zonal velocity
        v_grid,                 # Input: meridional velocity
        vor_grid,               # Input: relative vorticity
        temp_grid,              # Input: temperature anomaly
        humid_grid,             # Input: humidity
        dpres_dx,               # Input: zonal gradient of log surface pressure
        dpres_dy,               # Input: meridional gradient of log surface pressure
        Tₖ,                     # Input: reference temperature profile
        f,                      # Input: coriolis parameter
        scale,                  # Input: Scaling of vorticity, apply to Coriolis on the fly too
        coslat⁻¹,               # Input: 1/cos(latitude) for scaling
        whichring,              # Input: mapping from grid point to latitude ring
        atmosphere,             # Input: atmosphere for R_dry and μ_virt_temp
    )
    ij, k = @index(Global, NTuple)
    j = whichring[ij]           # latitude ring index for this grid point
    coslat⁻¹j = coslat⁻¹[j]     # get coslat⁻¹ for this latitude
    f_j = f[j] * scale          # coriolis parameter for this latitude, scaled
    ω = vor_grid[ij, k] + f_j   # absolute vorticity

    # compute virtual temperature on the fly, temp_grid is anomaly
    (; R_dry) = atmosphere
    Tᵥ = virtual_temperature(temp_grid[ij, k] + Tₖ[k], humid_grid[ij, k], atmosphere)
    RTᵥ = R_dry * (Tᵥ - Tₖ[k])    # dry gas constant * virtual temperature anomaly
    u_tend_grid[ij, k] = (u_tend_grid[ij, k] + v_grid[ij, k] * ω - RTᵥ * dpres_dx[ij]) * coslat⁻¹j
    v_tend_grid[ij, k] = (v_tend_grid[ij, k] - u_grid[ij, k] * ω - RTᵥ * dpres_dy[ij]) * coslat⁻¹j
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
    TS = model.time_stepping

    # physics has filled the grid tendencies (u, v, temperature, humidity). The wind tendencies
    # need the 1/coslat scaling that the curl/divergence operators expect (as in grid_tendencies!).
    u_tend_grid = get_tendency_step(vars.tendencies.grid.u, TS, DummyParameterization())
    v_tend_grid = get_tendency_step(vars.tendencies.grid.v, TS, DummyParameterization())
    RingGrids._scale_lat!(u_tend_grid, coslat⁻¹)
    RingGrids._scale_lat!(v_tend_grid, coslat⁻¹)

    # One batched transform of the full, contiguous grid-tendency parent into the spectral-tendency
    # parent (mirrors dynamics_tendencies!). The fuse slots are aligned, so this maps
    # u → u_tendency, v → v_tendency, temperature → temperature, humidity → humidity. The
    # dynamics-only product slots (uT_anomaly, kinetic_energy, uq, ...) are never written on the
    # physics-only path (grid_tendencies! is not called), so they stay zero and their spectral
    # images are harmless. A per-variable transform of the individual slot views would pass
    # non-contiguous SubArrays to the GPU Legendre transform, whose `reinterpret` fails to compile
    # (InvalidIRError); the full contiguous parent transform is the GPU-safe path.
    transform!(get_tendency_step(parent(vars.fused.spectral_tendencies), TS, DynamicalCore()),
               get_tendency_step(parent(vars.fused.grid_tendencies), TS, DynamicalCore()),
               scratch_memory, S)

    # divergence and curl of the (now spectral) u, v tendencies for vor, div tendencies
    vor_tend = get_tendency_step(vars.tendencies.vorticity, TS, DynamicalCore())
    div_tend = get_tendency_step(vars.tendencies.divergence, TS, DynamicalCore())
    u_tend = get_tendency_step(vars.dynamics.u_tendency, TS, DynamicalCore())   # spectral u-tendency (fused slot, filled above)
    v_tend = get_tendency_step(vars.dynamics.v_tendency, TS, DynamicalCore())   # spectral v-tendency (fused slot, filled above)
    curl!(vor_tend, u_tend, v_tend, S)         # ∂ζ/∂t = ∇×(u_tend, v_tend)
    divergence!(div_tend, u_tend, v_tend, S)   # ∂D/∂t = ∇⋅(u_tend, v_tend)
    return nothing
end

"""$(TYPEDSIGNATURES)

Compute the gridded contribution to the temperature tendency:
* adds the adiabatic + T'D terms to `temp_tend_grid` via `_temperature_tendency_kernel!`
* writes `(uT_anomaly_grid, vT_anomaly_grid) = (u·T', v·T')` for the flux divergence."""
function temperature_grid_tendency!(vars::Variables, model::PrimitiveEquation)
    (; adiabatic_conversion, atmosphere, implicit, time_stepping) = model
    temp_tend_grid = get_tendency_step(vars.tendencies.grid.temperature, time_stepping, DynamicalCore())
    div_grid = get_prognostic_step(vars.grid.divergence, time_stepping, DynamicalCore())
    temp = get_prognostic_step(vars.grid.temperature, time_stepping, DynamicalCore())

    # use scratch array with zeros in case humidity doesn't exist
    humid = haskey(vars.grid, :humidity) ?
        get_prognostic_step(vars.grid.humidity, time_stepping, DynamicalCore()) :
        fill!(vars.scratch.grid.a, 0)

    (; pres_flux, pres_flux_sum_above, div_sum_above) = vars.dynamics
    (; temp_profile) = implicit

    # semi-implicit: terms here are explicit+implicit evaluated at time step i
    # implicit_correction! then calculated the implicit terms from Vi-1 minus Vi
    # to move the implicit terms to i-1 which is cheaper then the alternative below

    # Launch kernel to compute temperature tendency with adiabatic conversion
    arch = architecture(temp_tend_grid)
    launch!(
        arch, RingGridWorkOrder, size(temp_tend_grid), _temperature_tendency_kernel!,
        temp_tend_grid, temp, div_grid, humid, div_sum_above, pres_flux_sum_above,
        pres_flux, temp_profile, adiabatic_conversion.σ_lnp_A, adiabatic_conversion.σ_lnp_B, atmosphere
    )

    # write uT_anomaly_grid, vT_anomaly_grid (= u·T', v·T') for the flux divergence
    flux_grid_divergence!(get_step(vars.dynamics.grid.uT_anomaly), get_step(vars.dynamics.grid.vT_anomaly), temp, vars, model)
    return nothing
end

"""$(TYPEDSIGNATURES)

Compute the temperature tendency.

    ∂T/∂t += -∇⋅((u, v)*T') + T'D + κTᵥ*Dlnp/Dt

`+=` because the tendencies already contain parameterizations and vertical advection.
`T'` is the anomaly with respect to the reference/average temperature. Tᵥ is the virtual
temperature used in the adiabatic term κTᵥ*Dlnp/Dt.

Here, the previosuly computed gridded tendency contribution are accumulated and the spectral
tendency computed with divergence is computed via the `flux_spectral_divergence!` function."""
function temperature_spectral_tendency!(vars::Variables, S::AbstractSpectralTransform, time_stepping::AbstractTimeStepper)
    temp_tend = get_tendency_step(vars.tendencies.temperature, time_stepping, DynamicalCore())
    uT_spec = get_step(vars.dynamics.uT_anomaly)
    vT_spec = get_step(vars.dynamics.vT_anomaly)
    flux_spectral_divergence!(temp_tend, uT_spec, vT_spec, S; add = true, flipsign = true)
    return nothing
end

temperature_spectral_tendency!(vars::Variables, model::PrimitiveEquation) =
    temperature_spectral_tendency!(vars, model.spectral_transform, model.time_stepping)

@kernel inbounds = true function _temperature_tendency_kernel!(
        temp_tend_grid,             # Input/Output: temperature tendency
        temp_grid,                  # Input: temperature anomaly
        div_grid,                   # Input: divergence
        humid_grid,                 # Input: humidity
        div_sum_above,              # Input: sum of div from layers above
        pres_flux_sum_above,        # Input: sum of pres_flux from layers above
        pres_flux,                  # Input: (u,v)⋅∇lnp term
        temp_profile,               # Input: reference temperature profile
        σ_lnp_A,                    # Input: adiabatic conversion coefficient A
        σ_lnp_B,                    # Input: adiabatic conversion coefficient B
        atmosphere,                 # Input: atmosphere for κ and μ_virt_temp
    )

    ij, k = @index(Global, NTuple)
    Tₖ = temp_profile[k]    # average layer temperature from reference profile

    # coefficients from Simmons and Burridge 1981
    σ_lnp_A_k = σ_lnp_A[k]   # eq. 3.12, -1/Δσₖ*ln(σ_k+1/2/σ_k-1/2)
    σ_lnp_B_k = σ_lnp_B[k]   # eq. 3.12 -αₖ

    # Adiabatic conversion term following Simmons and Burridge 1981 but for σ coordinates
    # += as tend already contains parameterizations + vertical advection
    Tᵥ = virtual_temperature(temp_grid[ij, k] + Tₖ, humid_grid[ij, k], atmosphere)
    (; κ) = atmosphere
    temp_tend_grid[ij, k] +=
        temp_grid[ij, k] * div_grid[ij, k] +                # +T'D term of hori advection
        κ * Tᵥ * (                                          # +κTᵥ*Dlnp/Dt, adiabatic term
        σ_lnp_A_k * (div_sum_above[ij, k] + pres_flux_sum_above[ij, k]) +  # eq. 3.12 1st term
            σ_lnp_B_k * (div_grid[ij, k] + pres_flux[ij, k]) +             # eq. 3.12 2nd term
            pres_flux[ij, k]
    )                                                        # eq. 3.13

end

# no humidity tendency for dry core
humidity_tendency!(::Variables, ::PrimitiveDry) = nothing

"""$(TYPEDSIGNATURES)

Computes the gridded contributation to the humidity tendency `humid_tend_grid` via the `horizontal_grid_advection!`
Grid half of `humidity_tendency!`. Adds the `+q·div` advection term to `humid_tend_grid` and
writes the `(uq, vq)` flux intermediates to the grid-side named slots — no transform."""
function humidity_grid_tendency!(vars::Variables, model::PrimitiveWet)
    (; time_stepping) = model
    humid_tend_grid = get_tendency_step(vars.tendencies.grid.humidity, time_stepping, DynamicalCore())
    humid_grid = get_prognostic_step(vars.grid.humidity, time_stepping, DynamicalCore())
    horizontal_grid_advection!(humid_tend_grid, humid_grid, vars, model; add = true,
                               uA_grid = get_step(vars.dynamics.grid.uq),
                               vA_grid = get_step(vars.dynamics.grid.vq))
    return nothing
end
humidity_grid_tendency!(::Variables, ::PrimitiveDry) = nothing

"""$(TYPEDSIGNATURES)

Computes the spectral humidity tendency via the `horizontal_spectral_advection!`
Adds `-∇⋅(uq, vq)` to the previously computed gridded and transformed tendency."""
function humidity_spectral_tendency!(vars::Variables, model::PrimitiveWet)
    S = model.spectral_transform
    humid_tend = get_tendency_step(vars.tendencies.humidity, model.time_stepping, DynamicalCore())
    uq_spec = get_step(vars.dynamics.uq)
    vq_spec = get_step(vars.dynamics.vq)
    horizontal_spectral_advection!(humid_tend, uq_spec, vq_spec, S)
    return nothing
end
humidity_spectral_tendency!(::Variables, ::PrimitiveDry) = nothing

function tracer_advection!(
        vars::Variables,
        model::AbstractModel,
    )
    TS = model.time_stepping

    for (name, tracer) in model.tracers
        tracer_tend = get_tendency_step(vars.tendencies.tracers[name], TS, DynamicalCore())
        tracer_tend_grid = get_tendency_step(vars.tendencies.grid_tracers[name], TS, DynamicalCore())
        tracer_grid = get_prognostic_step(vars.grid.tracers[name], TS, DynamicalCore(), model)

        # add horizontal advection to parameterization + vertical advection + forcing/drag tendencies
        tracer.active && horizontal_advection!(tracer_tend, tracer_tend_grid, tracer_grid, vars, model, add = true)
    end
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
        model::AbstractModel;
        add::Bool = true,                   # use muladd (true) or overwrite (false) for the +A*div term
        uA_grid = get_step(vars.dynamics.grid.uT_anomaly),   # caller picks the correct named slot
        vA_grid = get_step(vars.dynamics.grid.vT_anomaly),
    )

    # barotropic model doesn't have divergence, the +A*div term is then zero
    if haskey(vars.grid, :divergence)
        div_grid = get_prognostic_step(vars.grid.divergence, model.time_stepping, DynamicalCore(), model)

        kernel_func = add ? muladd : @inline (a, b, c) -> a * b

        # Launch kernel to compute +A*div term of the advection operator
        arch = architecture(A_tend_grid)
        launch!(
            arch, RingGridWorkOrder, size(A_tend_grid), _horizontal_advection_kernel!,
            A_tend_grid, A_grid, div_grid, kernel_func
        )
    end

    # write u*A and v*A on grid
    flux_grid_divergence!(uA_grid, vA_grid, A_grid, vars, model)
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
        model::AbstractModel;
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
    flux_grid_divergence!(uA_grid, vA_grid, A_grid, vars, model)

    scratch_memory = vars.scratch.transform_memory
    S = model.spectral_transform
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
        model::AbstractModel,
    )
    u = get_prognostic_step(vars.grid.u, model.time_stepping, DynamicalCore(), model)
    v = get_prognostic_step(vars.grid.v, model.time_stepping, DynamicalCore(), model)
    (; coslat⁻¹) = model.geometry
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
        uA_grid,        # Output: u*A on grid
        vA_grid,        # Output: v*A on grid
        A_grid,         # Input: field to be advected
        u_grid,         # Input: zonal velocity
        v_grid,         # Input: meridional velocity
        coslat⁻¹,       # Input: 1/cos(latitude) for scaling
        whichring,      # Input: mapping from grid point to latitude ring
    )
    I = @index(Global, Cartesian)

    j = whichring[I[1]]             # latitude ring index for this grid point
    coslat⁻¹j = coslat⁻¹[j]         # get coslat⁻¹ for this latitude
    Acoslat⁻¹j = A_grid[I] * coslat⁻¹j
    uA_grid[I] = u_grid[I] * Acoslat⁻¹j
    vA_grid[I] = v_grid[I] * Acoslat⁻¹j
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
        model::AbstractModel,
    )
    (; f) = model.coriolis
    scale = vars.prognostic.scale[]     # used to scale Coriolis f on the fly, as it's being added to a scaled vorticity
    (; coslat⁻¹) = model.geometry
    time_stepping = model.time_stepping
    u_tend_grid = get_tendency_step(vars.tendencies.grid.u, time_stepping, DynamicalCore())         # already contains forcing
    v_tend_grid = get_tendency_step(vars.tendencies.grid.v, time_stepping, DynamicalCore())         # already contains forcing
    u = get_prognostic_step(vars.grid.u, time_stepping, DynamicalCore(), model)
    v = get_prognostic_step(vars.grid.v, time_stepping, DynamicalCore(), model)
    vor = get_prognostic_step(vars.grid.vorticity, time_stepping, DynamicalCore(), model)

    (; whichring) = u_tend_grid.grid

    arch = architecture(u_tend_grid)
    launch!(
        architecture(u), RingGridWorkOrder, size(u), _vorticity_flux_kernel!,
        u_tend_grid, v_tend_grid, u, v, vor, f, scale, coslat⁻¹, whichring
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
        S::AbstractSpectralTransform,
        time_stepper::AbstractTimeStepper;
        div::Bool = true,
        add::Bool = false,
    )
    vor_tend = get_tendency_step(vars.tendencies.vorticity, time_stepper, DynamicalCore())
    u_tend = get_tendency_step(vars.dynamics.u_tendency, time_stepper, DynamicalCore())
    v_tend = get_tendency_step(vars.dynamics.v_tendency, time_stepper, DynamicalCore())

    curl!(vor_tend, u_tend, v_tend, S; add)                   # ∂ζ/∂t = ∇×(u_tend, v_tend)

    if div                                                   # not needed/available in barotropic model
        div_tend = get_tendency_step(vars.tendencies.divergence, time_stepper, DynamicalCore())
        divergence!(div_tend, u_tend, v_tend, S; add)        # ∂D/∂t = ∇⋅(u_tend, v_tend)
    end
    return nothing
end

@kernel inbounds = true function _vorticity_flux_kernel!(
        u_tend_grid, v_tend_grid, u, v, vor, f, scale, coslat⁻¹, whichring
    )
    # Get indices
    ij, k = @index(Global, NTuple)
    j = whichring[ij]

    # Get the coriolis parameter and cosine latitude factor for this latitude
    f_j = f[j] * scale              # scaled as is vorticity
    coslat⁻¹j = coslat⁻¹[j]

    # Calculate absolute vorticity, scale f as is vorticity too on the fly
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
    vorticity_flux_curldiv!(vars, model, div = true, add = true)

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
    vorticity_flux_curldiv!(vars, model, div = false, add = true)

#TODO: OLD VERSION / SEQUENTIAL VERSION: MIGHT BE DELETED
function bernoulli_potential!(vars::Variables, model::ShallowWater)
    S = model.spectral_transform
    bernoulli_grid_potential!(vars, model, model.time_stepping)

    bernoulli = get_step(vars.dynamics.kinetic_energy)
    bernoulli_grid = get_step(vars.dynamics.grid.kinetic_energy)
    scratch_memory = vars.scratch.transform_memory
    transform!(bernoulli, bernoulli_grid, scratch_memory, S)

    bernoulli_spectral_potential!(vars, model, model.time_stepping)
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
        TS::AbstractTimeStepper,
    )
    bernoulli_grid_potential!(vars, S, TS)

    bernoulli = get_step(vars.dynamics.kinetic_energy)
    bernoulli_grid = get_step(vars.dynamics.grid.kinetic_energy)
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

    bernoulli_spectral_potential!(vars, S, TS)
    return nothing
end

"""$(TYPEDSIGNATURES)
Gridded contribution of `bernoulli_potential!` for ShallowWater: writes `kinetic_energy_grid = ½(u²+v²)+Φ`
(geopotential Φ is included because in ShallowWater geopotential lives only on the grid)."""
function bernoulli_grid_potential!(vars::Variables, model::ShallowWater, time_stepping::AbstractTimeStepper)
    u = get_prognostic_step(vars.grid.u, time_stepping, DynamicalCore(), model)
    v = get_prognostic_step(vars.grid.v, time_stepping, DynamicalCore(), model)
    Φ = vars.dynamics.geopotential
    bernoulli_grid = get_step(vars.dynamics.grid.kinetic_energy)
    half = convert(eltype(bernoulli_grid), 0.5)
    @. bernoulli_grid = half * (u^2 + v^2) + Φ
    return nothing
end

"""$(TYPEDSIGNATURES)
Gridded contribution of `bernoulli_potential!` for PrimitiveEquation: writes `kinetic_energy_grid = ½(u²+v²)`."""
function bernoulli_grid_potential!(vars::Variables, ::Union{PrimitiveEquation, AbstractSpectralTransform}, time_stepping::AbstractTimeStepper)
    u = get_prognostic_step(vars.grid.u, time_stepping, DynamicalCore())
    v = get_prognostic_step(vars.grid.v, time_stepping, DynamicalCore())
    bernoulli_grid = get_step(vars.dynamics.grid.kinetic_energy)
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
    _bernoulli_spectral_potential!(vars, model.spectral_transform, model.time_stepping)
    return nothing
end

function bernoulli_spectral_potential!(vars::Variables, model::PrimitiveEquation)
    bernoulli_spectral_potential!(vars, model.spectral_transform, model.time_stepping)
    return nothing
end

function bernoulli_spectral_potential!(vars::Variables, S::AbstractSpectralTransform, time_stepper::AbstractTimeStepper)
    # PrimitiveEquation path: add spectral geopotential to KE before the Laplacian.
    bernoulli = get_step(vars.dynamics.kinetic_energy)
    geopot = vars.dynamics.spectral_geopotential
    bernoulli .+= geopot
    _bernoulli_spectral_potential!(vars, S, time_stepper)
    return nothing
end

function _bernoulli_spectral_potential!(vars::Variables, S::AbstractSpectralTransform, time_stepper::AbstractTimeStepper)
    bernoulli = get_step(vars.dynamics.kinetic_energy)
    div_tend = get_tendency_step(vars.tendencies.divergence, time_stepper, DynamicalCore()) 
    ∇²!(div_tend, bernoulli, S, add = true, flipsign = true)
    return nothing
end


"""$(TYPEDSIGNATURES)
Computes the (negative) divergence of the volume fluxes `uh, vh` for the continuity equation, -∇⋅(uh, vh)."""
function volume_flux_divergence!(
        vars::Variables,
        model::ShallowWater,
    )

    η = get_prognostic_step(vars.grid.η, model.time_stepping, ContinuityEquation(), model)
    η_tend = get_tendency_step(vars.tendencies.η, model.time_stepping, ContinuityEquation())
    (; orography) = model.orography
    H = model.atmosphere.layer_thickness

    # compute dynamic layer thickness h on the grid
    # η is the interface displacement, update to layer thickness h = η + H - Hb
    # H is the layer thickness at rest without mountains, Hb the orography
    # TODO this leaves η <- h between here and the transforms after the time stepping
    # change to h = η + H - Hb here using a scratch array for h?
    η .+= H .- orography

    # now do -∇⋅(uh, vh) and store in η_tend
    flux_divergence!(η_tend, η, vars, model, add = true, flipsign = true)
    return nothing
end


"""$(TYPEDSIGNATURES)
Calculates the average temperature of a layer from the l=m=0 harmonic
and stores the result in `diagn.temp_average`"""
function temperature_average!(
        vars::Variables,
        temp::LowerTriangularArray,
        S::AbstractSpectralTransform,
    )
    # average from l=m=0 harmonic divided by norm of the sphere
    @. vars.dynamics.average_temperature_profile = real(temp[1, :]) / S.norm_sphere
    return nothing
end

function reset_tendencies!(vars::Variables, time_stepping::AbstractTimeStepper; value = 0)
    _reset_tendencies!(vars.tendencies, time_stepping, value)
    return vars
end

# recursively fill all arrays in a NamedTuple, unpacking nested NamedTuples
# this avoids Union-typed iteration which Enzyme cannot differentiate
@inline _reset_tendencies!(nt::NamedTuple, ts, value) = _reset_tendencies_inner!(values(nt), ts, value)
@inline _reset_tendencies_inner!(::Tuple{}, _, _) = nothing
@inline function _reset_tendencies_inner!(t::Tuple, ts, value)
    _reset_tendency!(first(t), ts, value)
    _reset_tendencies_inner!(Base.tail(t), ts, value)
    return nothing
end

# dispatch on element type: nested NamedTuple vs array
@inline _reset_tendency!(nt::NamedTuple, ts, value) = _reset_tendencies_inner!(values(nt), ts, value)
@inline function _reset_tendency!(a::AbstractArray, time_stepping, value)
    a_step = get_tendency_step(a, time_stepping, ResetTendencies())
    return fill!(a_step, value)
end
