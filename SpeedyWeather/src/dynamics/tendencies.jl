"""$(TYPEDSIGNATURES)
Calculate all tendencies for the BarotropicModel."""
function dynamics_tendencies!(
        vars::Variables,
        lf::Integer,                    # leapfrog index to evaluate tendencies at
        model::Barotropic,
    )
    forcing!(vars, lf, model)           # = (Fᵤ, Fᵥ) forcing for u, v
    drag!(vars, lf, model)              # drag term for u, v
    vorticity_flux!(vars, model)        # = ∇×(v(ζ+f) + Fᵤ, -u(ζ+f) + Fᵥ)
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

    # = ∇×(v(ζ+f) + Fᵤ, -u(ζ+f) + Fᵥ), tendency for vorticity
    # = ∇⋅(v(ζ+f) + Fᵤ, -u(ζ+f) + Fᵥ), tendency for divergence
    vorticity_flux!(vars, model)

    geopotential!(vars, planet)             # geopotential Φ = gη in shallow water
    bernoulli_potential!(vars, model)       # = -∇²(E+Φ), tendency for divergence

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

    # ∂ln(pₛ)/∂t = -(ū, v̄)⋅∇ln(pₛ) - D̄
    surface_pressure_tendency!(vars, spectral_transform)

    # calculate vertical velocity σ̇ in sigma coordinates for the vertical mass flux M = pₛ * σ̇
    vertical_velocity!(vars, geometry)

    # add the RTₖlnpₛ term to geopotential
    linear_pressure_gradient!(vars, lf_implicit, atmosphere, implicit)

    # use σ̇ for the vertical advection of u, v, T, q
    vertical_advection!(vars, model)

    # vorticity advection, pressure gradient term
    vordiv_tendencies!(vars, model)

    # hor. advection + adiabatic term
    temperature_tendency!(vars, model)

    # horizontal advection of humidity (nothing for wetcore)
    humidity_tendency!(vars, model)

    # add -∇²(E + ϕ + RTₖlnpₛ) term to div tendency
    bernoulli_potential!(vars, spectral_transform)

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
    dpres_dx_spec = vars.scratch.a_2D           # reuse 2D work arrays for gradients
    dpres_dy_spec = vars.scratch.b_2D           # in spectral space
    (; dpres_dx, dpres_dy) = vars.dynamics      # but store in grid space

    ∇!(dpres_dx_spec, dpres_dy_spec, pres, S)                                       # CALCULATE ∇ln(pₛ)
    transform!(dpres_dx, dpres_dx_spec, scratch_memory, S, unscale_coslat = true)   # transform to grid: zonal gradient
    transform!(dpres_dy, dpres_dy_spec, scratch_memory, S, unscale_coslat = true)   # meridional gradient

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
    pres_tend = vars.tendencies.pressure
    pres_tend_grid = vars.tendencies.grid.pressure
    (; dpres_dx, dpres_dy, u_mean_grid, v_mean_grid, div_mean) = vars.dynamics
    scratch_memory = vars.scratch.transform_memory

    # in grid-point space the the (ū, v̄)⋅∇lnpₛ term (swap sign in spectral)
    # += to allow for forcing contributions already in pres_tend_grid
    @. pres_tend_grid += u_mean_grid * dpres_dx + v_mean_grid * dpres_dy

    ūv̄∇lnpₛ = vars.scratch.a_2D             # reuse 2D work array
    transform!(ūv̄∇lnpₛ, pres_tend_grid, scratch_memory, S)

    # for semi-implicit div_mean is calc at time step i-1 in vertical_integration!
    @. pres_tend -= ūv̄∇lnpₛ + div_mean      # add the -div_mean term in spectral, swap sign

    pres_tend.data[1:1] .= 0                # for mass conservation
    return nothing
end

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
        implicit::ImplicitPrimitiveEquation,
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
        implicit::ImplicitPrimitiveEquation,
        S::AbstractSpectralTransform,
    )
    (; f) = coriolis                            # coriolis parameter
    Tₖ = implicit.temp_profile                  # reference temperature profile
    (; coslat⁻¹) = geometry

    # tendencies already contain parameterizations + advection, therefore accumulate
    u_tend_grid = vars.tendencies.grid.u
    v_tend_grid = vars.tendencies.grid.v
    (; u, v) = vars.grid                                    # velocities, vorticity, temperature
    vor = vars.grid.vorticity
    temp = vars.grid.temperature
    vars.scratch.grid.a .= 0
    humid = haskey(vars.grid, :humidity) ? vars.grid.humidity : vars.scratch.grid.a
    (; dpres_dx, dpres_dy) = vars.dynamics              # zonal/meridional gradient of logarithm of surface pressure
    scratch_memory = vars.scratch.transform_memory

    # Launch kernel to compute u_tend and v_tend with vorticity flux and pressure gradient
    (; whichring) = u_tend_grid.grid            # precomputed ring indices
    arch = architecture(u_tend_grid)
    launch!(
        arch, RingGridWorkOrder, size(u_tend_grid), _vordiv_tendencies_kernel!,
        u_tend_grid, v_tend_grid, u, v, vor, temp, humid,
        dpres_dx, dpres_dy, Tₖ, f, coslat⁻¹, whichring, atmosphere,
    )
    # divergence and curl of that u, v_tend vector for vor, div tendencies
    vor_tend = vars.tendencies.vorticity
    div_tend = vars.tendencies.divergence
    u_tend = vars.scratch.a
    v_tend = vars.scratch.b

    transform!(u_tend, u_tend_grid, scratch_memory, S)
    transform!(v_tend, v_tend_grid, scratch_memory, S)

    curl!(vor_tend, u_tend, v_tend, S, add = true)            # ∂ζ/∂t += ∇×(u_tend, v_tend)
    divergence!(div_tend, u_tend, v_tend, S, add = true)      # ∂D/∂t += ∇⋅(u_tend, v_tend)
    return nothing
end

@kernel inbounds = true function _vordiv_tendencies_kernel!(
        u_tend_grid,            # Input/Output: zonal wind tendency
        v_tend_grid,            # Input/Output: meridional wind tendency
        u_grid,                 # Input: zonal velocity
        v_grid,                 # Input: meridional velocity
        vor_grid,               # Input: relative vorticity
        temp_grid,              # Input: temperature anomaly
        humid_grid,             # Input: humidity
        dpres_dx,                 # Input: zonal gradient of log surface pressure
        dpres_dy,                 # Input: meridional gradient of log surface pressure
        Tₖ,                     # Input: reference temperature profile
        f,              # Input: coriolis parameter
        coslat⁻¹,       # Input: 1/cos(latitude) for scaling
        whichring,      # Input: mapping from grid point to latitude ring
        atmosphere,             # Input: atmosphere for R_dry and μ_virt_temp
    )
    ij, k = @index(Global, NTuple)
    j = whichring[ij]           # latitude ring index for this grid point
    coslat⁻¹j = coslat⁻¹[j]     # get coslat⁻¹ for this latitude
    f_j = f[j]                  # coriolis parameter for this latitude

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
    u_tend = vars.scratch.a
    v_tend = vars.scratch.b

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
Compute the temperature tendency.

    ∂T/∂t += -∇⋅((u, v)*T') + T'D + κTᵥ*Dlnp/Dt

`+=` because the tendencies already contain parameterizations and vertical advection.
`T'` is the anomaly with respect to the reference/average temperature. Tᵥ is the virtual
temperature used in the adiabatic term κTᵥ*Dlnp/Dt."""
function temperature_tendency!(
        vars::Variables,
        adiabatic_conversion::AbstractAdiabaticConversion,
        atmosphere::AbstractAtmosphere,
        implicit::ImplicitPrimitiveEquation,
        G::Geometry,
        S::AbstractSpectralTransform,
    )
    temp_tend = vars.tendencies.temperature
    temp_tend_grid = vars.tendencies.grid.temperature
    div_grid = vars.grid.divergence
    temp = vars.grid.temperature

    # use scratch array with zeros in case humidity doesn't exist
    vars.scratch.grid.a .= 0
    humid = haskey(vars.grid, :humidity) ? vars.grid.humidity : vars.scratch.grid.a

    (; pres_flux, pres_flux_sum_above, div_sum_above) = vars.dynamics
    scratch_memory = vars.scratch.transform_memory
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

    transform!(temp_tend, temp_tend_grid, scratch_memory, S)

    # now add the -∇⋅((u, v)*T') term
    flux_divergence!(temp_tend, temp, vars, G, S, add = true, flipsign = true)
    return nothing
end

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

function humidity_tendency!(
        vars::Variables,
        model::PrimitiveWet
    )
    G = model.geometry
    S = model.spectral_transform

    humid_tend = vars.tendencies.humidity
    humid_tend_grid = vars.tendencies.grid.humidity
    humid = vars.grid.humidity

    # add horizontal advection to parameterization + vertical advection tendencies
    horizontal_advection!(humid_tend, humid_tend_grid, humid, vars, G, S, add = true)

    return nothing
end

# no humidity tendency for dry core
humidity_tendency!(::Variables, ::PrimitiveDry) = nothing

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
Compute the horizontal advection"""
function horizontal_advection!(
        A_tend::LowerTriangularArray,       # Output: tendency to write into
        A_tend_grid::AbstractField,         # Input: tendency incl prev terms
        A_grid::AbstractField,              # Input: grid field to be advected
        vars::Variables,
        G::Geometry,
        S::AbstractSpectralTransform;
        add::Bool = true,                   # add/overwrite A_tend_grid?
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

    scratch_memory = vars.scratch.transform_memory
    transform!(A_tend, A_tend_grid, scratch_memory, S)  # for +A*div in spectral space

    # now add the -∇⋅((u, v)*A) term
    flux_divergence!(A_tend, A_grid, vars, G, S, add = true, flipsign = true)

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
    )
    (; u, v) = vars.grid
    scratch_memory = vars.scratch.transform_memory
    (; coslat⁻¹) = G

    # reuse general work arrays a, b, a_grid, b_grid
    uA = vars.scratch.a                 # = u*A in spectral
    vA = vars.scratch.b                 # = v*A in spectral
    uA_grid = vars.scratch.grid.a       # = u*A on grid
    vA_grid = vars.scratch.grid.b       # = v*A on grid

    # Launch kernel to compute u*A and v*A with coslat scaling
    (; whichring) = A_grid.grid         # precomputed ring indices
    arch = architecture(A_grid)
    launch!(
        arch, RingGridWorkOrder, size(A_grid), _flux_divergence_kernel!,
        uA_grid, vA_grid, A_grid, u, v, coslat⁻¹, whichring
    )

    transform!(uA, uA_grid, scratch_memory, S)
    transform!(vA, vA_grid, scratch_memory, S)

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

    (; f) = coriolis
    (; coslat⁻¹) = geometry

    u_tend_grid = vars.tendencies.grid.u                # already contains forcing
    v_tend_grid = vars.tendencies.grid.v                # already contains forcing
    (; u, v) = vars.grid                                  # velocities and vorticity on grid
    vor = vars.grid.vorticity
    (; whichring) = u.grid                              # precomputed ring indices
    scratch_memory = vars.scratch.transform_memory      # scratch memory for transforms

    # Launch the kernel for vorticity flux calculation
    arch = S.architecture

    launch!(
        arch, RingGridWorkOrder, size(u), _vorticity_flux_kernel!,
        u_tend_grid, v_tend_grid, u, v, vor, f, coslat⁻¹, whichring
    )

    # divergence and curl of that u, v_tend vector for vor, div tendencies
    vor_tend = vars.tendencies.vorticity
    u_tend = vars.scratch.a
    v_tend = vars.scratch.b

    transform!(u_tend, u_tend_grid, scratch_memory, S)
    transform!(v_tend, v_tend_grid, scratch_memory, S)

    curl!(vor_tend, u_tend, v_tend, S; add)                 # ∂ζ/∂t = ∇×(u_tend, v_tend)

    if div                                                  # not needed/availbel in barotropic model
        div_tend = vars.tendencies.divergence
        divergence!(div_tend, u_tend, v_tend, S; add)       # ∂D/∂t = ∇⋅(u_tend, v_tend)
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

function bernoulli_potential!(vars::Variables, model::ShallowWater)
    S = model.spectral_transform
    scratch_memory = vars.scratch.transform_memory
    (; u, v) = vars.grid
    Φ = vars.grid.geopotential
    bernoulli = vars.scratch.a                                  # reuse work arrays a, a_grid
    bernoulli_grid = vars.scratch.grid.a
    div_tend = vars.tendencies.divergence

    half = convert(eltype(bernoulli_grid), 0.5)
    @. bernoulli_grid = half * (u^2 + v^2) + Φ
    transform!(bernoulli, bernoulli_grid, scratch_memory, S)    # to spectral space
    ∇²!(div_tend, bernoulli, S, add = true, flipsign = true)    # add -∇²(½(u² + v²) + ϕ)
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
    (; u, v) = vars.grid
    scratch_memory = vars.scratch.transform_memory
    geopot = vars.dynamics.geopotential
    bernoulli = vars.scratch.a                              # reuse work arrays a, a_grid
    bernoulli_grid = vars.scratch.grid.a
    div_tend = vars.tendencies.divergence

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

    bernoulli_grid .= 1 // 2 .* (u.^2 + v.^2)                    # = ½(u² + v²) on grid
    transform!(bernoulli, bernoulli_grid, scratch_memory, S)    # to spectral space
    bernoulli .+= geopot                                        # add geopotential Φ
    ∇²!(div_tend, bernoulli, S, add = true, flipsign = true)    # add -∇²(½(u² + v²) + ϕ)
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
