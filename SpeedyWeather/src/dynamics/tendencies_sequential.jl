# SEQUENTIAL versions of the dynamical core tendency methods 
# Currently the dycore on all architectures fuses the tendency computations, 
# here we provide sequential implementations of the tendency methods. 
# Those might be used for testing and debugging and potentially in the future 
# for a more memory optimized CPU-version. 


"""
$(TYPEDSIGNATURES)

SEQUENTIAL VERSION 

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
        time_stepping::AbstractTimeStepper,
    )
    surface_pressure_grid_tendency!(vars, time_stepping)

    pres_tend = get_tendency_step(vars.tendencies.pressure, time_stepping, DynamicalCore())
    pres_tend_grid = get_tendency_step(vars.tendencies.grid.pressure, time_stepping, DynamicalCore())
    scratch_memory = vars.scratch.transform_memory
    transform!(pres_tend, pres_tend_grid, scratch_memory, S)

    surface_pressure_spectral_tendency!(vars)
    return nothing
end


"""$(TYPEDSIGNATURES)

SEQUENTIAL VERSION

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
        time_stepping::AbstractTimeStepper,
    )
    vordiv_grid_tendencies!(vars, coriolis, atmosphere, geometry, implicit, time_stepping)

    u_tend_grid = get_tendency_step(vars.tendencies.grid.u, time_stepping, DynamicalCore())
    v_tend_grid = get_tendency_step(vars.tendencies.grid.v, time_stepping, DynamicalCore()) 
    u_tend = get_tendency_step(vars.dynamics.u_tendency, time_stepping, DynamicalCore()) 
    v_tend = get_tendency_step(vars.dynamics.v_tendency, time_stepping, DynamicalCore()) 
    scratch_memory = vars.scratch.transform_memory
    transform!(u_tend, u_tend_grid, scratch_memory, S)
    transform!(v_tend, v_tend_grid, scratch_memory, S)

    vordiv_spectral_tendencies!(vars, time_stepping, S)
    return nothing
end


"""$(TYPEDSIGNATURES)

SEQUENTIAL VERSION

Compute the temperature tendency.

    ∂T/∂t += -∇⋅((u, v)*T') + T'D + κTᵥ*Dlnp/Dt

`+=` because the tendencies already contain parameterizations and vertical advection.
`T'` is the anomaly with respect to the reference/average temperature. Tᵥ is the virtual
temperature used in the adiabatic term κTᵥ*Dlnp/Dt."""
function temperature_tendency!(
        vars::Variables,
        model::PrimitiveEquation,
    )
    (; spectral_transform, time_stepping) = model

    temperature_grid_tendency!(vars, model)

    temp_tend = get_tendency_step(vars.tendencies.temperature, time_stepping, DynamicalCore())
    temp_tend_grid = get_tendency_step(vars.tendencies.grid.temperature, time_stepping, DynamicalCore())
    uT_grid = get_step(vars.dynamics.grid.uT_anomaly)
    vT_grid = get_step(vars.dynamics.grid.vT_anomaly)
    uT_spec = get_step(vars.dynamics.uT_anomaly)
    vT_spec = get_step(vars.dynamics.vT_anomaly)
    scratch_memory = vars.scratch.transform_memory
    transform!(temp_tend, temp_tend_grid, scratch_memory, spectral_transform)
    transform!(uT_spec, uT_grid, scratch_memory, spectral_transform)
    transform!(vT_spec, vT_grid, scratch_memory, spectral_transform)

    temperature_spectral_tendency!(vars, spectral_transform, time_stepping)
    return nothing
end


#TODO: OLD VERSION / SEQUENTIAL VERSION: MIGHT BE DELETED
function humidity_tendency!(
        vars::Variables,
        model::PrimitiveWet
    )
    (; spectral_transform, time_stepping) = model

    humidity_grid_tendency!(vars, model)

    humid_tend = get_tendency_step(vars.tendencies.humidity, time_stepping, DynamicalCore())
    humid_tend_grid = get_tendency_step(vars.tendencies.grid.humidity, time_stepping, DynamicalCore())
    uq_grid = get_step(vars.dynamics.grid.uq)
    vq_grid = get_step(vars.dynamics.grid.vq)
    uq_spec = get_step(vars.dynamics.uq)
    vq_spec = get_step(vars.dynamics.vq)
    scratch_memory = vars.scratch.transform_memory
    transform!(humid_tend, humid_tend_grid, scratch_memory, spectral_transform)
    transform!(uq_spec, uq_grid, scratch_memory, spectral_transform)
    transform!(vq_spec, vq_grid, scratch_memory, spectral_transform)

    humidity_spectral_tendency!(vars, model)
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
        model::AbstractModel;
        add::Bool = true,                   # add/overwrite A_tend_grid?
    # Forwarded to `flux_divergence!`; defaults to the unfused scratch.{a,b} for unfused callers.
        uA = vars.scratch.a,
        vA = vars.scratch.b,
        uA_grid = vars.scratch.grid.a,
        vA_grid = vars.scratch.grid.b,
    )
    (; spectral_transform) = model

    horizontal_grid_advection!(A_tend_grid, A_grid, vars, model; add, uA_grid, vA_grid)

    scratch_memory = vars.scratch.transform_memory
    transform!(A_tend, A_tend_grid, scratch_memory, spectral_transform)  # for +A*div in spectral space
    transform!(uA, uA_grid, scratch_memory, spectral_transform)
    transform!(vA, vA_grid, scratch_memory, spectral_transform)

    horizontal_spectral_advection!(A_tend, uA, vA, spectral_transform)
    return nothing
end


"""
$(TYPEDSIGNATURES)

SEQUENTIAL VERSION

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
        model::AbstractModel;
        div::Bool = true,     # also calculate div of vor flux?
        add::Bool = false,    # accumulate in vor/div tendencies?
    )
    vorticity_flux_grid_tendencies!(vars, model)

    time_stepping = model.time_stepping
    S = model.spectral_transform
    u_tend_grid = get_tendency_step(vars.tendencies.grid.u, time_stepping, DynamicalCore())
    v_tend_grid = get_tendency_step(vars.tendencies.grid.v, time_stepping, DynamicalCore())
    u_tend = get_tendency_step(vars.dynamics.u_tendency, time_stepping, DynamicalCore())
    v_tend = get_tendency_step(vars.dynamics.v_tendency, time_stepping, DynamicalCore())
    scratch_memory = vars.scratch.transform_memory
    transform!(u_tend, u_tend_grid, scratch_memory, S)
    transform!(v_tend, v_tend_grid, scratch_memory, S)

    vorticity_flux_spectral_tendencies!(vars, S, time_stepping; div, add)
    return nothing
end

"""$(TYPEDSIGNATURES)

SEQUENTIAL VERSION

Computes the (negative) divergence of the volume fluxes `uh, vh` for the continuity
equation, `η_tend -= ∇⋅(uh, vh)`, with its own grid→spectral transforms of the fluxes."""
function volume_flux_divergence!(vars::Variables, model::ShallowWater)
    volume_flux_divergence_grid!(vars, model)

    S = model.spectral_transform
    scratch_memory = vars.scratch.transform_memory
    transform!(get_step(vars.dynamics.uh), get_step(vars.dynamics.grid.uh), scratch_memory, S)
    transform!(get_step(vars.dynamics.vh), get_step(vars.dynamics.grid.vh), scratch_memory, S)

    volume_flux_divergence_spectral!(vars, model)
    return nothing
end

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
