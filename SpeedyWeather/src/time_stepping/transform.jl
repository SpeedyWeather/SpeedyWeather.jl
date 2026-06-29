# dispatch over time stepping to decide which prognostic steps to transform to grid space
SpeedyTransforms.transform!(vars::Variables, model::AbstractModel; kwargs...) =
    transform!(vars, model.time_stepping, model; kwargs...)

"""$(TYPEDSIGNATURES)
Propagate the spectral state of the prognostic variables of `vars` to the
grid variables `vars` for the barotropic vorticity model.
Updates grid vorticity, spectral stream function and spectral and grid velocities u, v."""
function SpeedyTransforms.transform!(
        vars::Variables,
        time_stepping::AbstractTimeStepper,
        model::Barotropic;
        initialize::Bool = false,
    )
    initialize && initialize!(vars, time_stepping, model)

    S = model.spectral_transform
    # model is passed on to get_prognostic_step as for 2D models Leapfrog choses the 1st step
    # not the non-existing 2nd step though both represent the current time step (not the previous one)
    u_grid = get_prognostic_step(vars.grid.u, time_stepping, S, model)
    v_grid = get_prognostic_step(vars.grid.v, time_stepping, S, model)
    vor_grid = get_prognostic_step(vars.grid.vorticity, time_stepping, S, model)

    # U = u*coslat, V=v*coslat
    U = vars.scratch.a      # reuse scratch arrays for velocities in spectral
    V = vars.scratch.b      # reuse scratch arrays for velocities in spectral
    vor = get_prognostic_step(vars.prognostic.vorticity, time_stepping, S)

    # get vorticity on grid from spectral vor
    scratch_memory = vars.scratch.transform_memory
    transform!(vor_grid, vor, scratch_memory, S)

    # get spectral U, V from spectral vorticity via stream function Ψ
    # U = u*coslat = -coslat*∂Ψ/∂lat
    # V = v*coslat = ∂Ψ/∂lon, radius omitted in both cases
    UV_from_vor!(U, V, vor, S)

    # transform from U, V in spectral to u, v on grid (U, V = u, v*coslat)
    transform!(u_grid, U, scratch_memory, S, unscale_coslat = true)
    transform!(v_grid, V, scratch_memory, S, unscale_coslat = true)

    for (name, tracer) in model.tracers
        tracer_var = get_prognostic_step(vars.prognostic.tracers[name], time_stepping, S, model)
        tracer_grid = get_prognostic_step(vars.grid.tracers[name], time_stepping, S, model)
        tracer.active && transform!(tracer_grid, tracer_var, scratch_memory, S)
    end

    # transform random pattern for random process unless random_process=nothing
    transform!(vars, model.random_process, S)

    return nothing
end

"""
$(TYPEDSIGNATURES)
Propagate the spectral state of the prognostic variables of `vars` to the
grid variables in `vars` for the shallow water model. Updates grid vorticity,
grid divergence, grid interface displacement (`pres_grid`) and the velocities
u, v."""
function SpeedyTransforms.transform!(
        vars::Variables,
        time_stepping::AbstractTimeStepper,
        model::ShallowWater;
        initialize::Bool = false,
    )
    initialize && initialize!(vars, time_stepping, model)

    S = model.spectral_transform
    # model is passed on to get_prognostic_step as for 2D models Leapfrog choses the 1st step
    # not the non-existing 2nd step though both represent the current time step (not the previous one)
    u_grid = get_prognostic_step(vars.grid.u, time_stepping, S, model)
    v_grid = get_prognostic_step(vars.grid.v, time_stepping, S, model)
    vor_grid = get_prognostic_step(vars.grid.vorticity, time_stepping, S, model)
    div_grid = get_prognostic_step(vars.grid.divergence, time_stepping, S, model)
    η_grid = get_prognostic_step(vars.grid.η, time_stepping, S, model)

    vor = get_prognostic_step(vars.prognostic.vorticity, time_stepping, S)
    div = get_prognostic_step(vars.prognostic.divergence, time_stepping, S)
    η = get_prognostic_step(vars.prognostic.η, time_stepping, S)

    # U = u*coslat, V=v*coslat
    U = vars.scratch.a
    V = vars.scratch.b

    scratch_memory = vars.scratch.transform_memory

    transform!(vor_grid, vor, scratch_memory, S)    # get vorticity on grid from spectral vor
    transform!(div_grid, div, scratch_memory, S)    # get divergence on grid from spectral div
    transform!(η_grid, η, scratch_memory, S)        # get η on grid from spectral η

    # get spectral U, V from vorticity and divergence via stream function Ψ and vel potential ϕ
    # U = u*coslat = -coslat*∂Ψ/∂lat + ∂ϕ/dlon
    # V = v*coslat =  coslat*∂ϕ/∂lat + ∂Ψ/dlon
    UV_from_vordiv!(U, V, vor, div, S)

    # transform from U, V in spectral to u, v on grid (U, V = u, v*coslat)
    transform!(u_grid, U, scratch_memory, S, unscale_coslat = true)
    transform!(v_grid, V, scratch_memory, S, unscale_coslat = true)

    for (name, tracer) in model.tracers
        tracer_var = get_prognostic_step(vars.prognostic.tracers[name], time_stepping, S, model)
        tracer_grid = get_prognostic_step(vars.grid.tracers[name], time_stepping, S, model)
        tracer.active && transform!(tracer_grid, tracer_var, scratch_memory, S)
    end

    # transform random pattern for random process unless random_process=nothing
    transform!(vars, model.random_process, S)

    return nothing
end

"""$(TYPEDSIGNATURES)
Propagate the spectral state of the prognostic variables of `vars` to the
grid variables in `vars` for primitive equation models. Updates grid vorticity,
grid divergence, grid temperature, pressure (`pres_grid`) and the velocities
u, v."""
function SpeedyTransforms.transform!(
        vars::Variables,
        time_stepping::AbstractTimeStepper,
        model::PrimitiveEquation;
        initialize::Bool = false,
    )

    # used to copy 1st step to 2nd step for leapfrog to always transform
    # the 2nd step to grid
    initialize && initialize!(vars, time_stepping, model)

    S = model.spectral_transform

    vor = get_prognostic_step(vars.prognostic.vorticity, time_stepping, S)
    div = get_prognostic_step(vars.prognostic.divergence, time_stepping, S)
    temp = get_prognostic_step(vars.prognostic.temperature, time_stepping, S)

    scratch_memory = vars.scratch.transform_memory
    U = vars.scratch.a                              # reuse work arrays
    V = vars.scratch.b                              # U = u*coslat, V=v*coslat

    # retain previous time step for vertical advection and parameterizations
    # if not initial step do before transforms i.e before that step is overwritten
    initialize || move_prognostic_grid_variables_back!(vars, time_stepping, model)

    # get spectral U, V from vorticity and divergence via stream function Ψ and vel potential ϕ
    # U = u*coslat = -coslat*∂Ψ/∂lat + ∂ϕ/dlon
    # V = v*coslat =  coslat*∂ϕ/∂lat + ∂Ψ/dlon
    UV_from_vordiv!(U, V, vor, div, S)

    # Batched spec→grid for the prognostic state: one call covers vorticity, divergence,
    # temperature, pressure (and humidity for PrimitiveWet).
    prog_parent = parent(vars.fused.prognostic)
    grid_parent = parent(vars.fused.grid)
    transform!(get_prognostic_step(grid_parent, time_stepping, S), get_prognostic_step(prog_parent, time_stepping, S), scratch_memory, S)

    if model isa PrimitiveWet
        humid_grid = get_prognostic_step(vars.grid.humidity, time_stepping, S)
        hole_filling!(humid_grid, model.hole_filling, model)  # remove negative humidity
    end

    # Batched spec→grid for the velocities: the general-purpose `:spectral_scratch` fuse packs
    # `(:a, :b)` (here holding U, V) into one Spectral3D parent, and `:uv_grid` packs `(:u, :v)`
    # TODO: theoretically we could merge this with the other big transform and then unscale coslat
    # seperately, but the dimensions don't quite align, shall we still do that in a hacky way?    
    transform!(get_prognostic_step(parent(vars.fused.uv_grid), time_stepping, S),
               parent(vars.fused.spectral_scratch), scratch_memory, S;
               unscale_coslat = true)

    # at initial step copy 2nd step (current) to 1st (prev) to retain those fields
    # only do after transforms to avoid copying uninitialized zeros
    initialize && move_prognostic_grid_variables_back!(vars, time_stepping, model)

    # convert the logarithm of surface pressure to actual surface pressure in Pascal for parameterizations
    # dispatch over DummyParameterization (= any parameterization) to let time steppers decide the step
    log_pₛ = get_prognostic_step(vars.grid.pressure, time_stepping, DummyParameterization())    # log Pa
    vars.parameterizations.surface_pressure .= exp.(log_pₛ)                                     # in Pa

    # include humidity effect into temp for everything stability-related
    temperature_average!(vars, temp, S)
    geopotential!(vars, model)                  # calculate geopotential
        
    for (name, tracer) in model.tracers
        tracer_var = get_prognostic_step(vars.prognostic.tracers[name], time_stepping, S)
        tracer_grid = get_prognostic_step(vars.grid.tracers[name], time_stepping, S)
        tracer.active && transform!(tracer_grid, tracer_var, scratch_memory, S)
    end

    # transform random pattern for random process unless random_process=nothing
    transform!(vars, model.random_process, S)

    return nothing
end

# only needed for Leapfrog to retain a copy of the current time step moved to the previous
# fallback is nothing for other time steppers
move_prognostic_grid_variables_back!(::Variables, ::AbstractTimeStepper, ::AbstractModel) = nothing