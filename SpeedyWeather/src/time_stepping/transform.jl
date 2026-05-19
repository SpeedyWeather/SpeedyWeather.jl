"""$(TYPEDSIGNATURES)
Propagate the spectral state of the prognostic variables of `vars` to the
grid variables `vars` for the barotropic vorticity model.
Updates grid vorticity, spectral stream function and spectral and grid velocities u, v."""
function SpeedyTransforms.transform!(
        vars::Variables,
        lf::Integer,
        model::Barotropic;
        kwargs...
    )
    u_grid = vars.grid.u
    v_grid = vars.grid.v
    vor_grid = vars.grid.vorticity

    # U = u*coslat, V=v*coslat
    U = vars.scratch.a                          # reuse work arrays for velocities in spectral
    V = vars.scratch.b                          # reuse work arrays for velocities in spectral
    vor = get_step(vars.prognostic.vorticity, lf)     # relative vorticity at leapfrog step lf

    scratch_memory = vars.scratch.transform_memory
    S = model.spectral_transform
    transform!(vor_grid, vor, scratch_memory, S)    # get vorticity on grid from spectral vor

    # get spectral U, V from spectral vorticity via stream function Ψ
    # U = u*coslat = -coslat*∂Ψ/∂lat
    # V = v*coslat = ∂Ψ/∂lon, radius omitted in both cases
    UV_from_vor!(U, V, vor, S)

    # transform from U, V in spectral to u, v on grid (U, V = u, v*coslat)
    transform!(u_grid, U, scratch_memory, S, unscale_coslat = true)
    transform!(v_grid, V, scratch_memory, S, unscale_coslat = true)

    for (name, tracer) in model.tracers
        tracer_var = get_step(vars.prognostic.tracers[name], lf)  # tracer at leapfrog step lf
        tracer.active && transform!(vars.grid.tracers[name], tracer_var, scratch_memory, S)
    end

    # transform random pattern for random process unless random_process=nothing
    transform!(vars, lf, model.random_process, S)

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
        lf::Integer,
        model::ShallowWater;
        kwargs...
    )

    vor_grid = vars.grid.vorticity
    u_grid = vars.grid.u
    v_grid = vars.grid.v
    div_grid = vars.grid.divergence
    η_grid = vars.grid.η

    vor = get_step(vars.prognostic.vorticity, lf)     # relative vorticity at leapfrog step lf
    div = get_step(vars.prognostic.divergence, lf)     # divergence at leapfrog step lf
    η = get_step(vars.prognostic.η, lf)         # interface displacement η at leapfrog step lf

    # U = u*coslat, V=v*coslat
    U = vars.scratch.a
    V = vars.scratch.b

    scratch_memory = vars.scratch.transform_memory
    S = model.spectral_transform

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
        tracer_var = get_step(vars.prognostic.tracers[name], lf)  # tracer at leapfrog step lf
        tracer.active && transform!(vars.grid.tracers[name], tracer_var, scratch_memory, S)
    end

    # transform random pattern for random process unless random_process=nothing
    transform!(vars, lf, model.random_process, S)

    return nothing
end

"""$(TYPEDSIGNATURES)
Propagate the spectral state of the prognostic variables of `vars` to the
grid variables in `vars` for primitive equation models. Updates grid vorticity,
grid divergence, grid temperature, pressure (`pres_grid`) and the velocities
u, v."""
function SpeedyTransforms.transform!(
        vars::Variables,
        lf::Integer,
        model::PrimitiveEquation;
        initialize::Bool = false,
    )

    pres_grid = vars.grid.pressure
    u_grid = vars.grid.u
    v_grid = vars.grid.v
    temp_grid = vars.grid.temperature

    pres_grid_prev = vars.grid.pressure_prev
    u_grid_prev = vars.grid.u_prev
    v_grid_prev = vars.grid.v_prev
    temp_grid_prev = vars.grid.temperature_prev

    vor = get_step(vars.prognostic.vorticity, lf)         # relative vorticity at leapfrog step lf
    div = get_step(vars.prognostic.divergence, lf)         # divergence at leapfrog step lf
    temp = get_step(vars.prognostic.temperature, lf)       # temperature at leapfrog step lf

    if model isa PrimitiveWet                       # dry model don't have humidity variables
        humid_grid = vars.grid.humidity
        humid_grid_prev = vars.grid.humidity_prev
    end

    scratch_memory = vars.scratch.transform_memory

    U = vars.scratch.a                              # reuse work arrays
    V = vars.scratch.b                              # U = u*coslat, V=v*coslat
    S = model.spectral_transform

    # retain previous time step for vertical advection and parameterizations
    if initialize == false                          # only store prev after initial step
        @. u_grid_prev = u_grid
        @. v_grid_prev = v_grid
        @. temp_grid_prev = temp_grid
        @. pres_grid_prev = exp(pres_grid)

        if model isa PrimitiveWet
            @. humid_grid_prev = humid_grid
        end

        for (name, tracer) in model.tracers
            if tracer.active
                name_prev = Symbol(name, :_prev)
                vars.grid.tracers[name_prev] .= vars.grid.tracers[name]
            end
        end
    end

    # Mega-batched spec→grid for the prognostic state: one call covers vorticity, divergence,
    # temperature, pressure (and humidity for PrimitiveWet). The fused parents pack their
    # members along axis 2 in matching declaration order, so slot k of `vars.fused.prognostic`
    # maps to slot k of `vars.fused.grid` (enforced by `_assert_fuse_alignment` at construction).
    prog_parent = parent(vars.fused.prognostic)
    grid_parent = parent(vars.fused.grid)
    transform!(grid_parent, get_step(prog_parent, lf), scratch_memory, S)

    if model isa PrimitiveWet
        hole_filling!(humid_grid, model.hole_filling, model)  # remove negative humidity
    end

    # get spectral U, V from vorticity and divergence via stream function Ψ and vel potential ϕ
    # U = u*coslat = -coslat*∂Ψ/∂lat + ∂ϕ/dlon
    # V = v*coslat =  coslat*∂ϕ/∂lat + ∂Ψ/dlon
    UV_from_vordiv!(U, V, vor, div, S)

    # Batched spec→grid for the velocities: the general-purpose `:spectral_scratch` fuse packs
    # `(:a, :b)` (here holding U, V) into one Spectral3D parent, and `:uv_grid` packs `(:u, :v)`
    # TODO: theoretically we could merge this with the other big transform and then unscale coslat
    # seperately, shall we do that?
    transform!(parent(vars.fused.uv_grid), parent(vars.fused.spectral_scratch), scratch_memory, S;
               unscale_coslat = true)

    # include humidity effect into temp for everything stability-related
    temperature_average!(vars, temp, S)
    geopotential!(vars, model)                  # calculate geopotential

    if initialize   # at initial step store prev <- current
        @. u_grid_prev = u_grid
        @. v_grid_prev = v_grid
        @. temp_grid_prev = temp_grid
        @. pres_grid_prev = exp(pres_grid)      # store pressure in Pa not log(Pa) for parameterizations

        if model isa PrimitiveWet
            @. humid_grid_prev = humid_grid
        end

        for (name, tracer) in model.tracers
            if tracer.active
                name_prev = Symbol(name, :_prev)
                vars.grid.tracers[name_prev] .= vars.grid.tracers[name]
            end
        end
    end

    for (name, tracer) in model.tracers
        tracer_var = get_step(vars.prognostic.tracers[name], lf)  # tracer at leapfrog step lf
        tracer.active && transform!(vars.grid.tracers[name], tracer_var, scratch_memory, S)
    end

    # transform random pattern for random process unless random_process=nothing
    transform!(vars, lf, model.random_process, S)

    return nothing
end
