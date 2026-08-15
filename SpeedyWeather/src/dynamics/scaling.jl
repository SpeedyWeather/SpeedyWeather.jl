"""$(TYPEDSIGNATURES)
Scales the prognostic variables vorticity and divergence with
the Earth's radius which is used in the dynamical core."""
@propagate_inbounds function scale_prognostic!(vars::Variables, scale::Real)
    progn = vars.prognostic             # for convenience
    new_scale = scale / progn.scale[]   # undo previous scale and new scale in one go
    haskey(progn, :vorticity) && (progn.vorticity .*= new_scale)
    haskey(progn, :divergence) && (progn.divergence .*= new_scale)
    # no need to actually scale the diagnostic variables as they will be
    # overwritten by the transform of the prognostic variables anyway
    progn.scale[] = scale               # store scaling information
    return vars
end

function scale_tendencies!(vars::Variables, model::AbstractModel)
    scale = vars.prognostic.scale[]
    (; tendencies) = vars
    TS = model.time_stepping

    # Scale each fused tendency parent exactly ONCE, then the standalone tendencies below.
    # Fused tendencies (e.g. the grid u/v/temperature/humidity/pressure tendencies) are
    # `SubArray` views into a shared fused parent buffer. Scaling them per-view would issue
    # several in-place broadcasts into the same buffer, which Reactant mis-handles within a
    # single compiled trace: multiple in-place updates to distinct views of one buffer corrupt
    # the data (this broke the reactant correctness tests). Scaling the whole parent once (a
    # single broadcast) is equivalent — the scale factor is uniform across members — and
    # Reactant-safe. Non-tendency intermediates that share these parents (uT_anomaly, uq,
    # kinetic_energy, the spectral u/v tendencies, …) are overwritten later in
    # `grid_tendencies!`/`spectral_tendencies!` before being read, so scaling them is harmless.
    haskey(vars.fused, :spectral_tendencies) &&
        scale!(get_tendency_step(parent(vars.fused.spectral_tendencies), TS, DummyParameterization()), scale)
    haskey(vars.fused, :grid_tendencies) &&
        scale!(get_tendency_step(parent(vars.fused.grid_tendencies), TS, DummyParameterization()), scale)

    # Scale the standalone (non-fused) tendencies individually — unrolled per name
    _scale_tendencies_unrolled!(vars, TS, scale)
    return nothing
end

# scale one standalone tendency; the fused-view members are skipped (their fuse parent is
# scaled contiguously above). With a concrete `var` the `is_view_entry` check constant-folds.
@inline function _scale_one_tendency!(var, TS, scale)
    is_view_entry(var) || scale!(get_tendency_step(var, TS, DummyParameterization()), scale)
    return nothing
end

@generated function _scale_tendencies_unrolled!(vars::Variables{Po, G, T}, TS, scale) where {Po, G, T}
    calls = Expr[]
    for name in _tendency_names(T)                  # spectral
        push!(calls, :(_scale_one_tendency!(getfield(vars.tendencies, $(QuoteNode(name))), TS, scale)))
    end
    for name in _tendency_and_uv_names(T)           # grid (+ u, v)
        push!(calls, :(_scale_one_tendency!(getfield(vars.tendencies.grid, $(QuoteNode(name))), TS, scale)))
    end
    for name in _namespace_names(T, :tracers)       # grid tracers
        push!(calls, :(_scale_one_tendency!(getfield(vars.tendencies.grid_tracers, $(QuoteNode(name))), TS, scale)))
    end
    return Expr(:block, calls..., :(return nothing))
end

function unscale_tendencies!(vars::Variables)
    scale = vars.prognostic.scale[]
    (; tendencies) = vars

    # Mirror `scale_tendencies!`: unscale each fused tendency parent once, then the standalone
    # tendencies individually (skipping fused members, which are views into the parents). This
    # avoids multiple in-place broadcasts into a shared buffer, which Reactant mis-handles in a
    # single compiled trace — see [`scale_tendencies!`](@ref) for details.
    inv_scale = inv(scale)
    haskey(vars.fused, :spectral_tendencies) && (parent(vars.fused.spectral_tendencies).data .*= inv_scale)
    haskey(vars.fused, :grid_tendencies) && (parent(vars.fused.grid_tendencies).data .*= inv_scale)

    # spectral
    for varname in tendency_names(vars)
        var = getfield(tendencies, varname)
        is_view_entry(var) || unscale!(var, scale)
    end

    # grid
    for varname in tendency_and_uv_names(vars)
        var = getfield(tendencies.grid, varname)
        is_view_entry(var) || unscale!(var, scale)
    end

    # tracers
    for varname in tracer_tendency_names(vars)
        var = getfield(tendencies.grid_tracers, varname)
        is_view_entry(var) || unscale!(var, scale)
    end
    return nothing
end


"""$(TYPEDSIGNATURES)
Undo the radius-scaling of vorticity and divergence from `scale_prognostic!(vars, scale::Real)`."""
function unscale!(vars::Variables)
    progn = vars.prognostic             # for convenience
    inv_scale = inv(progn.scale[])
    haskey(progn, :vorticity) && (progn.vorticity .*= inv_scale)
    haskey(progn, :divergence) && (progn.divergence .*= inv_scale)

    # and the corresponding grid variables if they exist
    haskey(vars.grid, :vorticity) && (vars.grid.vorticity .*= inv_scale)
    haskey(vars.grid, :divergence) && (vars.grid.divergence .*= inv_scale)

    # also unscale tendencies
    unscale_tendencies!(vars)

    progn.scale[] = 1                   # set scale back to 1=unscaled
    return vars
end

"""$(TYPEDSIGNATURES)
Scale the variable `var` with scalar `scale`."""
@inline function scale!(
        variable::Union{LowerTriangularArray, AbstractField},
        scale::Real
    )
    variable.data .*= scale
    return variable
end

"""$(TYPEDSIGNATURES)
Undo the scaling of the variable `var` with scalar `scale`."""
@inline function unscale!(
        variable::Union{LowerTriangularArray, AbstractField},
        scale::Real
    )
    return scale!(variable, inv(scale))
end
