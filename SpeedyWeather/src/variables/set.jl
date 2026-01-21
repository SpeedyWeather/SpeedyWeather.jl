export set!

"""
$(TYPEDSIGNATURES)
Sets new values for the keyword arguments (velocities, vorticity, divergence, etc..) into the
prognostic variable struct `progn` at timestep index `lf`. If `add==true` they are added to the 
current value instead. If a `AbstractSpectralTransform` S is provided, it is used when needed to set 
the variable, otherwise it is recomputed. In case `u` and `v` are provied, actually the divergence
and vorticity are set and `coslat_scaling_included` specficies whether or not the 1/cos(lat) 
scaling is already included in the arrays or not (default: `false`). If a function or callable 
object is provided, `static_func` specficies whether or not the function is static (i.e. does not 
contain any dynamic code) or not (default: `true`). On GPU, only static functions are executed 
efficiently.

The input may be:
* A function or callable object `f(lond, latd, σ) -> value` (multilevel variables) 
* A function or callable object `f(lond, latd) -> value` (surface level variables)
* An instance of `AbstractField` 
* An instance of `LowerTriangularArray` 
* A scalar `<: Number` (interpreted as a constant field in grid space)
"""
function set!(
        progn::PrognosticVariables,
        geometry::Geometry;
        u = nothing,
        v = nothing,
        vor = nothing,
        div = nothing,
        temp = nothing,
        humid = nothing,
        pres = nothing,
        sea_surface_temperature = nothing,
        sea_ice_concentration = nothing,
        soil_temperature = nothing,
        snow_depth = nothing,
        soil_moisture = nothing,
        lf::Integer = 1,
        add::Bool = false,
        spectral_transform::Union{Nothing, AbstractSpectralTransform} = nothing,
        coslat_scaling_included::Bool = false,
        static_func::Bool = true,
        kwargs...
    )
    # ATMOSPHERE
    isnothing(vor)   || set!(get_step(progn.vor, lf), vor, geometry, spectral_transform; add, static_func)
    isnothing(div)   || set!(get_step(progn.div, lf), div, geometry, spectral_transform; add, static_func)
    isnothing(temp)  || set!(get_step(progn.temp, lf), temp, geometry, spectral_transform; add, static_func)
    isnothing(humid) || set!(get_step(progn.humid, lf), humid, geometry, spectral_transform; add, static_func)
    isnothing(pres)  || set!(get_step(progn.pres, lf), pres, geometry, spectral_transform; add, static_func)

    # or provide u, v instead of vor, div
    isnothing(u) | isnothing(v) || set_vordiv!(get_step(progn.vor, lf), get_step(progn.div, lf), u, v, geometry, spectral_transform; add, coslat_scaling_included, static_func)

    # OCEAN
    isnothing(sea_surface_temperature)  || set!(progn.ocean.sea_surface_temperature, sea_surface_temperature, geometry, spectral_transform; add, static_func)
    isnothing(sea_ice_concentration)    || set!(progn.ocean.sea_ice_concentration, sea_ice_concentration, geometry, spectral_transform; add, static_func)

    # LAND
    isnothing(soil_temperature)         || set!(progn.land.soil_temperature, soil_temperature, geometry, spectral_transform; add, static_func)
    isnothing(snow_depth)               || set!(progn.land.snow_depth, snow_depth, geometry, spectral_transform; add, static_func)
    isnothing(soil_moisture)            || set!(progn.land.soil_moisture, soil_moisture, geometry, spectral_transform; add, static_func)

    # TRACERS
    for varname in keys(kwargs)
        if varname in keys(progn.tracers)
            tracer_var = get_step(progn.tracers[varname], lf)
            set!(tracer_var, kwargs[varname], geometry, spectral_transform; add, static_func)
        else
            throw(UndefVarError(varname))
        end
    end
    return
end


# set LTA <- LTA
function set!(
        var::LowerTriangularArray,
        L::LowerTriangularArray,
        varargs...;
        add::Bool = false,
        kwargs...,
    )
    if add
        if size(var) == size(L)
            var .+= L
        else
            L_var = SpeedyTransforms.spectral_truncation(L, size(var, 1, as = Matrix), size(var, 2, as = Matrix))
            var .+= L_var
        end
    else
        size(var) != size(L) || fill!(var, 0) # copyto! copies over the largest subset, when size(var) > size(L), the copyto! isn't enough by itself
        copyto!(var, L)
    end
    return var
end

# set LTA <- Grid
function set!(
        var::LowerTriangularArray,
        field::AbstractField,
        geometry::Union{Geometry, Nothing} = nothing,
        S::Union{Nothing, AbstractSpectralTransform} = nothing;
        add::Bool = false,
        kwargs...,
    )
    if isnothing(S)
        specs = transform(field)
    else
        # convert to number format in S, needed for FFTW
        field = convert.(eltype(S), field)
        specs = transform(field, S)
    end
    return set!(var, specs; add, kwargs...)
end

# set LTA <- func
function set!(
        var::LowerTriangularArray,
        f::Function,
        geometry::Geometry,
        S::Union{AbstractSpectralTransform, Nothing} = nothing;
        add::Bool = false,
        kwargs...,
    )
    (; grid, nlayers, NF) = geometry.spectral_grid
    field = ndims(var) == 1 ? zeros(NF, grid) : zeros(NF, grid, nlayers)
    set!(field, f, geometry, S; add = false, kwargs...)
    return set!(var, field, geometry, S; add, kwargs...)
end

# set LTA <- number
function set!(
        var::LowerTriangularArray,
        s::Number,
        geometry::Geometry,
        S::Union{AbstractSpectralTransform, Nothing} = nothing;
        add::Bool = false,
        kwargs...,
    )
    # appropiate normalization, assume standard 2√π normalisation if no transform is given
    norm_sphere = isnothing(S) ? 2sqrt(π) : S.norm_sphere

    # all elements are zero except for the 0,0 one
    var_new = zero(var)

    for k in eachmatrix(var_new)
        var_new[1, k] = norm_sphere * s
    end

    return set!(var, var_new, geometry, S; add, kwargs...)
end

# set Field <- Field
function set!(
        var::AbstractField,
        field::AbstractField,
        geometry::Geometry,
        S::Union{Nothing, AbstractSpectralTransform} = nothing;
        add::Bool = false,
        kwargs...,
    )
    if add
        if fields_match(var, field)
            var .+= field
        else
            var .+= interpolate(var.grid, field; NF = eltype(var))
        end
    else
        interpolate!(var, field; NF = eltype(var))
    end
    return var
end

# set Field <- LTA
function set!(
        var::AbstractField,
        specs::LowerTriangularArray,
        geometry::Geometry,
        S::Union{Nothing, AbstractSpectralTransform} = nothing;
        add::Bool = false,
        kwargs...,
    )
    field = isnothing(S) ? transform(specs) : transform(specs, S)
    return set!(var, field, geometry, S; add, kwargs...)
end

# set Field <- Func
function set!(
        var::AbstractField3D,
        f::Function,
        geometry::Geometry,
        S::Union{Nothing, AbstractSpectralTransform} = nothing;
        add::Bool = false,
        static_func = true,
    )
    (; londs, latds, σ_levels_full) = geometry

    # on GPU no dynamically generated function are allowd in kernels, transfer them to CPU and back
    if typeof(architecture(var)) <: GPU && static_func == false
        arch_cpu = CPU()

        var_cpu = on_architecture(arch_cpu, var)
        _set_function_3d!(var_cpu, f, adapt(Array, londs), adapt(Array, latds), adapt(Array, σ_levels_full); add = add)
        var.data .= on_architecture(architecture(var), var_cpu.data)
    else
        _set_function_3d!(var, f, londs, latds, σ_levels_full; add = add)
    end

    return var
end

function _set_function_3d!(var::AbstractField, f::Function, londs::AbstractVector, latds::AbstractVector, σ_levels_full::AbstractVector; add::Bool = false)
    kernel_func = add ? (a, b) -> a + b : (a, b) -> b

    @boundscheck size(var) == (length(londs), length(σ_levels_full)) || throw(DimensionMismatch())
    return launch!(
        architecture(var), RingGridWorkOrder, size(var), set_field_3d_kernel!,
        var, londs, latds, σ_levels_full, f, kernel_func
    )
end

@kernel function set_field_3d_kernel!(var, londs, latds, σ_levels_full, f, kernel_func)
    ij, k = @index(Global, NTuple)
    var[ij, k] = kernel_func(var[ij, k], f(londs[ij], latds[ij], σ_levels_full[k]))
end

# if geometry available
function set!(
        var::AbstractField2D,
        f::Function,
        geometry::Geometry,
        S::Union{Nothing, AbstractSpectralTransform} = nothing;
        kwargs...
    )

    (; londs, latds) = geometry     # use coordinates from geometry
    return _set!(var, f, londs, latds; kwargs...)
end

# otherwise recompute longitude, latitude vectors
function set!(
        var::AbstractField2D,
        f::Function,
        S::Union{Nothing, AbstractSpectralTransform} = nothing;
        kwargs...
    )
    # otherwise recompute longitude, latitude vectors
    londs, latds = RingGrids.get_londlatds(var)
    londs = on_architecture(architecture(var), londs)
    latds = on_architecture(architecture(var), latds)
    return _set!(var, f, londs, latds; kwargs...)
end

# set Grid (surface/single level) <- Func
function _set!(
        var::AbstractField2D,
        f::Function,
        londs::AbstractVector,
        latds::AbstractVector;
        add::Bool = false,
        kwargs...,
    )
    kernel_func = add ? (a, b) -> a + b : (a, b) -> b
    var.data .= kernel_func.(var.data, f.(londs, latds))
    return var
end

# set Grid <- Number
function set!(
        var::AbstractField,
        s::Number,
        geometry::Union{Geometry, Nothing} = nothing,
        S::Union{Nothing, AbstractSpectralTransform} = nothing;
        add::Bool = false,
        kwargs...,
    )
    kernel = add ? (a, b) -> a + b : (a, b) -> b
    s = convert(eltype(var), s)
    return var .= kernel.(var, s)
end

# set vor_div <- func
function set_vordiv!(
        vor::LowerTriangularArray,
        div::LowerTriangularArray,
        u_func,
        v_func,
        geometry::Geometry,
        S::Union{Nothing, AbstractSpectralTransform} = nothing;
        add::Bool = false,
        coslat_scaling_included::Bool = false,
        kwargs...,
    )
    u_L = similar(vor)
    set!(u_L, u_func, geometry, S; kwargs...)
    v_L = similar(vor)
    set!(v_L, v_func, geometry, S; kwargs...)

    return set_vordiv!(vor, div, u_L, v_L, geometry, S; add, coslat_scaling_included, kwargs...)
end

# set vor_div <- grid
function set_vordiv!(
        vor::LowerTriangularArray,
        div::LowerTriangularArray,
        u::AbstractField,
        v::AbstractField,
        geometry::Geometry,
        S::AbstractSpectralTransform = SpectralTransform(geometry.spectral_grid);
        add::Bool = false,
        coslat_scaling_included::Bool = false,
        kwargs...,
    )
    u_ = coslat_scaling_included ? u : RingGrids.scale_coslat⁻¹(u)
    v_ = coslat_scaling_included ? v : RingGrids.scale_coslat⁻¹(v)

    # convert to number format of spectral transform, otherwise FFTW complains
    u_ = eltype(S) == eltype(u_) ? u_ : convert.(eltype(S), u_)
    v_ = eltype(S) == eltype(v_) ? v_ : convert.(eltype(S), v_)

    u_spec = transform(u_, S)
    v_spec = transform(v_, S)

    return set_vordiv!(vor, div, u_spec, v_spec, geometry, S; add, coslat_scaling_included = true, kwargs...)
end

# set vor_div <- LTA
function set_vordiv!(
        vor::LowerTriangularArray,
        div::LowerTriangularArray,
        u::LowerTriangularArray,
        v::LowerTriangularArray,
        geometry::Geometry,
        S::AbstractSpectralTransform = SpectralTransform(geometry.spectral_grid);
        add::Bool = false,
        coslat_scaling_included::Bool = false,
        kwargs...,
    )
    u_ = coslat_scaling_included ? u : transform(RingGrids.scale_coslat⁻¹(transform(u, S)), S)
    v_ = coslat_scaling_included ? v : transform(RingGrids.scale_coslat⁻¹(transform(u, S)), S)
    radius = geometry.radius[]

    return if size(vor) != size(u_) != size(v_)
        u_new = zero(vor)
        copyto!(u_new, u_)

        v_new = zero(vor)
        copyto!(v_new, v_)

        curl!(vor, u_new, v_new, S; add, radius)
        divergence!(div, u_new, v_new, S; add, radius)
    else
        curl!(vor, u_, v_, S; add, radius)
        divergence!(div, u_, v_, S; add, radius)
    end
end

"""
$(TYPEDSIGNATURES)

Sets properties of the simuluation `S`. Convenience wrapper to call the other concrete 
`set!` methods. All `kwargs` are forwarded to these methods, which are documented 
seperately. See their documentation for possible `kwargs`. 
"""
function set!(S::AbstractSimulation; kwargs...)
    return set!(S.prognostic_variables, S.model.geometry; spectral_transform = S.model.spectral_transform, kwargs...)
end

function set!(progn::PrognosticVariables, model::AbstractModel; kwargs...)
    progn.scale[] != 1 && @warn "Prognostic variables are scaled with $(progn.scale[]), but `set!` assumes unscaled variables."
    return set!(progn, model.geometry; spectral_transform = model.spectral_transform, kwargs...)
end
