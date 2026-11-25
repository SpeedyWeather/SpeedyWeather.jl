export set!

"""
$(TYPEDSIGNATURES)
Sets new values for variables defined through keyword matching the keys in a NamedTuple.
Spectral variables can be set at timestep index `step`. If `add==true` they are added to the 
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

Specify the namespace as a symbol in case the `vars::NamedTuple` contains them, e.g.

    set!(vars, sea_surface_temperature = 1, namespace=:ocean)
   
The keyword argument `step` specifies the step index to set, default is 1.
For variables without a step dimension use `step = nothing` to skip the `get_step`.

With `truncate = true` (default) spectral variables are truncated after they are set, removing the
extra degree ``l_{max}+1`` that all spectral variables carry but only vector quantities use
(see [`truncate_extra_degree!`](@ref)). All prognostic variables are scalars, so this degree is
otherwise filled by the transform from grid space (or by copying from a larger spectral field)
but never changed by the time stepping. Use `truncate = false` to set spectral vector quantities,
e.g. the spectral `u_tendency` in `group = :dynamics`.
"""
function set!(
        vars::NamedTuple,
        geometry::Geometry;
        step = 1,
        add::Bool = false,
        spectral_transform::Union{Nothing, AbstractSpectralTransform} = nothing,
        coslat_scaling_included::Bool = false,
        static_func::Bool = true,
        namespace::Union{Nothing, Symbol} = nothing,
        scratch_field::Union{Nothing, AbstractField} = nothing,
        truncate::Bool = true,
        kwargs...
    )
    # special case for u,v setting vor, div
    if :u in keys(kwargs) && :v in keys(kwargs)
        vor = get_step(vars.vorticity, step)
        div = get_step(vars.divergence, step)
        set_vordiv!(
            vor, div, kwargs[:u], kwargs[:v],
            geometry, spectral_transform; add, coslat_scaling_included, static_func
        )
        truncate && truncate_extra_degree!(vor)
        truncate && truncate_extra_degree!(div)
    elseif :u in keys(kwargs) || :v in keys(kwargs)
        @warn "Only one of `u` and `v` provided, but both are needed to set `vor` and `div`. Skipping."
    end

    # normal case, search for varname in prognostic variables
    for varname in keys(kwargs)
        if varname in (:u, :v)  # already handled in special case above
            nothing
        elseif varname in keys(vars)
            var = ArrayDimensions.hastime(vars[varname]) ? get_step(vars[varname], step) : vars[varname]
            set!(var, kwargs[varname], geometry, spectral_transform; add, static_func, scratch_field)
            truncate && truncate_extra_degree!(var)
        elseif namespace in keys(vars)
            if varname in keys(vars[namespace])
                var = ArrayDimensions.hastime(vars[namespace][varname]) ?
                    get_step(vars[namespace][varname], step) : vars[namespace][varname]
                set!(var, kwargs[varname], geometry, spectral_transform; add, static_func, scratch_field)
                truncate && truncate_extra_degree!(var)
            else
                # throw error if varname can't be found and print existing variables
                @warn "`$varname` not defined in NamedTuple with keys = $(keys(vars[namespace])). Skipping."
            end
        else
            # throw error if varname can't be found and print existing variables
            @warn "`$varname` not defined in NamedTuple with keys = $(keys(vars)). Skipping."
        end
    end
    return nothing
end

"""$(TYPEDSIGNATURES) set the `variables`, default is in the group `prognostic`."""
set!(variables::Variables, args...; group::Symbol = :prognostic, kwargs...) =
    set!(getfield(variables, group), args...; kwargs...)

export truncate_extra_degree!

"""$(TYPEDSIGNATURES)
Set the extra degree of the spectral variable(s) `vars` to zero, i.e. all harmonics with a degree
``l`` larger than the truncation (``l > m_{max}``, a no-op for square `LowerTriangularArray`s).
All spectral variables carry one more degree ``l_{max} = m_{max}+1`` (1-based) as vector quantities
(e.g. ``u, v`` or gradients) need it for the meridional derivative. Scalar quantities, which all
prognostic variables are, must not use it: the dynamical core never changes this degree (its
tendencies are zero there), but it would still enter the grid-point fields in every time step.
Recurses into (namespaced) `NamedTuple`s, e.g. `vars.prognostic`, and leaves non-spectral
variables untouched."""
truncate_extra_degree!(L::LowerTriangularArray) = SpeedyTransforms.spectral_truncation!(L)
truncate_extra_degree!(vars::NamedTuple) = (foreach(truncate_extra_degree!, values(vars)); vars)
truncate_extra_degree!(var) = var     # grid variables, clock, ...

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
        # convert to number format in S, needed for FFTW, but skipped if that's already the case
        # as broadcasting the conversion would allocate (and is not traceable by Reactant)
        field = eltype(field) == eltype(S) ? field : convert.(eltype(S), field)

        # only allocate temporary array if we need to transform to a different truncation
        if var.spectrum == S.spectrum && !add
            transform!(var, field, S)
            return nothing
        else
            specs = transform(field, S)
        end
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
        scratch_field::Union{AbstractField, Nothing} = nothing,
        kwargs...,
    )
    (; grid, nlayers, NF) = geometry.spectral_grid
    field = isnothing(scratch_field) ? (ndims(var) == 1 ? zeros(NF, grid) : zeros(NF, grid, nlayers)) : scratch_field
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

    # set the [1, :] row in a GPU-safe way via a view on the underlying data
    data = var_new.data
    norm_s = convert(eltype(data), norm_sphere * s)
    if ndims(data) == 1
        # 1D LTA: data is a Vector, only one matrix, so a single element
        view(data, 1:1) .= norm_s
    else
        view(data, 1, ntuple(_ -> :, ndims(data) - 1)...) .= norm_s
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
        kwargs...
    )
    (; londs, latds, σ_levels_full) = geometry

    # Reactant: use broadcasting (no KA kernel, no scalar indexing)
    if typeof(architecture(var)) <: ReactantDevice
        _set_function_3d_broadcast!(var, f, londs, latds, σ_levels_full; add)
        # on GPU no dynamically generated function are allowd in kernels, transfer them to CPU and back
    elseif typeof(architecture(var)) <: GPU && static_func == false
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
    launch!(
        architecture(var), RingGridWorkOrder, size(var), set_field_3d_kernel!,
        var, londs, latds, σ_levels_full, f, kernel_func
    )
    return var
end

@kernel function set_field_3d_kernel!(var, londs, latds, σ_levels_full, f, kernel_func)
    ij, k = @index(Global, NTuple)
    var[ij, k] = kernel_func(var[ij, k], f(londs[ij], latds[ij], σ_levels_full[k]))
end

# Reactant-compatible broadcasting version: no KA kernel, no scalar indexing
function _set_function_3d_broadcast!(
        var::AbstractField, f::Function,
        londs::AbstractVector, latds::AbstractVector, σ_levels_full::AbstractVector;
        add::Bool = false,
    )
    # reshape σ for broadcasting: londs/latds are (npoints,), σ_row is (1, nlayers)
    NF = eltype(var)
    σ_row = reshape(σ_levels_full, 1, :)
    val_raw = f.(londs, latds, σ_row)
    val = eltype(val_raw) === NF ? val_raw : NF.(val_raw)
    if add
        var.data .= var.data .+ val
    else
        var.data .= val
    end
    return var
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

# set vor_div <- func, evaluated on the grid to avoid a spectral round trip of u, v before
# the 1/coslat scaling (which can only be applied in grid space)
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
    (; grid, NF) = geometry.spectral_grid
    u = ndims(vor) == 1 ? zeros(NF, grid) : zeros(NF, grid, size(vor, 2))
    v = similar(u)
    set!(u, u_func, geometry, S; kwargs...)
    set!(v, v_func, geometry, S; kwargs...)

    return set_vordiv!(vor, div, u, v, geometry, S; add, coslat_scaling_included)
end

# set vor_div <- grid
function set_vordiv!(
        vor::LowerTriangularArray,
        div::LowerTriangularArray,
        u::AbstractField,
        v::AbstractField,
        geometry::Geometry,
        S::Union{Nothing, AbstractSpectralTransform} = nothing;
        add::Bool = false,
        coslat_scaling_included::Bool = false,
        kwargs...,
    )
    S = isnothing(S) ? SpectralTransform(geometry.spectral_grid) : S
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
        S::Union{Nothing, AbstractSpectralTransform} = nothing;
        add::Bool = false,
        coslat_scaling_included::Bool = false,
        kwargs...,
    )
    S = isnothing(S) ? SpectralTransform(geometry.spectral_grid) : S

    # first bring u, v to the size of vor, div (and S), copyto! copies the largest matching
    # subset of harmonics, i.e. truncates larger or zero-pads smaller u, v
    u_ = size(u) == size(vor) ? u : copyto!(zero(vor), u)
    v_ = size(v) == size(vor) ? v : copyto!(zero(vor), v)

    # the 1/coslat scaling can only be applied in grid space
    if !coslat_scaling_included
        u_ = transform(RingGrids.scale_coslat⁻¹(transform(u_, S)), S)
        v_ = transform(RingGrids.scale_coslat⁻¹(transform(v_, S)), S)
    end

    radius = geometry.radius[]
    if add  # not curl!(...; add, radius) as that also divides the existing vor, div by the radius
        vor .+= SpeedyTransforms.curl(u_, v_, S; radius)
        div .+= SpeedyTransforms.divergence(u_, v_, S; radius)
    else
        curl!(vor, u_, v_, S; radius)
        divergence!(div, u_, v_, S; radius)
    end
    return vor, div
end

# if number provided for u, v then vor = div = 0
function set_vordiv!(
        vor::LowerTriangularArray,
        div::LowerTriangularArray,
        u::Number,  # this is a constant field, curl and div are zero
        v::Number,  # this is a constant field, curl and div are zero
        geometry::Geometry,
        S::Union{Nothing, SpectralTransform} = nothing;
        add::Bool = false,
        kwargs...
    )
    add && return vor, div      # adding zero vorticity and divergence changes nothing
    vor .= 0    # curl of a constant is zero
    div .= 0    # divergence of a constant is zero
    return vor, div
end

"""$(TYPEDSIGNATURES)
Sets properties of the simuluation `S`. Convenience wrapper to call the other concrete 
`set!` methods. All `kwargs` are forwarded to these methods, which are documented 
seperately. See their documentation for possible `kwargs`. """
function set!(S::AbstractSimulation; kwargs...)
    return set!(S.variables, S.model.geometry; spectral_transform = S.model.spectral_transform, kwargs...)
end

function set!(vars::Variables, model::AbstractModel; kwargs...)
    vars.prognostic.scale[] != 1 && @warn "Prognostic variables are scaled with $(vars.prognostic.scale[]), but `set!` assumes unscaled variables."
    return set!(vars, model.geometry; spectral_transform = model.spectral_transform, kwargs...)
end
