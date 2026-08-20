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
"""
function set!(
        vars::NamedTuple,
        geometry::Geometry;
        kwargs...
    )
    @nospecialize kwargs
    return set_variables!(vars, geometry, _set_pairs(kwargs), _set_options(kwargs))
end

# Options of `set!` that control *how* variables are set; everything else in the kwargs is a
# variable to set. Collected into a fully concrete `SetOptions` so that `set_variables!` has
# exactly one signature.
const SET_OPTIONS = (:step, :add, :spectral_transform, :coslat_scaling_included, :static_func, :namespace)

# `step = nothing` is allowed to skip `get_step` for variables without a step dimension
_set_step(step) = isnothing(step) ? nothing : Int(step)

"""$(TYPEDSIGNATURES) Options controlling how `set!` sets a variable, see `set!`."""
struct SetOptions
    step::Union{Nothing, Int}
    add::Bool
    spectral_transform::Union{Nothing, AbstractSpectralTransform}
    coslat_scaling_included::Bool
    static_func::Bool
    namespace::Union{Nothing, Symbol}
end

# Type-erase the variables to set into a runtime Vector of pairs, and the options into a
# concrete struct. This is deliberate: the `Base.Pairs` kwargs type encodes the keyword names
# AND the value types, so passing them on makes `set_variables!` and everything it calls
# (`transform!`, the Legendre/FFT kernels, ...) re-specialize for every distinct combination
# used at a call site. There are eight such combinations among the initial conditions alone,
# each costing seconds of inference. Note the options must be passed *positionally* as a
# concrete struct: splatting them back out as keywords would reintroduce one specialization
# per distinct option NamedTuple and defeat the erasure. `set!` is a setup-path function
# called a handful of times per simulation, so the dynamic dispatch this introduces is
# irrelevant next to the compile time it saves.
function _set_pairs(@nospecialize(kwargs))
    return Pair{Symbol, Any}[k => v for (k, v) in kwargs if !(k in SET_OPTIONS)]
end

function _set_options(@nospecialize(kwargs))
    return SetOptions(
        _set_step(get(kwargs, :step, 1)),
        get(kwargs, :add, false),
        get(kwargs, :spectral_transform, nothing),
        get(kwargs, :coslat_scaling_included, false),
        get(kwargs, :static_func, true),
        get(kwargs, :namespace, nothing),
    )
end

# guard the coupling between `SET_OPTIONS` and the `SetOptions` fields: adding an option to
# one but not the other would silently reinterpret it as a variable name to set
@assert Set(SET_OPTIONS) == Set(fieldnames(SetOptions)) "SET_OPTIONS and SetOptions fields must match"

# as above but defaulting the spectral transform to the model's when not given explicitly
function _set_options(@nospecialize(kwargs), spectral_transform::AbstractSpectralTransform)
    return SetOptions(
        _set_step(get(kwargs, :step, 1)),
        get(kwargs, :add, false),
        get(kwargs, :spectral_transform, spectral_transform),
        get(kwargs, :coslat_scaling_included, false),
        get(kwargs, :static_func, true),
        get(kwargs, :namespace, nothing),
    )
end

"""$(TYPEDSIGNATURES)
Core of `set!`: sets the variables given as `name => value` pairs in `vars_to_set`, with
`options` controlling how. Deliberately free of keyword-dependent specialization so that all
call sites share a single compiled method, see `_set_pairs`."""
function set_variables!(
        vars::NamedTuple,
        geometry::Geometry,
        vars_to_set::Vector{Pair{Symbol, Any}},
        options::SetOptions,
    )
    (; step, add, spectral_transform, coslat_scaling_included, static_func, namespace) = options
    varnames = map(first, vars_to_set)

    # special case for u,v setting vor, div
    has_u = :u in varnames
    has_v = :v in varnames
    if has_u && has_v
        (; vorticity, divergence) = vars
        u = vars_to_set[findfirst(p -> first(p) === :u, vars_to_set)].second
        v = vars_to_set[findfirst(p -> first(p) === :v, vars_to_set)].second
        set_vordiv!(
            get_step(vorticity, step), get_step(divergence, step), u, v,
            geometry, spectral_transform; add, coslat_scaling_included, static_func
        )
    elseif has_u || has_v
        @warn "Only one of `u` and `v` provided, but both are needed to set `vor` and `div`. Skipping."
    end

    # normal case, search for varname in prognostic variables
    for (varname, value) in vars_to_set
        if varname in (:u, :v)  # already handled in special case above
            nothing
        elseif varname in keys(vars)
            var = ArrayDimensions.hastime(vars[varname]) ? get_step(vars[varname], step) : vars[varname]
            set!(var, value, geometry, spectral_transform; add, static_func)
        elseif namespace in keys(vars)
            if varname in keys(vars[namespace])
                var = ArrayDimensions.hastime(vars[namespace][varname]) ?
                    get_step(vars[namespace][varname], step) : vars[namespace][varname]
                set!(var, value, geometry, spectral_transform; add, static_func)
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

    if size(vor) != size(u_) != size(v_)
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
        kwargs...
    )
    vor .= 0    # curl of a constant is zero
    div .= 0    # divergence of a constant is zero
    return vor, div
end

"""$(TYPEDSIGNATURES)
Sets properties of the simuluation `S`. Convenience wrapper to call the other concrete 
`set!` methods. All `kwargs` are forwarded to these methods, which are documented 
seperately. See their documentation for possible `kwargs`. """
function set!(S::AbstractSimulation; group::Symbol = :prognostic, kwargs...)
    @nospecialize kwargs
    options = _set_options(kwargs, S.model.spectral_transform)
    return set_variables!(getfield(S.variables, group), S.model.geometry, _set_pairs(kwargs), options)
end

function set!(vars::Variables, model::AbstractModel; group::Symbol = :prognostic, kwargs...)
    @nospecialize kwargs
    vars.prognostic.scale[] != 1 && @warn "Prognostic variables are scaled with $(vars.prognostic.scale[]), but `set!` assumes unscaled variables."
    options = _set_options(kwargs, model.spectral_transform)
    return set_variables!(getfield(vars, group), model.geometry, _set_pairs(kwargs), options)
end
