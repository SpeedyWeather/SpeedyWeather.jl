module SpeedyTransformsEnzymeExt

using Enzyme
using Enzyme.EnzymeCore
import .EnzymeRules: forward, reverse, augmented_primal
using .EnzymeRules

# import all functions for which we define rules
using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays

import SpeedyTransforms: _fourier!, wrapped_view

# NOTE: differentiate the gradient operators with `set_runtime_activity(Reverse)`. They destructure
# the coefficient arrays out of `S.gradients` and pass them as bare arrays into a kernel, so Enzyme
# computes ∂L/∂coefficients regardless of how `S` is annotated. Under STATIC activity analysis that
# gradient has nowhere to go and is written into the PRIMAL, silently corrupting every later call
# with the same transform. Measured on Julia 1.10 for `divergence!` — primal corrupted?
#
#     mode                   Const(S)   Duplicated(S, make_zero(S))
#     Reverse                yes        no
#     set_runtime_activity   no         no
#
# so runtime activity is what makes this correct; the argument annotation alone is not enough.
# `SpectralTransform` is also deliberately not declared `EnzymeRules.inactive_type`: that suppressed
# the shadow without suppressing the accumulation, so it corrupted the primal under both annotations.

# Rules for SpeedyTransforms

# _fourier!

# Computes the scale for the adjoint/pullback of all discrete Fourier transforms.
function adjoint_scale(S::SpectralTransform)
    (; nlons) = S
    (; nlat_half) = S.grid
    # The serial (K=1) plan vector (always built) is sufficient to read each ring's nfreq.
    rfft_plans_1D = S.rfft_plan_serial
    nfreqs = [rfft_plan.osz[1] for rfft_plan in rfft_plans_1D]

    scale = zeros(Int, maximum(nfreqs), 1, nlat_half) # the scratch memory is (Freq x lvl x lat), so we insert
    # an additional dimension here for easier matrix multiply

    for i in 1:nlat_half
        scale[1:nfreqs[i], 1, i] = rfft_adjoint_scale(nfreqs[i], nlons[i])
    end

    # TODO: transfer array to GPU in case we are on GPU?
    return scale
end

# Computes the scale for the adjoint/pullback of a real discrete fourier transform.
function rfft_adjoint_scale(n_freq::Int, n_real::Int)
    if iseven(n_real)
        return [1 < i < n_freq ? 2 : 1 for i in 1:n_freq]
    else
        return [1 < i ? 2 : 1 for i in 1:n_freq]
    end
end

### Custom rule for _fourier!(f_north, f_north, grid, S)
function augmented_primal(
        config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_fourier!)}, ::Type{<:Const},
        f_north::Duplicated, f_south::Duplicated, grids::Duplicated{<:AbstractField}, S::Union{Const, MixedDuplicated}
    )

    func.val(f_north.val, f_south.val, grids.val, S.val) # forward pass

    # save grids in tape if grids will be overwritten
    if overwritten(config)[4] # TODO: Not sure this is really necessary because grids won't ever get overwritten by this _fourier!
        tape = copy(grids.val)
    else
        tape = nothing
    end

    return AugmentedReturn(nothing, nothing, tape) # because the function actually returns nothing

end

function reverse(
        config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_fourier!)}, ::Type{<:Const}, tape,
        f_north::Duplicated, f_south::Duplicated, grids::Duplicated{<:AbstractField}, S::Union{Const, MixedDuplicated}
    )

    # adjoint/jvp of FFT has a different scaling, compute it, apply it later to f_north, f_south
    scale = adjoint_scale(S.val)

    # retrieve grids value, either from original grids or from tape if grids may have been overwritten.
    gridsval = overwritten(config)[4] ? tape : grids.val

    # compute the adjoint
    dgridval = zero(gridsval)
    _fourier!(dgridval, f_north.dval ./ scale, f_south.dval ./ scale, S.val) # inverse FFT (w/o normalization)
    grids.dval .+= dgridval

    # no derivative wrt the f_north and f_south that were input because they are overwritten.
    # `fill!` and not `make_zero!`: because make_zero! would zero the whole fused buffer
    fill!(f_north.dval, 0)
    fill!(f_south.dval, 0)

    # the function has no return values, so we also return nothing here
    return (nothing, nothing, nothing, nothing)
end

### Custom rule for _fourier!(grid, f_north, f_south, S)
function augmented_primal(
        config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_fourier!)}, ::Type{<:Const},
        grids::Duplicated{<:AbstractField}, f_north::Duplicated, f_south::Duplicated, S::Union{Const, MixedDuplicated}
    )

    func.val(grids.val, f_north.val, f_south.val, S.val) # forward pass

    # TODO: make an overwritten check here?

    return AugmentedReturn(nothing, nothing, nothing) # because the function actually returns nothing

end

function reverse(
        config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_fourier!)}, ::Type{<:Const}, tape,
        grids::Duplicated{<:AbstractField}, f_north::Duplicated, f_south::Duplicated, S::Union{Const, MixedDuplicated}
    )

    # adjoint/vjp of FFT has a different scaling, compute it, apply it later to f_north, f_south
    scale = adjoint_scale(S.val)

    # TODO: retrieve from tape here if overwritten?

    # compute the adjoint # TODO: could we reuse the f_north.val for that a well? and not allocate here
    dfnorthval = zero(f_north.val)
    dfsouthval = zero(f_south.val)

    _fourier!(dfnorthval, dfsouthval, grids.dval, S.val) # inverse FFT (w/o normalization)

    f_north.dval .+= scale .* dfnorthval
    f_south.dval .+= scale .* dfsouthval

    # no derivative wrt the grids that were input because they are overwritten.
    # `fill!` on the view's data, not `make_zero!` on the whole parent — see above.
    fill!(grids.dval.data, 0)

    # the function has no return values, so we also return nothing here
    return (nothing, nothing, nothing, nothing)
end

### FORWARD RULES
#
# Both the FFT and the full spectral transform are LINEAR in their input array, and `S` is fixed
# geometry carrying no tangent. The forward-mode tangent of a linear map is therefore the map
# itself applied to the tangent: run the primal on `.val`, then the identical call on each tangent
# `.dval`.

# `nlayer`-th tangent of an annotated argument, `nothing` for inactive (Const) arguments.
# Homogeneous-tuple indexing keeps this type stable for width > 1.
@inline _dval(x::Union{Duplicated, DuplicatedNoNeed}, ::Int) = x.dval
@inline _dval(x::Union{BatchDuplicated, BatchDuplicatedNoNeed}, b::Int) = x.dval[b]
@inline _dval(::Const, ::Int) = nothing

# Assemble what a forward rule must return for the mutated output argument `out`,
# following `EnzymeRules.forward_rule_return_type(config, RT)`.
@inline function _forward_return(config, out::Annotation)
    return if needs_shadow(config)
        if needs_primal(config)
            width(config) == 1 ? Duplicated(out.val, out.dval) : BatchDuplicated(out.val, out.dval)
        else
            out.dval
        end
    else
        needs_primal(config) ? out.val : nothing
    end
end

# One tangent pass of a linear in-place `op(dout, din, ...)`. An inactive (Const) input carries a
# zero tangent, so the output tangent is zeroed rather than transformed; an inactive output means
# Enzyme asserted the derivative is not needed, so there is nothing to propagate.
@inline function _linear_tangent!(op::F, dout, din, args...) where {F}
    dout === nothing && return nothing              # output inactive, nothing to propagate
    din === nothing && return fill!(dout.data, 0)   # input inactive ⇒ zero tangent (view-safe)
    op(dout, din, args...)
    return nothing
end

### Custom forward rule for _fourier!(f_north, f_south, field, S) — grid to Fourier space
function forward(
        config::EnzymeRules.FwdConfig, func::Const{typeof(_fourier!)}, ::Type{RT},
        f_north::Annotation{<:AbstractArray{<:Complex, 3}}, f_south::Annotation{<:AbstractArray{<:Complex, 3}},
        field::Annotation{<:AbstractField}, S::Annotation,
    ) where {RT <: Annotation}

    func.val(f_north.val, f_south.val, field.val, S.val)                    # primal
    for b in 1:width(config)                                                # one FFT per tangent
        dfield = _dval(field, b)
        dfnorth, dfsouth = _dval(f_north, b), _dval(f_south, b)
        if dfield === nothing                                               # inactive input ⇒ zero tangents
            dfnorth === nothing || fill!(dfnorth, 0)    # view-safe, see the reverse rules
            dfsouth === nothing || fill!(dfsouth, 0)
        elseif !(dfnorth === nothing || dfsouth === nothing)
            func.val(dfnorth, dfsouth, dfield, S.val)
        end
    end

    return nothing      # `_fourier!` returns nothing, so RT is Const and no shadow is requested
end

### Custom forward rule for _fourier!(field, f_north, f_south, S) — Fourier space to grid
function forward(
        config::EnzymeRules.FwdConfig, func::Const{typeof(_fourier!)}, ::Type{RT},
        field::Annotation{<:AbstractField}, f_north::Annotation{<:AbstractArray{<:Complex, 3}},
        f_south::Annotation{<:AbstractArray{<:Complex, 3}}, S::Annotation,
    ) where {RT <: Annotation}

    func.val(field.val, f_north.val, f_south.val, S.val)                    # primal
    for b in 1:width(config)                                                # one inverse FFT per tangent
        dfield = _dval(field, b)
        dfnorth, dfsouth = _dval(f_north, b), _dval(f_south, b)
        if dfnorth === nothing || dfsouth === nothing                       # inactive input ⇒ zero tangent
            dfield === nothing || fill!(dfield.data, 0)   # view-safe
        elseif dfield !== nothing
            func.val(dfield, dfnorth, dfsouth, S.val)
        end
    end

    return nothing      # `_fourier!` returns nothing, so RT is Const and no shadow is requested
end


### Analytic-adjoint rules for the whole `transform!` (positional cores `_transform_grid!` /
### `_transform_spec!`, which every `transform!` call routes through — chunked and batched alike).
#
# These rules make `transform!` a single AD boundary: the forward pass runs the primal transform
# unchanged (chunked or batched, allocation-free); the reverse pass applies the analytic adjoint
# of the (linear) spectral transform directly. This is preferred over letting Enzyme differentiate
# the internals because:
#  * differentiating the chunk loop is unsafe — Enzyme mis-constructs the per-iteration view
#    shadows (degenerate (0, 1) shadow → OOB in the `_fourier!` reverse / GC corruption on Julia
#    1.10) and reuses the last iteration's shadow for every chunk (silently zeroing all but the
#    last chunk's gradient);
#  * differentiating the batched path natively works but goes through `_legendre!` with an
#    S-derived loop bound loaded from the (mutable) model — the type-analysis failure that needs
#    `maxtypeoffset!` on Julia ≥ 1.11 — whereas an analytic rule sidesteps it and compiles cheaper;
#  * nested `autodiff`/`autodiff_deferred` inside a rule is not an option (Enzyme compilation
#    reentrancy → stack overflow).
# Note the rules treat `S` as inactive (`Const`-like, only `.val` read): there is no gradient
# w.r.t. the transform geometry itself (Legendre polynomials, quadrature weights). State AD and
# parameter AD for physical parameters are unaffected (those do not flow through `S`)
#
# Adjoint derivation (transform is linear; see docs/src/spectral_transform.md and legendre.jl):
#   synthesis (spec→grid) = inverse Legendre ∘ inverse FFT
#   analysis  (grid→spec) = forward FFT ∘ forward Legendre (the latter carries the solid-angle
#                           quadrature weight ΔΩ = sinθ Δθ Δϕ that the synthesis lacks).
#   The FFT adjoint reuses `adjoint_scale` exactly as the `_fourier!` reverse rules do; the
#   Legendre adjoint is the opposite-direction Legendre with the ΔΩ weight removed/added.
# Both pullbacks are FD-validated in the unit tests with EnzymeTestUtils: rel err ~1e-6, batched and chunked.

import SpeedyTransforms:
    _transform_grid!, _transform_spec!, _largest_planned_batch, _legendre!

# adjoint of spec→grid transform! w.r.t coeffs: field_bar → coeffs_bar (accumulates into coeffs_bar).
# Allocation-free apart from the FFT plan outputs (inherent to the primal too): the freq-space
# intermediates reuse the passed `scratch` (.north/.south/.column, which the forward pass no longer
# needs by the time the reverse runs), and the forward Legendre accumulates straight into
# coeffs_bar (`add=true`) instead of into a temporary spectral array.
function spec2grid_pullback!(coeffs_bar, field_bar, scratch, S; unscale_coslat::Bool = false)
    (; nlat_half) = S.grid
    K = size(field_bar, 2)
    K_batched = _largest_planned_batch(K, S)
    scale = adjoint_scale(S)
    dOmega = reshape(view(S.solid_angles, 1:nlat_half), 1, 1, :)
    clat = reshape(view(S.coslat⁻¹, 1:nlat_half), 1, 1, :)
    c = 1
    while c <= K
        c_end = min(c + K_batched - 1, K)
        chunk = c:c_end
        Kc = c_end - c + 1
        dg_n = view(scratch.north, :, 1:Kc, :)                          # reuse scratch, no allocation
        dg_s = view(scratch.south, :, 1:Kc, :)
        _fourier!(dg_n, dg_s, wrapped_view(field_bar, :, chunk), S)      # adjoint of inverse FFT: fwd FFT
        dg_n .*= scale                                                  # (also zeros stale padding rows: scale=0 there)
        dg_s .*= scale
        if unscale_coslat                                               # adjoint of the coslat unscaling
            dg_n .*= clat
            dg_s .*= clat
        end
        dg_n ./= dOmega                                                 # cancel the ΔΩ the fwd Legendre applies
        dg_s ./= dOmega
        # adjoint of inverse Legendre = fwd Legendre; accumulate onto the coeffs cotangent (add=true)
        _legendre!(wrapped_view(coeffs_bar, :, chunk), dg_n, dg_s, scratch.column, S; add = true)
        c = c_end + 1
    end
    return coeffs_bar
end

# adjoint of grid→spec transform! w.r.t field: coeffs_bar → field_bar (accumulates into field_bar).
# Allocation-free apart from the FFT plan outputs: freq-space intermediates reuse `scratch`, and the
# inverse FFT accumulates straight into field_bar (`add=true`) instead of into a temporary field.
function grid2spec_pullback!(field_bar, coeffs_bar, scratch, S)
    (; nlat_half) = S.grid
    K = size(coeffs_bar, 2)
    K_batched = _largest_planned_batch(K, S)
    scale = adjoint_scale(S)
    dOmega = reshape(view(S.solid_angles, 1:nlat_half), 1, 1, :)
    c = 1
    while c <= K
        c_end = min(c + K_batched - 1, K)
        chunk = c:c_end
        Kc = c_end - c + 1
        df_n = view(scratch.north, :, 1:Kc, :)                          # reuse scratch, no allocation
        df_s = view(scratch.south, :, 1:Kc, :)
        _legendre!(df_n, df_s, wrapped_view(coeffs_bar, :, chunk), scratch.column, S)   # adjoint of fwd Legendre: inv Legendre
        df_n .*= dOmega                                                 # re-apply the ΔΩ weight
        df_s .*= dOmega
        df_n ./= scale                                                  # adjoint of fwd FFT: inv FFT of df/scale
        df_s ./= scale
        # accumulate onto the field cotangent (add=true) instead of into a temporary field
        _fourier!(wrapped_view(field_bar, :, chunk), df_n, df_s, S; add = true)
        c = c_end + 1
    end
    return field_bar
end

function augmented_primal(
        config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_transform_grid!)}, ::Type{<:Annotation},
        field::Duplicated, coeffs::Duplicated, scratch::Union{Const, Duplicated},
        S::Union{Const, Duplicated, MixedDuplicated}, unscale_coslat::Const,
    )
    func.val(field.val, coeffs.val, scratch.val, S.val, unscale_coslat.val)
    primal = needs_primal(config) ? field.val : nothing
    shadow = needs_shadow(config) ? field.dval : nothing
    return AugmentedReturn(primal, shadow, nothing)
end

function reverse(
        ::EnzymeRules.RevConfigWidth{1}, ::Const{typeof(_transform_grid!)}, ::Type{<:Annotation}, tape,
        field::Duplicated, coeffs::Duplicated, scratch::Union{Const, Duplicated},
        S::Union{Const, Duplicated, MixedDuplicated}, unscale_coslat::Const,
    )
    spec2grid_pullback!(coeffs.dval, field.dval, scratch.val, S.val; unscale_coslat = unscale_coslat.val)
    # Zero the output cotangent: the primal OVERWRITES this argument, so its input-side cotangent is
    # zero (pinned against finite differences by the `test_reverse` checks in the unit tests).
    # `fill!(x.data, 0)` and NOT `make_zero!(x)`: these arguments are usually views into a fused
    # parent buffer shared with other variables, and `make_zero!` zeroes the whole parent, wiping
    # cotangents that belong to its neighbours. That silently made any reverse gradient seeded on a
    # grid-space variable (e.g. `vars.grid.u`) come back identically zero.
    fill!(field.dval.data, 0)
    return (nothing, nothing, nothing, nothing, nothing)
end

function augmented_primal(
        config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_transform_spec!)}, ::Type{<:Annotation},
        coeffs::Duplicated, field::Duplicated, scratch::Union{Const, Duplicated},
        S::Union{Const, Duplicated, MixedDuplicated},
    )
    func.val(coeffs.val, field.val, scratch.val, S.val)
    primal = needs_primal(config) ? coeffs.val : nothing
    shadow = needs_shadow(config) ? coeffs.dval : nothing
    return AugmentedReturn(primal, shadow, nothing)
end

function reverse(
        ::EnzymeRules.RevConfigWidth{1}, ::Const{typeof(_transform_spec!)}, ::Type{<:Annotation}, tape,
        coeffs::Duplicated, field::Duplicated, scratch::Union{Const, Duplicated},
        S::Union{Const, Duplicated, MixedDuplicated},
    )
    grid2spec_pullback!(field.dval, coeffs.dval, scratch.val, S.val)
    fill!(coeffs.dval.data, 0)
    return (nothing, nothing, nothing, nothing)
end

### Forward rules for the whole `transform!` (same positional cores as the reverse rules above).
# The transform is linear and `S` is inactive, so the tangent is the transform of the tangent.
# Having these makes `transform!` an AD boundary in forward mode as well, so Enzyme never differentiates
# the chunk loop or `_legendre!` internals (both problematic, see the comment above the reverse rules).
# `scratch` is write-before-read workspace, so every pass can share `scratch.val` and no shadow
# scratch is needed. The tangent passes run BEFORE the primal so that `scratch` is left holding the
# primal's content on return, making the rule's observable side effects identical to the primal's
# (the transform is linear, so no tangent pass depends on the primal's output).

function forward(
        config::EnzymeRules.FwdConfig, func::Const{typeof(_transform_grid!)}, ::Type{RT},
        field::Annotation{<:AbstractField}, coeffs::Annotation{<:LowerTriangularArray},
        scratch::Annotation, S::Annotation, unscale_coslat::Const,
    ) where {RT <: Annotation}

    for b in 1:width(config)                                                     # one transform per tangent
        _linear_tangent!(func.val, _dval(field, b), _dval(coeffs, b), scratch.val, S.val, unscale_coslat.val)
    end
    func.val(field.val, coeffs.val, scratch.val, S.val, unscale_coslat.val)      # primal last, see above

    return _forward_return(config, field)
end

function forward(
        config::EnzymeRules.FwdConfig, func::Const{typeof(_transform_spec!)}, ::Type{RT},
        coeffs::Annotation{<:LowerTriangularArray}, field::Annotation{<:AbstractField},
        scratch::Annotation, S::Annotation,
    ) where {RT <: Annotation}

    for b in 1:width(config)                                                     # one transform per tangent
        _linear_tangent!(func.val, _dval(coeffs, b), _dval(field, b), scratch.val, S.val)
    end
    func.val(coeffs.val, field.val, scratch.val, S.val)                          # primal last, see above

    return _forward_return(config, coeffs)
end

end
