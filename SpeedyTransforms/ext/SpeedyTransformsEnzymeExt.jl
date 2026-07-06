module SpeedyTransformsEnzymeExt

using Enzyme
using Enzyme.EnzymeCore
import .EnzymeRules: reverse, augmented_primal
using .EnzymeRules

# import all functions for which we define rules
using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays

import SpeedyTransforms: _fourier!, wrapped_view

# The spectral transform `S` is fixed geometry (Legendre polynomials, quadrature weights, FFT
# plans) and is never a differentiation target. Marking it inactive stops Enzyme from building a
# shadow of the large `SpectralTransform` aggregate when it is loaded out of a `Duplicated`
# (mutable) model — the "cannot deduce type of copy" type-analysis failure on Julia >= 1.11 that
# otherwise requires `Enzyme.API.maxtypeoffset!`. The transform! adjoint rules below read only
# `S.val`, so treating `S` as a constant everywhere is consistent.
EnzymeRules.inactive_type(::Type{<:SpectralTransform}) = true

# Rules for SpeedyTransforms

# _fourier!

# Computes the scale for the adjoint/pullback of all discrete Fourier transforms.
function adjoint_scale(S::SpectralTransform)
    (; nlons) = S
    (; nlat_half) = S.grid
    # The K=1 plan vector (always built) is sufficient to read each ring's nfreq.
    rfft_plans_1D = S.rfft_plans[1]
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

    # no derivative wrt the f_north and f_south that were input because they are overwritten
    make_zero!(f_north.dval)
    make_zero!(f_south.dval)

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

    # no derivative wrt the grids that were input because they are overwritten
    make_zero!(grids.dval)

    # the function has no return values, so we also return nothing here
    return (nothing, nothing, nothing, nothing)
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
# Both pullbacks are FD-validated (`_adjoint_check.jl`): rel err ~1e-6, batched and chunked.
# The pullbacks chunk over layers (Kc ≤ largest planned batch ≤ S.nlayers) so the internal
# `_fourier!`/`_legendre!` calls stay within the planned/serial FFT limits, matching the primal.

import SpeedyTransforms:
    _transform_grid!, _transform_spec!, _largest_planned_batch, _legendre!

# adjoint of spec→grid transform! w.r.t coeffs: field_bar → coeffs_bar (accumulates into coeffs_bar).
# Allocation-free apart from the FFT plan outputs (inherent to the primal too): the freq-space
# intermediates reuse the passed `scratch` (.north/.south/.column, which the forward pass no longer
# needs by the time the reverse runs), and the forward Legendre accumulates straight into
# coeffs_bar (`reset=false`) instead of into a temporary spectral array.
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
        # adjoint of inverse Legendre = fwd Legendre; accumulate onto the coeffs cotangent (reset=false)
        _legendre!(wrapped_view(coeffs_bar, :, chunk), dg_n, dg_s, scratch.column, S; reset = false)
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
    make_zero!(field.dval)      # the output cotangent has been propagated to coeffs
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
    make_zero!(coeffs.dval)     # the output cotangent has been propagated to field
    return (nothing, nothing, nothing, nothing)
end

end
