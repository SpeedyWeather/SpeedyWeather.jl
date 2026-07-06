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

### Custom rules for the CHUNKED transforms (_chunked_spec2grid!/_chunked_grid2spec!)
#
# The chunk loop passes `wrapped_view(field, :, chunk)` — a Field/LowerTriangularArray wrapping
# a SubArray — into `transform!`. Enzyme cannot be allowed to differentiate through that loop:
#  (a) it mis-constructs the shadow of the nested view wrapper when the parent is itself
#      view-backed (degenerate (0, 1) shadow → OOB in the `_fourier!` reverse; GC corruption
#      on Julia 1.10), and
#  (b) even for plain parents it reuses the LAST iteration's view shadows for ALL chunk
#      reverses, silently zeroing every chunk's gradient except the last.
# These rules replace the whole chunked transform: the forward pass runs the primal chunk loop
# unchanged (allocation-free views); the reverse pass applies the analytic adjoint of the
# (linear) spectral transform directly — no differentiation of the loop, no nested autodiff
# (nested `autodiff`/`autodiff_deferred` inside a rule triggers Enzyme compilation reentrancy /
# stack overflow). The batched (non-chunked) path is left to native Enzyme + the `_fourier!`
# rules above, which is correct and cheap there; these rules only cover the chunked path.
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
    _chunked_spec2grid!, _chunked_grid2spec!, _largest_planned_batch, _legendre!, ColumnScratchMemory

# adjoint of spec→grid transform! w.r.t coeffs: field_bar → coeffs_bar (accumulates into coeffs_bar)
function spec2grid_pullback!(coeffs_bar, field_bar, S; unscale_coslat::Bool = false)
    NF = eltype(S)
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
        dg_n = zeros(Complex{NF}, S.nfreq_max, Kc, nlat_half)
        dg_s = zeros(Complex{NF}, S.nfreq_max, Kc, nlat_half)
        _fourier!(dg_n, dg_s, wrapped_view(field_bar, :, chunk), S)      # adjoint of inverse FFT: fwd FFT
        dg_n .*= scale
        dg_s .*= scale
        if unscale_coslat                                               # adjoint of the coslat unscaling
            dg_n .*= clat
            dg_s .*= clat
        end
        dg_n ./= dOmega                                                 # cancel the ΔΩ the fwd Legendre applies
        dg_s ./= dOmega
        col = ColumnScratchMemory(zeros(Complex{NF}, Kc), zeros(Complex{NF}, Kc))
        cbar = zeros(Complex{NF}, S.spectrum, Kc)
        _legendre!(cbar, dg_n, dg_s, col, S)                            # adjoint of inverse Legendre: fwd Legendre
        wrapped_view(coeffs_bar, :, chunk).data .+= cbar.data
        c = c_end + 1
    end
    return coeffs_bar
end

# adjoint of grid→spec transform! w.r.t field: coeffs_bar → field_bar (accumulates into field_bar)
function grid2spec_pullback!(field_bar, coeffs_bar, S)
    NF = eltype(S)
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
        df_n = zeros(Complex{NF}, S.nfreq_max, Kc, nlat_half)
        df_s = zeros(Complex{NF}, S.nfreq_max, Kc, nlat_half)
        col = ColumnScratchMemory(zeros(Complex{NF}, Kc), zeros(Complex{NF}, Kc))
        _legendre!(df_n, df_s, wrapped_view(coeffs_bar, :, chunk), col, S)   # adjoint of fwd Legendre: inv Legendre
        df_n .*= dOmega                                                 # re-apply the ΔΩ weight
        df_s .*= dOmega
        df_n ./= scale                                                  # adjoint of fwd FFT: inv FFT of df/scale
        df_s ./= scale
        fbar = zeros(NF, S.grid, Kc)
        _fourier!(fbar, df_n, df_s, S)
        wrapped_view(field_bar, :, chunk).data .+= fbar.data
        c = c_end + 1
    end
    return field_bar
end

function augmented_primal(
        config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_chunked_spec2grid!)}, ::Type{<:Annotation},
        field::Duplicated, coeffs::Duplicated, scratch::Union{Const, Duplicated},
        S::Union{Const, MixedDuplicated}, unscale_coslat::Const,
    )
    func.val(field.val, coeffs.val, scratch.val, S.val, unscale_coslat.val)
    primal = needs_primal(config) ? field.val : nothing
    shadow = needs_shadow(config) ? field.dval : nothing
    return AugmentedReturn(primal, shadow, nothing)
end

function reverse(
        ::EnzymeRules.RevConfigWidth{1}, ::Const{typeof(_chunked_spec2grid!)}, ::Type{<:Annotation}, tape,
        field::Duplicated, coeffs::Duplicated, scratch::Union{Const, Duplicated},
        S::Union{Const, MixedDuplicated}, unscale_coslat::Const,
    )
    spec2grid_pullback!(coeffs.dval, field.dval, S.val; unscale_coslat = unscale_coslat.val)
    make_zero!(field.dval)      # the output cotangent has been propagated to coeffs
    return (nothing, nothing, nothing, nothing, nothing)
end

function augmented_primal(
        config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_chunked_grid2spec!)}, ::Type{<:Annotation},
        coeffs::Duplicated, field::Duplicated, scratch::Union{Const, Duplicated},
        S::Union{Const, MixedDuplicated},
    )
    func.val(coeffs.val, field.val, scratch.val, S.val)
    primal = needs_primal(config) ? coeffs.val : nothing
    shadow = needs_shadow(config) ? coeffs.dval : nothing
    return AugmentedReturn(primal, shadow, nothing)
end

function reverse(
        ::EnzymeRules.RevConfigWidth{1}, ::Const{typeof(_chunked_grid2spec!)}, ::Type{<:Annotation}, tape,
        coeffs::Duplicated, field::Duplicated, scratch::Union{Const, Duplicated},
        S::Union{Const, MixedDuplicated},
    )
    grid2spec_pullback!(field.dval, coeffs.dval, S.val)
    make_zero!(coeffs.dval)     # the output cotangent has been propagated to field
    return (nothing, nothing, nothing, nothing)
end

end
