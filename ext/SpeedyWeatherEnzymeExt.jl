module SpeedyWeatherEnzymeExt

using SpeedyWeather
using Enzyme
using Enzyme.EnzymeCore
using SpeedyWeather.ProgressMeter
import .EnzymeRules: reverse, augmented_primal
using .EnzymeRules

# import all functions for which we define rules
import SpeedyWeather.SpeedyTransforms: _fourier!, ColumnScratchMemory

# Rules for SpeedyTransforms

# _fourier! 

# Computes the scale for the adjoint/pullback of all discrete Fourier transforms. 
function adjoint_scale(S::SpectralTransform)
    (; nlons, rfft_plans) = S
    (; nlat_half) = S.grid
    nfreqs = [rfft_plan.osz[1] for rfft_plan in rfft_plans] # TODO: This works with FFTW, but does it with cuFFT as well?

    scale = zeros(Int, maximum(nfreqs), 1, nlat_half) # the scratch memory is (Freq x lvl x lat), so we insert 
                                                      # an additional dimension here for easier matrix multiply

    for i=1:nlat_half
        scale[1:nfreqs[i],1,i] = rfft_adjoint_scale(nfreqs[i], nlons[i])
    end 

    # TODO: transfer array to GPU in case we are on GPU?
    return scale
end 

# Computes the scale for the adjoint/pullback of a real discrete fourier transform.
function rfft_adjoint_scale(n_freq::Int, n_real::Int)
    if iseven(n_real)
        return [1 < i < n_freq ? 2 : 1 for i=1:n_freq]
    else 
        return [1 < i ? 2 : 1 for i=1:n_freq]
    end 
end 

### Custom rule for _fourier!(f_north, f_north, grid, S)
function augmented_primal(config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_fourier!)}, ::Type{<:Const}, 
    f_north::Duplicated, f_south::Duplicated, grids::Duplicated{<:AbstractField}, scratch_memory::Duplicated{<:ColumnScratchMemory}, S::Union{Const, MixedDuplicated}) 

    func.val(f_north.val, f_south.val, grids.val, scratch_memory.val, S.val) # forward pass

    # save grids in tape if grids will be overwritten
    if overwritten(config)[4] # TODO: Not sure this is really necessary because grids won't ever get overwritten by this _fourier!
        tape = copy(grids.val)
    else
        tape = nothing
    end

    return AugmentedReturn(nothing, nothing, tape) # because the function actually returns nothing

end 

function reverse(config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_fourier!)}, ::Type{<:Const}, tape,
    f_north::Duplicated, f_south::Duplicated, grids::Duplicated{<:AbstractField}, scratch_memory::Duplicated{<:ColumnScratchMemory}, S::Union{Const, MixedDuplicated})

    # adjoint/jvp of FFT has a different scaling, compute it, apply it later to f_north, f_south
    scale = adjoint_scale(S.val)
    
    # retrieve grids value, either from original grids or from tape if grids may have been overwritten.
    gridsval = overwritten(config)[4] ? tape : grids.val

    # compute the adjoint
    dgridval = zero(gridsval)
    _fourier!(dgridval, f_north.dval ./ scale, f_south.dval ./ scale, scratch_memory.dval, S.val) # inverse FFT (w/o normalization)
    grids.dval .+= dgridval 

    # no derivative wrt the f_north and f_south that were input because they are overwritten
    make_zero!(f_north.dval) 
    make_zero!(f_south.dval)

    # the function has no return values, so we also return nothing here
    return (nothing, nothing, nothing, nothing, nothing)
end

### Custom rule for _fourier!(grid, f_north, f_south, S)
function augmented_primal(config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_fourier!)}, ::Type{<:Const}, 
    grids::Duplicated{<:AbstractField}, f_north::Duplicated, f_south::Duplicated, scratch_memory::Duplicated{<:ColumnScratchMemory}, S::Union{Const, MixedDuplicated}) 

    func.val(grids.val, f_north.val, f_south.val, scratch_memory.val, S.val) # forward pass

    # TODO: make an overwritten check here? 

    return AugmentedReturn(nothing, nothing, nothing) # because the function actually returns nothing

end 

function reverse(config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_fourier!)}, ::Type{<:Const}, tape,
    grids::Duplicated{<:AbstractField}, f_north::Duplicated, f_south::Duplicated, scratch_memory::Duplicated{<:ColumnScratchMemory}, S::Union{Const, MixedDuplicated})

    # adjoint/vjp of FFT has a different scaling, compute it, apply it later to f_north, f_south
    scale = adjoint_scale(S.val)
    
    # TODO: retrieve from tape here if overwritten? 

    # compute the adjoint # TODO: could we reuse the f_north.val for that a well? and not allocate here
    dfnorthval = zero(f_north.val)
    dfsouthval = zero(f_south.val)

    _fourier!(dfnorthval, dfsouthval, grids.dval, scratch_memory.dval, S.val) # inverse FFT (w/o normalization)

    f_north.dval .+= scale .* dfnorthval
    f_south.dval .+= scale .* dfsouthval 

    # no derivative wrt the grids that were input because they are overwritten
    make_zero!(grids.dval) 

    # the function has no return values, so we also return nothing here
    return (nothing, nothing, nothing, nothing, nothing)
end

###
# implement make_zero where the default one fails

# this lock is part of the ProgressMeter that's part of the Feedback of all models
@inline function Enzyme.make_zero(
    ::Type{ProgressMeter.ProgressCore}, 
    seen::IdDict, 
    prev::ProgressMeter.ProgressCore, 
    ::Val{copy_if_inactive} = Val(false),
)::ProgressMeter.ProgressCore where {copy_if_inactive} 
    return prev
end

end