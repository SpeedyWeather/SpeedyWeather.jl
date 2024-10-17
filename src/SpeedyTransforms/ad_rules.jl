using Enzyme
import .EnzymeRules: reverse, augmented_primal
using .EnzymeRules

"""
$(TYPEDSIGNATURES)
Computes the scale for the adjoint/pullback of all discrete Fourier transforms. 
"""
function adjoint_scale(S::SpectralTransform)
    (; nlat_half, nlons, rfft_plans) = S
    nfreqs = [rfft_plan.osz[1] for rfft_plan in rfft_plans] # TODO: This works with FFTW, but does it with CUFFT as well?

    scale = zeros(Int, maximum(nfreqs), nlat_half)

    for i=1:nlat_half
        scale[1:nfreqs[i],i] = rfft_adjoint_scale(nfreqs[i], nlons[i])
    end 
    return reshape(scale, maximum(nfreqs), 1, nlat_half) # the scratch memory is (Freq x lvl x lat), so we insert 
                                                         # an additional dimension here for easier matrix multiply
end 

"""
$(TYPEDSIGNATURES)
Computes the scale for the adjoint/pullback of a real discrete fourier transform.
"""
function rfft_adjoint_scale(n_freq::Int, n_real::Int)
    if iseven(n_real)
        return [1; [2 for i=2:(n_freq-1)]; 1]
    else 
        return [1; [2 for i=2:n_freq]]
    end 
end 

function augmented_primal(config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_fourier!)}, ::Type{<:Const}, 
    f_north::Duplicated, f_south::Duplicated, grids::Duplicated, S::Const) 

    println("Augmented Primal Used") # TODO: remove 

    func.val(f_north.val, f_south.val, grids.val, S.val) # forward pass

    # save grids in tape if grids will be overwritten
    if overwritten(config)[4] # TODO: Not sure this is really necessary because grids won't ever get overwritten by this _fourier!
        tape = copy(grids.val)
    else
        tape = nothing
    end

    return AugmentedReturn(nothing, nothing, tape) # because the function actually returns nothing

end 

function reverse(config::EnzymeRules.RevConfigWidth{1}, func::Const{typeof(_fourier!)}, ::Type{<:Const}, tape,
    f_north::Duplicated, f_south::Duplicated, grids::Duplicated, S::Const)

    println("Custom Reverse Used") # TODO: remove

    # adjoint/vjp of FFT has a different scaling, compute it, apply it later to f_north, f_south
    scale = adjoint_scale(S.val)
    
    # retrieve grids value, either from original grids or from tape if grids may have been overwritten.
    gridsval = overwritten(config)[4] ? tape : grids.val

    # compute the actual vjp
    dgridval = zero(gridsval)
    _fourier!(dgridval, scale .* f_north.val, scale .* f_south.val, S.val) # inverse FFT (w/o normalization)
    grids.dval .+= dgridval 

    # no derivative wrt the f_north and f_south that were input because they are overwritten
    make_zero!(f_north.dval) 
    make_zero!(f_south.dval)

    # the function has no return values, so we also return nothing here
    return (nothing, nothing, nothing, nothing)
end
