"""
$(TYPEDSIGNATURES)
Computes the scale for the adjoint/pullback of all discrete Fourier transforms. 
"""
function adjoint_scale(S::SpectralTransform)
    (; nlat_half, nlons, rfft_plans) = S
    nfreqs = [rfft_plan.osz[1] for rfft_plan in rfft_plans] # TODO: This works with FFTW, but does it with CUFFT as well?

    scale = zeros(Int, maximum(nfreqs), nlat_half)

    for i=1:nlat_half
        scale[1:nfreqs[i],:] = rfft_adjoint_scale(nfreqs[i], nlons[i])
    end 
    return scale
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