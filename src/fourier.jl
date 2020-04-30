"""
Inverse Fourier transform in zonal direction.
"""
function fourier_inverse(   input::Array{Complex{T},2},
                            geometry::Geometry{T},
                            scale::Bool=false) where {T<:AbstractFloat}

    @unpack nlon, nlat, mx, cosg⁻¹ = geometry

    # preallocate, TODO turn into plan_irfft for performance
    output = zeros(T, nlon, nlat)
    nlonT = T(nlon)

    for j in 1:nlat
        # Do inverse FFT then multiply by number of longitudes
        # add the truncated wavenumbers with zeros
        output[:,j] = nlonT*irfft(vcat(input[:,j], zeros(Complex{T}, div(96,2)+1-mx)), nlon)
    end

    # Scale by cosine(lat) if needed
    if scale
        for j in 1:nlat
            output[:,j] *= geometry.cosg⁻¹[j]
        end
    end

    return output
end

"""
Fourier transform in the zonal direction.
"""
function fourier(input::Array{T,2},geometry::Geometry{T}) where T

    @unpack nlon, nlat, mx = geometry

    # preallocate output
    #TODO pass on output array as argument, turn fourier into fourier!
    output = zeros(Complex{T}, mx, nlat)

    one_over_nlon = 1/T(nlon)

    # Copy grid-point data into working array
    for j in 1:nlat
        # Do direct FFT then divide by number of longitudes
        #TODO use plan_rfft and use a type-flexible FFT instead of FFTW
        output[:,j] = one_over_nlon*rfft(input[:,j])[1:mx]  # truncate to mx ?
    end

    return output
end
