"""
Inverse Fourier transform in zonal direction.
"""
function fourier_inverse(   input::Array{Complex{NF},2},
                            G::GeoSpectral{NF},
                            scale::Bool=false) where {NF<:AbstractFloat}

    @unpack nlon, nlat, nlon_half, cosg⁻¹ = G.geometry
    @unpack mx = G.spectral

    # preallocate, TODO turn into plan_irfft for performance
    output = zeros(NF, nlon, nlat)
    nlonNF = convert(NF,nlon)

    for j in 1:nlat
        # Do inverse FFT then multiply by number of longitudes
        # add the truncated wavenumbers with zeros
        output[:,j] = nlonNF*irfft(vcat(input[:,j], zeros(Complex{NF}, nlon_half+1-mx)), nlon)
    end

    # Scale by cosine(lat) if needed
    if scale
        for j in 1:nlat
            output[:,j] *= cosg⁻¹[j]
        end
    end

    return output
end

"""
Fourier transform in the zonal direction.
"""
function fourier(   input::Array{NF,2},
                    G::GeoSpectral{NF}) where {NF<:AbstractFloat}

    @unpack nlon, nlat = G.geometry
    @unpack mx = G.spectral

    # preallocate output
    #TODO pass on output array as argument, turn fourier into fourier!
    output = zeros(Complex{NF}, mx, nlat)

    one_over_nlon = convert(NF,1/nlon)

    # Copy grid-point data into working array
    for j in 1:nlat
        # Do direct FFT then divide by number of longitudes
        #TODO use plan_rfft and use a type-flexible FFT instead of FFTW
        output[:,j] = one_over_nlon*rfft(input[:,j])[1:mx]  # truncate to mx ?
    end

    return output
end
