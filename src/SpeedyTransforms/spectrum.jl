"""$(TYPEDSIGNATURES)
Compute the power spectrum of the spherical harmonic coefficients
`spec` (lower triangular matrix/array) of type `Complex{NF}`."""
function power_spectrum(
    spec::LowerTriangularArray;
    normalize::Bool=true,
)
    
    lmax, mmax = size(spec, OneBased, as=Matrix)    # 1-based max degree l, order m
    trunc = min(lmax, mmax)                         # consider only the triangle
                                                    # ignore higher degrees if lmax > mmax
    spectrum = zeros(NF, trunc)   

    # zonal modes m = 0, *1 as not mirrored at -m
    @inbounds for l in 1:trunc
        spectrum[l] = abs(spec[l, 1])^2
    end

    # other modes m > 0 *2 as complex conj at -m
    @inbounds for m in 2:trunc
        for l in m:trunc
            spectrum[l] += 2*abs(spec[l, m])^2
        end
    end

    # divide by number of orders m at l for normalization, "average power at l"
    if normalize
        @inbounds for l in 1:trunc  # 1-based degree, hence:
            spectrum[l] /= 2l-1     # 1/(2l + 1) but l → l-1
        end
    end

    return spectrum
end