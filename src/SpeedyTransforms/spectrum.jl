"""$(TYPEDSIGNATURES)
Compute the power spectrum of the spherical harmonic coefficients
`spec` (lower triangular matrix/array) of type `Complex{NF}`.
For any additional dimensions in `spec`, the power spectrum is
computed along the first/spherical harmonic dimension."""
function power_spectrum(
    spec::LowerTriangularArray;
    normalize::Bool=true,
)
    
    lmax, mmax = size(spec, OneBased, as=Matrix)    # 1-based max degree l, order m
    trunc = min(lmax, mmax)                         # consider only the triangle
                                                    # ignore higher degrees if lmax > mmax
    spectrum = zeros(real(eltype(spec)), trunc, size(spec)[2:end]...)   

    @inbounds for k in eachmatrix(spec)
        # zonal modes m = 0, *1 as not mirrored at -m
        for l in 1:trunc    # use flat/vector indexing
            spectrum[l, k] = abs(spec[l, k])^2
        end

        # other modes m > 0 *2 as complex conj at -m
        for m in 2:trunc
            for l in m:trunc
                spectrum[l, k] += 2*abs(spec[l, m, k])^2
            end
        end
    end

    # divide by number of orders m at l for normalization, "average power at l"
    if normalize
        @inbounds for k in eachmatrix(spec)
            for l in 1:trunc            # 1-based degree, hence:
                spectrum[l, k] /= 2l-1  # 1/(2l + 1) but l → l-1
            end
        end
    end

    return spectrum
end