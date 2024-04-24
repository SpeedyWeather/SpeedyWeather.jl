function power_spectrum(alms::LowerTriangularMatrix{Complex{NF}};
                        normalize::Bool=true) where NF
    
    lmax, mmax = size(alms)      # 1-based max degree l, order m
    trunc = min(lmax, mmax)      # consider only the triangle
                                # ignore higher degrees if lmax > mmax
    spectrum = zeros(NF, trunc)   

    # zonal modes m = 0, *1 as not mirrored at -m
    @inbounds for l in 1:trunc
        spectrum[l] = abs(alms[l, 1])^2
    end

    # other modes m > 0 *2 as complex conj at -m
    @inbounds for m in 2:trunc
        for l in m:trunc
            spectrum[l] += 2*abs(alms[l, m])^2
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