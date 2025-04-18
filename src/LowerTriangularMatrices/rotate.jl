import LinearAlgebra: rotate!
export rotate!

function rotate!(L::LowerTriangularArray, degree::Real)
    lmax, mmax = size(L, OneBased, as=Matrix)

    for k in eachmatrix(L)
        lm = 0
        for m in 1:mmax
            # complex rotation, exp(im*2π*degree/360) but more accurate
            o = convert(complex(eltype(L)), cispi(-(m-1)*degree/180))
            for l in m:lmax
                lm += 1
                L[lm, k] *= o
            end
        end
    end
    return L
end