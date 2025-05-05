import LinearAlgebra: rotate!
export rotate!

"""$(TYPEDSIGNATURES)
Rotate the field(s) represented by a LowerTriangularArray zonally by `degree`
by multiplication of the spherical harmonics by exp(-i*m*2π*degree/360),
with `m` the order of the spherical harmonic."""
function rotate!(L::LowerTriangularArray, degree::Real)
    lmax, mmax = size(L, OneBased, as=Matrix)

    for k in eachmatrix(L)
        lm = 0
        for m in 1:mmax
            # complex rotation, exp(-i*m*2π*degree/360) but more accurate
            o = convert(complex(eltype(L)), cispi(-(m-1)*degree/180))
            for l in m:lmax
                lm += 1
                L[lm, k] *= o
            end
        end
    end
    return L
end

# also allow long names longitude, latitude
Base._reverse!(L::LowerTriangularArray, ::Val{:latitude}) = Base._reverse!(L, Val{:lat}())
Base._reverse!(L::LowerTriangularArray, ::Val{:longitude}) = Base._reverse!(L, Val{:lon}())

# turn into Val for multiple dispatch
Base._reverse!(L::LowerTriangularArray, sym::Symbol) = Base._reverse!(L, Val(sym))

"""$(TYPEDSIGNATURES)
Reverse the field represented by `L` in latitude when the element in `L`
represent the coefficients of the spherical harmonics. Reversal
in latitude direction is obtained by flipping the sign of the
odd harmonics. Reverses `L` in place."""
function Base._reverse!(L::LowerTriangularArray, ::Val{:lat})
    lmax, mmax = size(L, OneBased, as=Matrix)

    for k in eachmatrix(L)      # loop over any additional dimensions
        lm = 0
        for m in 1:mmax
            for l in m:lmax
                lm += 1
                if isodd(l+m)   # odd harmonics only
                    L[lm, k] = -L[lm, k]
                end
            end
        end
    end
    return L
end

"""$(TYPEDSIGNATURES)
Reverse the field represented by `L` in longitude (mirror at 0˚)
when the element in `L` represent the coefficients of the spherical harmonics.
Reversal in longitude direction is obtained by the complex harmonics.
Reverses `L` in place."""
function Base._reverse!(L::LowerTriangularArray, ::Val{:lon})
    L .= conj.(L)
    return L
end