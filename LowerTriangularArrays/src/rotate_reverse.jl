import LinearAlgebra: rotate!
export rotate!

"""$(TYPEDSIGNATURES)
Rotate the field(s) represented by a LowerTriangularArray zonally by `degree`
by multiplication of the spherical harmonics by exp(-i*m*2π*degree/360),
with `m` the order of the spherical harmonic."""
function rotate!(L::LowerTriangularArray, degree::Real)
    mmax = size(L, 2, OneBased, as = Matrix)

    # complex rotation per order m, exp(-i*m*2π*degree/360) but more accurate;
    # precomputed in the precision of `degree` before conversion for accuracy
    NF = complex(eltype(L))

    # allocates a vector here but keeps the cispi out of kernels
    # which is not supported on some architectures (e.g. Metal)
    rotation = [convert(NF, cispi(-(m - 1) * degree / 180)) for m in 1:mmax]

    arch = architecture(L)
    launch!(
        arch, SpectralWorkOrder, size(L), rotate_kernel!,
        L.data, on_architecture(arch, rotation), L.spectrum.m_indices,
    )
    return L
end

@kernel inbounds = true function rotate_kernel!(data, @Const(rotation), @Const(m_indices))
    I = @index(Global, NTuple)
    m = m_indices[I[1]]     # order m of the harmonic at running index lm = I[1]
    data[I...] *= rotation[m]
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
    (; l_indices, m_indices) = L.spectrum
    arch = architecture(L)
    launch!(arch, SpectralWorkOrder, size(L), reverse_lat_kernel!, L.data, l_indices, m_indices)
    return L
end

@kernel inbounds = true function reverse_lat_kernel!(data, @Const(l_indices), @Const(m_indices))
    I = @index(Global, NTuple)
    lm = I[1]
    if isodd(l_indices[lm] + m_indices[lm])     # odd harmonics only
        data[I...] = -data[I...]
    end
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
