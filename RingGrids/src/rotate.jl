abstract type AbstractLongitudeRotation{degree} end
struct LongitudeRotation{degree} <: AbstractLongitudeRotation{degree} end

function LongitudeRotation(degree::Integer)
    m = mod(degree, 360)
    m == 0   && return LongitudeRotation{0}()
    m == 90  && return LongitudeRotation{90}()
    m == 180 && return LongitudeRotation{180}()
    m == 270 && return LongitudeRotation{270}()
    throw(ArgumentError("Rotation degree must be multiples of 0, 90, 180, or 270; $degree given."))
end

import LinearAlgebra: rotate!
export rotate, rotate!

"""$(TYPEDSIGNATURES)
Rotate `field` eastward in longitude by `degree`, a multiple of 90˚.
Returns a rotated copy, see `rotate!` for the in-place version."""
rotate(field::AbstractField, args...) = rotate!(deepcopy(field), args...)

"""$(TYPEDSIGNATURES)
Rotate `field` eastward in longitude by `degree` in-place, a multiple of 90˚
(negative degrees rotate westward). Every ring is circularly shifted by the
corresponding number of quarter rings."""
rotate!(field::AbstractField, degree::Integer) = rotate!(field, LongitudeRotation(degree))
rotate!(field::AbstractField, ::LongitudeRotation{0}) = field
rotate!(field::AbstractField, ::LongitudeRotation{90}) = _rotate!(field, 1)
rotate!(field::AbstractField, ::LongitudeRotation{180}) = _rotate!(field, 2)
rotate!(field::AbstractField, ::LongitudeRotation{270}) = _rotate!(field, 3)

function _rotate!(field::AbstractField, quarter::Integer)
    @assert quarter in (0, 1, 2, 3) "Can only shift by 0, 1, 2, or 3 quarters of a ring"
    quarter == 0 && return field

    arch = architecture(field)
    ring_first, ring_length = eachring_on_architecture(arch, field.grid)

    # one thread per ring j (and layer k) as a circular shift is sequential within a ring
    worksize = (get_nlat(field), size(field)[2:end]...)
    launch!(arch, ArrayWorkOrder, worksize, rotate_kernel!, field, quarter, ring_first, ring_length)
    return field
end

# in-place reverse of field[a:b, k], used for the in-place circular shift below
@inline function reverse_ring!(field, k, a::Integer, b::Integer)
    @inbounds while a < b
        field[a, k], field[b, k] = field[b, k], field[a, k]
        a += 1
        b -= 1
    end
end

@kernel inbounds = true function rotate_kernel!(field, quarter, ring_first, ring_length)
    I = @index(Global, Cartesian)
    j = I[1]                                    # ring index
    k = CartesianIndex(Base.tail(Tuple(I)))     # all non-horizontal dimensions
    i0 = ring_first[j]                          # index of first grid point in ring j
    n = ring_length[j]
    shift = quarter * (n ÷ 4)

    if shift != 0                       # circshift! by `shift` via triple reversal, in-place
        reverse_ring!(field, k, i0, i0 + n - shift - 1)
        reverse_ring!(field, k, i0 + n - shift, i0 + n - 1)
        reverse_ring!(field, k, i0, i0 + n - 1)
    end
end
