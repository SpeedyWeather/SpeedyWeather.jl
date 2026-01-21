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

import LinearAlgebra: rotate!, circshift!
export rotate, rotate!

rotate(field::AbstractField, args...) = rotate!(deepcopy(field), args...)
rotate!(field::AbstractField, degree::Integer) = rotate!(field, LongitudeRotation(degree))
rotate!(field::AbstractField, ::LongitudeRotation{0}) = field
rotate!(field::AbstractField, ::LongitudeRotation{90}) = _rotate!(field, 1)
rotate!(field::AbstractField, ::LongitudeRotation{180}) = _rotate!(field, 2)
rotate!(field::AbstractField, ::LongitudeRotation{270}) = _rotate!(field, 3)

function _rotate!(field::AbstractField, quarter::Integer)
    @assert quarter in (0, 1, 2, 3) "Can only shift by 0, 1, 2, or 3 quarters of a ring"

    for k in eachlayer(field)
        for ring in eachring(field)
            v = view(field, ring, k)
            shift = quarter * (length(ring) ÷ 4)
            circshift!(v, shift)
        end
    end

    return field
end
