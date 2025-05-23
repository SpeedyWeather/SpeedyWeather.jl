import LinearAlgebra: rotate!
export rotate!, rotate

abstract type AbstractLongitudeRotation{degree} end
struct LongitudeRotation{degree} <: AbstractLongitudeRotation{degree} end

function LongitudeRotation(degree::Integer)
    mod(degree, 360) == 0   && return LongitudeRotation{0}()
    mod(degree, 360) == 90  && return LongitudeRotation{90}()
    mod(degree, 360) == 180 && return LongitudeRotation{180}()
    mod(degree, 360) == 270 && return LongitudeRotation{270}()
    throw(ArgumentError("Rotation degree must be multiples of 0, 90, 180, or 270"))
end

rotate( field::AbstractField, args...) = rotate!(deepcopy(field), args...)
rotate!(field::AbstractField, degree::Integer) = rotate!(field, LongitudeRotation(degree))
rotate!(field::AbstractField, ::LongitudeRotation{0}) = field
rotate!(field::AbstractField, ::LongitudeRotation{90}) =  _rotate!(field, 1)
rotate!(field::AbstractField, ::LongitudeRotation{180}) = _rotate!(field, 2)
rotate!(field::AbstractField, ::LongitudeRotation{270}) = _rotate!(field, 3)

function _rotate!(field::AbstractField, quarter::Integer)
    @assert quarter in (1, 2, 3) "Can only shift by 1, 2, or 3 quarters of a ring"

    for k in eachlayer(field)
        for ring in eachring(field)
            v = view(field, ring, k)
            shift = quarter*(length(ring) ÷ 4)
            LinearAlgebra.circshift!(v, shift)
        end
    end

    return field
end