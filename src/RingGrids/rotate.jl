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

rotate(grid::AbstractGridArray, args...) = rotate!(deepcopy(grid), args...)
rotate!(grid::AbstractGridArray, degree::Integer) = rotate!(grid, LongitudeRotation(degree))
rotate!(grid::AbstractGridArray, ::LongitudeRotation{0}) = grid
rotate!(grid::AbstractGridArray, ::LongitudeRotation{90}) = _rotate!(grid, 1)
rotate!(grid::AbstractGridArray, ::LongitudeRotation{180}) = _rotate!(grid, 2)
rotate!(grid::AbstractGridArray, ::LongitudeRotation{270}) = _rotate!(grid, 3)

function _rotate!(grid::AbstractGridArray, quarter::Integer)
    @assert quarter in (1, 2, 3) "Can only shift by 1, 2, or 3 quarters of a ring"

    for k in eachgrid(grid)
        for ring in eachring(grid)
            v = view(grid, ring, k)
            shift = quarter*(length(ring) ÷ 4)
            LinearAlgebra.circshift!(v, shift)
        end
    end

    return grid
end