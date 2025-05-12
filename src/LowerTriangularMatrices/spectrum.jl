abstract type AbstractSpectrum end
struct CPU end
const DEFAULT_ARCHITECTURE = CPU

struct Spectrum{A, D, DO} <: AbstractSpectrum
    lmax::Int
    mmax::Int
    architecture::A
    degrees::D
    degrees_orders::DO
end

function Spectrum(
    lmax::Integer,
    mmax::Integer,
    architecture = DEFAULT_ARCHITECTURE(),
    degrees = [m:lmax for m in 1:mmax],
    degrees_orders = degrees_orders(lmax, mmax),
)
    return Spectrum(lmax, mmax, architecture, degrees, degrees_orders)
end

triangle_number(m::Integer) = m*(m+1)÷2
nonzeros(l::Integer, m::Integer) = l*m - triangle_number(m-1)

function degrees_orders(lmax::Integer, mmax::Integer)
    degrees_orders = Vector{Tuple{Int, Int}}(undef, nonzeros(lmax, mmax))
    lm = 0
    for m in 1:mmax
        for l in m:lmax
            lm += 1
            degrees_orders[lm] = (l, m)
        end
    end
    return degrees_orders
end

struct LowerTriangularArray{T, N, S, ArrayType <: AbstractArray{T,N}} <: AbstractArray{T,N}
    spectrum::S
    data::ArrayType
end

Base.size(L::LowerTriangularArray) = size(L.data)
Base.getindex(L::LowerTriangularArray, args...) = getindex(L.data, args...)