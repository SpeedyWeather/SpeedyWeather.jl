abstract type AbstractSpectrum end
struct CPU end
const DEFAULT_ARCHITECTURE = CPU

struct Spectrum{A, O, DO} <: AbstractSpectrum
    lmax::Int
    mmax::Int
    architecture::A
    orders::O
    degrees_orders::DO
end

function Spectrum(
    lmax::Integer,
    mmax::Integer;
    architecture = DEFAULT_ARCHITECTURE(),
    orders = [m:lmax for m in 1:mmax],
    degrees_orders = degrees_orders(lmax, mmax),
)
    return Spectrum(lmax, mmax, architecture, orders, degrees_orders)
end

triangle_number(m::Integer) = m*(m+1)÷2
nonzeros(l::Integer, m::Integer) = l*m - triangle_number(m-1)
nonzeros(s::Spectrum) = nonzeros(s.lmax, s.mmax)
resolution(s::Spectrum) = (s.lmax, s.mmax)
truncation(s::Spectrum) = s.mmax - 1

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

Base.:(==)(s1::Spectrum, s2::Spectrum) = 
    s1.lmax == s2.lmax && s1.mmax == s2.mmax

Base.show(io::IO, s::Spectrum) = print(io, "Spectrum(T$(s.mmax-1): (lmax=$(s.lmax), mmax=$(s.mmax)) on $(typeof(s.architecture)))")