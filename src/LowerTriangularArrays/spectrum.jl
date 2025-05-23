abstract type AbstractSpectrum end
struct CPU end
const DEFAULT_ARCHITECTURE = CPU

"""$(TYPEDSIGNATURES) 
Encodes the spectral trunction, orders and degrees of the spherical harmonics. 
Is used by every `LowerTriangularArray` and also defines the architecture on which the 
data of the `LowerTriangularArray` is stored.
"""
struct Spectrum{A, O, L, M} <: AbstractSpectrum
    lmax::Int
    mmax::Int
    architecture::A
    orders::O
    l_indices::L    # used by GPU kernels 
    m_indices::M    # used by GPU kernels
    lm_orders::O    # used by eachorder
end

"""
$(TYPEDSIGNATURES)
Create a `Spectrum` from the spectral truncation `lmax` and `mmax`. Both are 
assumed to be one-based, i.e. `lmax=5` and `mmax=5` will create a spectrum 
with T4 truncation. 
"""
function Spectrum(
    lmax::Integer,
    mmax::Integer;
    architecture = DEFAULT_ARCHITECTURE(),
  )
    return Spectrum(lmax, mmax, architecture, 
    [m:lmax for m in 1:mmax],   # orders 
    l_indices(lmax, mmax),      # l_indices 
    m_indices(lmax, mmax),      # m_indices
    lm_orders(lmax, mmax))      # lm_orders
end

"""
$(TYPEDSIGNATURES)
Create a `Spectrum` for the spectral truncation `trunc`. `trunc` is assumed to be 
zero-based, i.e. `trunc=4` will create a `Spectrum` with T4 truncation. 
"""
Spectrum(trunc::Integer; one_degree_more=true, kwargs...) = 
    Spectrum(trunc+1+one_degree_more, trunc+1; kwargs...)

Spectrum(; trunc::Integer, kwargs...) = Spectrum(trunc; kwargs...)

triangle_number(m::Integer) = m*(m+1)÷2
nonzeros(l::Integer, m::Integer) = l*m - triangle_number(m-1)
nonzeros(s::Spectrum) = nonzeros(s.lmax, s.mmax)
resolution(s::Spectrum) = (s.lmax, s.mmax)
truncation(s::Spectrum) = s.mmax - 1

function l_indices(lmax::Integer, mmax::Integer)
    l_vector = Vector{Int}(undef, nonzeros(lmax, mmax))
    lm = 0
    for m in 1:mmax
        for l in m:lmax
            lm += 1
            l_vector[lm] = l
        end
    end
    return l_vector
end

function m_indices(lmax::Integer, mmax::Integer)
    m_vector = Vector{Int}(undef, nonzeros(lmax, mmax))
    lm = 0
    for m in 1:mmax
        for l in m:lmax
            lm += 1
            m_vector[lm] = m
        end
    end
    return m_vector
end

function lm_orders(lmax::Integer, mmax::Integer)
    lm_orders = Vector{UnitRange{Int}}(undef, mmax)
    lm = 0
    for m in 1:mmax
        lm_orders[m] = lm+1:lm+(lmax-m+1)
        lm += lmax-m+1
    end
    return lm_orders
end

Base.:(==)(s1::Spectrum, s2::Spectrum) = 
    s1.lmax == s2.lmax && s1.mmax == s2.mmax

Base.show(io::IO, s::Spectrum) = print(io, "Spectrum(T$(s.mmax-1): (lmax=$(s.lmax), mmax=$(s.mmax)) on $(typeof(s.architecture)))")

eachorder(s::Spectrum) = s.lm_orders