const DEFAULT_ARCHITECTURE = CPU

abstract type AbstractSpectrum end

"""$(TYPEDSIGNATURES) 
Encodes the spectral trunction, orders and degrees of the spherical harmonics. 
Is used by every `LowerTriangularArray` and also defines the architecture on which the 
data of the `LowerTriangularArray` is stored.
"""
struct Spectrum{A, O, L} <: AbstractSpectrum
    lmax::Int
    mmax::Int
    architecture::A
    orders::O
    l_indices::L    # used by GPU kernels 
    m_indices::L    # used by GPU kernels
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

    orders = adapt(array_type(architecture), [m:lmax for m in 1:mmax])
    ls = adapt(array_type(architecture), l_indices(lmax, mmax))
    ms = adapt(array_type(architecture), m_indices(lmax, mmax))
    lm_orders_tuple = adapt(array_type(architecture), lm_orders(lmax, mmax))

    return Spectrum{typeof(architecture), 
                    typeof(orders), 
                    typeof(ls)}(
                    lmax, 
                    mmax, 
                    architecture, 
                    orders, 
                    ls, 
                    ms, 
                    lm_orders_tuple)
end

"""
$(TYPEDSIGNATURES)
Create a `Spectrum` for the spectral truncation `trunc`. `trunc` is assumed to be 
zero-based, i.e. `trunc=4` will create a `Spectrum` with T4 truncation. With
`one_degree_more==true` the `Spectrum` wil have an `lmax` increased by one, which 
is needed for spectral gradients. 
"""
Spectrum(trunc::Integer; one_degree_more=false, kwargs...) = 
    Spectrum(trunc+1+one_degree_more, trunc+1; kwargs...)

Spectrum(; trunc::Integer, kwargs...) = Spectrum(trunc; kwargs...)

"""
$(TYPEDSIGNATURES)
Create a `Spectrum` from another `Spectrum` but with a new architecture.
"""
Spectrum(spectrum::Spectrum; architecture::AbstractArchitecture=DEFAULT_ARCHITECTURE()) = 
    adapt(array_type(architecture), Spectrum(spectrum.lmax, 
            spectrum.mmax, 
            architecture, 
            spectrum.orders, 
            spectrum.l_indices, 
            spectrum.m_indices, 
            spectrum.lm_orders))

triangle_number(m::Integer) = m*(m+1)÷2
nonzeros(l::Integer, m::Integer) = l*m - triangle_number(m-1)
nonzeros(s::Spectrum) = nonzeros(s.lmax, s.mmax)
resolution(s::Spectrum) = (s.lmax, s.mmax)
truncation(s::Spectrum) = s.mmax - 1
orders(s::Spectrum) = s.orders
orders(s::Spectrum{<:GPU}) = Vector(s.orders) # on GPU transfer orders back to CPU first 

eachorder(s::Spectrum) = s.lm_orders

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

function Base.show(io::IO, S::Spectrum)
    println(io, "T$(S.mmax-1) Spectrum")
    println(io, "├ lmax=$(S.lmax) (degrees)")
    println(io, "├ mmax=$(S.mmax) (orders)")
    print(io,   "└ architecture: $(typeof(S.architecture))")
end

Architectures.ismatching(s::Spectrum, array_type::Type{<:AbstractArray}) = ismatching(s.architecture, array_type)
Architectures.ismatching(s::Spectrum, array::AbstractArray) = ismatching(s.architecture, typeof(array))

Adapt.@adapt_structure Spectrum

Architectures.architecture(s::Spectrum) = s.architecture
Architectures.on_architecture(architecture::AbstractArchitecture, s::Spectrum) = Spectrum(s; architecture)