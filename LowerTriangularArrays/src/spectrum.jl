abstract type AbstractSpectrum end

"""Encodes the spectral trunction, highested degree and order of the spherical harmonics (1-based).
Is used by every `LowerTriangularArray` and also defines the architecture on which the 
data of the `LowerTriangularArray` is stored. Fields are $(TYPEDFIELDS)"""
struct Spectrum{A, O, L, IntType} <: AbstractSpectrum
    "Highest order, meridional wavenumber, of the spherical harmonics (1-based)"
    lmax::IntType

    "Highest degree, zonal wavenumber, of the spherical harmonics (1-based)"
    mmax::IntType

    "Architecture used for LowerTriangularArrays created with this spectrum"
    architecture::A

    "[DERIVED] Precomputed to facilitate kernel launch or looping over elements"
    orders::O

    "[DERIVED] Precomputed to facilitate kernel launch or looping over elements"
    l_indices::L

    "[DERIVED] Precomputed to facilitate kernel launch or looping over elements"
    m_indices::L

    "[DERIVED] Precomputed to facilitate kernel launch or looping over elements"
    lm_orders::O
end

Adapt.@adapt_structure Spectrum

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

    orders = on_architecture(architecture, [m:lmax for m in 1:mmax])
    ls = on_architecture(architecture, l_indices(lmax, mmax))
    ms = on_architecture(architecture, m_indices(lmax, mmax))
    lm_orders_tuple = on_architecture(architecture, lm_orders(lmax, mmax))

    return Spectrum{
        typeof(architecture),
        typeof(orders),
        typeof(ls),
        typeof(lmax),
    }(
        lmax,
        mmax,
        architecture,
        orders,
        ls,
        ms,
        lm_orders_tuple
    )
end

"""$(TYPEDSIGNATURES)
Create a `Spectrum` for the spectral `truncation` (1-based).
E.g. `truncation = 4` will create a `Spectrum` with 4 zonal and meridional wavenumbers.
With kwarg `one_degree_more=true` the `Spectrum` wil have a meridional wavenumber `lmax`
increased by one, which is needed for spectral gradients."""
Spectrum(truncation::Integer; one_degree_more = false, kwargs...) =
    Spectrum(truncation + one_degree_more, truncation; kwargs...)

Spectrum(; truncation::Integer, kwargs...) = Spectrum(truncation; kwargs...)

"""$(TYPEDSIGNATURES)
Create a `Spectrum` from another `Spectrum` but with a new architecture."""
Spectrum(spectrum::Spectrum; architecture::AbstractArchitecture = DEFAULT_ARCHITECTURE()) =
    Spectrum(
    spectrum.lmax,
    spectrum.mmax,
    architecture,
    on_architecture(architecture, spectrum.orders),
    on_architecture(architecture, spectrum.l_indices),
    on_architecture(architecture, spectrum.m_indices),
    on_architecture(architecture, spectrum.lm_orders)
)

triangle_number(m::Integer) = m * (m + 1) ÷ 2
nonzeros(l::Integer, m::Integer) = l * m - triangle_number(m - 1)
nonzeros(s::Spectrum) = nonzeros(s.lmax, s.mmax)
resolution(s::Spectrum) = (s.lmax, s.mmax)
truncation(s::Spectrum, base = OneBased) = truncation(s, base)
truncation(s::Spectrum, ::Type{OneBased}) = s.mmax
truncation(s::Spectrum, ::Type{ZeroBased}) = s.mmax - 1
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
        lm_orders[m] = (lm + 1):(lm + (lmax - m + 1))
        lm += lmax - m + 1
    end
    return lm_orders
end

Base.:(==)(s1::Spectrum, s2::Spectrum) =
    s1.lmax == s2.lmax && s1.mmax == s2.mmax

function Base.show(io::IO, S::Spectrum)
    println(io, styled"T$(S.mmax) {warning:Spectrum}\{...\}")
    println(io, styled"├ {info:lmax} = $(S.lmax) {note:(degrees, 1-based)}")
    println(io, styled"├ {info:mmax} = $(S.mmax) {note:(orders, 1-based)}")
    print(io, styled"└ {info:architecture} = $(typeof(S.architecture))")
    return nothing
end

Architectures.ismatching(s::Spectrum, array_type::Type{<:AbstractArray}) = ismatching(s.architecture, array_type)
Architectures.ismatching(s::Spectrum, array::AbstractArray) = ismatching(s.architecture, typeof(array))
Architectures.architecture(s::Spectrum) = s.architecture
Architectures.on_architecture(architecture::AbstractArchitecture, s::Spectrum) = Spectrum(s; architecture)
Architectures.on_architecture(s::Spectrum, x) = on_architecture(architecture(s), x)

"""$(TYPEDSIGNATURES)
Iterator over all spherical harmonics in `S`, yielding `(l, m)` tuples of
degree `l` and order `m` (both 1-based) for every harmonic in the lower triangle.
To be used like

    for (l, m) in eachharmonic(S)
        L[l, m]
    end
"""
eachharmonic(S::Spectrum) = zip(S.l_indices, S.m_indices)
Base.eachindex(S::Spectrum) = Base.OneTo(length(S.l_indices))