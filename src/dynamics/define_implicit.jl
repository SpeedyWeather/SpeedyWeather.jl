# SHALLOW WATER MODEL
"""
    I = ImplicitShallowWater(   ξH₀::Vector,
                                g∇²::Vector,
                                ξg∇²::Vector,
                                S⁻¹::Vector)

Struct that holds various precomputed arrays for the semi-implicit correction to
prevent gravity waves from amplifying in the shallow water model."""
struct ImplicitShallowWater{NF<:AbstractFloat} <: AbstractImplicit{NF}
    ξH₀::Vector{NF}     # = 2αΔt*layer_thickness, store in vec for mutability
    g∇²::Vector{NF}     # = gravity*eigenvalues
    ξg∇²::Vector{NF}    # = 2αΔt*gravity*eigenvalues
    S⁻¹::Vector{NF}     # = 1 / (1-ξH₀*ξg∇²), implicit operator
end

# PRIMITIVE EQUATION MODEL
"""
    I = ImplicitPrimitiveEq(ξ::Vector,
                            R::Matrix,
                            U::Vector,
                            L::Matrix,
                            W::Vector,
                            S⁻¹::Matrix)

Struct that holds various precomputed arrays for the semi-implicit correction to
prevent gravity waves from amplifying in the primitive equation model."""
struct ImplicitPrimitiveEq{NF<:AbstractFloat} <: AbstractImplicit{NF}
    ξ::Base.RefValue{NF}# = 2α*Δt, packed in a RefValue for mutability

    # the following arrays all have ξ = 2αΔt absorbed L <- ξL, etc.
    R::Matrix{NF}       # divergence: used for δD = G_D + ξ(RδT + Uδlnpₛ)
    U::Vector{NF}       # divergence: see above
    L::Matrix{NF}       # tempereature: used for δT = G_T + ξLδD
    W::Vector{NF}       # surface pressure: used for δlnpₛ = G_lnpₛ + ξWδD

    # components of L
    L0::Vector{NF}
    L1::Matrix{NF}
    L2::Vector{NF}
    L3::Matrix{NF}
    L4::Vector{NF}

    # combined implicit operator S
    S::Matrix{NF}
    S⁻¹::Array{NF,3}    # used for δD = S⁻¹G, S = 1 - ξ²(RL + UW)
end