# SHALLOW WATER MODEL
"""
    I = ImplicitShallowWater(   ξH₀::Vector{NF}
                                g∇²::Vector{NF}
                                ξg∇²::Vector{NF}
                                div_impl::Vector{NF})

Struct that holds various precomputed arrays for the semi-implicit correction to
prevent gravity waves from amplifying in the shallow water model."""
struct ImplicitShallowWater{NF<:AbstractFloat} <: AbstractImplicit{NF}
    ξH₀::Vector{NF}         # = 2αΔt*layer_thickness, store in vec for mutability
    g∇²::Vector{NF}         # = gravity*eigenvalues
    ξg∇²::Vector{NF}        # = 2αΔt*gravity*eigenvalues
    div_impl::Vector{NF}    # = 1 / (1-ξH₀*ξg∇²)
end

# PRIMITIVE EQUATION MODEL
struct ImplicitPrimitiveEq{NF<:AbstractFloat} <: AbstractImplicit{NF}
    ξ::Vector{NF}       # = 2α*dt, packed in a vector for mutability

    # the following arrays all have ξ = 2αΔt absorbed L <- ξL, etc.
    L::Vector{NF}       # used for δT = G_T + ξLδD
    R::Matrix{NF}       # used for δD = G_D + ξ(RδT + Uδlnpₛ)
    U::Vector{NF}       # see above
    W::Vector{NF}       # used for δlnpₛ = G_lnpₛ + ξWδD

    S⁻¹::Array{NF,3}    # used for δD = S⁻¹G, S = 1 + ξ²(RL + UW)
end