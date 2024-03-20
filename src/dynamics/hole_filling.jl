abstract type AbstractHoleFilling end

export NoHoleFilling
struct NoHoleFilling <: AbstractHoleFilling end
NoHoleFilling(SG::SpectralGrid) = NoHoleFilling()
initialize!(::NoHoleFilling, ::PrimitiveWet) = nothing
hole_filling!(::AbstractGrid,::NoHoleFilling,::PrimitiveWet) = nothing

export ClipNegatives
struct ClipNegatives <: AbstractHoleFilling end
ClipNegatives(SG::SpectralGrid) = ClipNegatives()
initialize!(::ClipNegatives, ::PrimitiveWet) = nothing

# function barrier
function hole_filling!(
    A::AbstractGrid,
    H::ClipNegatives,
    model::PrimitiveWet
)
    hole_filling!(A, H)
end

function hole_filling!(A::AbstractGrid,::ClipNegatives)
    @inbounds for ij in eachgridpoint(A)
        A[ij] = max(A[ij], 0)
    end
end