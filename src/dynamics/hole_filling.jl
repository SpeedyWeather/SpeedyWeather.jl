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
    q::AbstractGridArray,
    H::ClipNegatives,
    model::PrimitiveWet
)
    hole_filling!(q, H)
end

function hole_filling!(q::AbstractGridArray,::ClipNegatives)
    @. q = max(q, 0)
end