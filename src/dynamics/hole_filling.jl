abstract type AbstractHoleFilling end

export ClipNegatives
struct ClipNegatives <: AbstractHoleFilling end
ClipNegatives(SG::SpectralGrid) = ClipNegatives()

# nothing to initialize
initialize!(::ClipNegatives, ::PrimitiveWet) = nothing

# function barrier
function hole_filling!(
    A::AbstractGrid,
    H::ClipNegatives,
    model::PrimitiveWet
)
    hole_filling!(A, H)
end

function hole_filling!(A::AbstractGrid, ::ClipNegatives)
    @inbounds for ij in eachgridpoint(A)
        A[ij] = max(A[ij], 0)
    end
end