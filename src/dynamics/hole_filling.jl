abstract type AbstractHoleFilling{NF} end

struct ClipNegatives{NF} <: AbstractHoleFilling{NF} end
ClipNegatives(SG::SpectralGrid) = ClipNegatives{SG.NF}()

# nothing to initialize
initialize!(H::ClipNegatives,::PrimitiveWet) = nothing

# function barrier
function hole_filling!(
    A::AbstractGrid,
    H::ClipNegatives,
    model::PrimitiveWet)

    hole_filling!(A,H)
end

function hole_filling!(A::AbstractGrid,H::ClipNegatives)
    @inbounds for ij in eachgridpoint(A)
        A[ij] = max(A[ij],0)
    end
end