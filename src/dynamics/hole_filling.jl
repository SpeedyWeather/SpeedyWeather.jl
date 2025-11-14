abstract type AbstractHoleFilling <: AbstractModelComponent end

# no hole filling
hole_filling!(::AbstractField, ::Nothing, ::AbstractModel) = nothing

export ClipNegatives
struct ClipNegatives <: AbstractHoleFilling end
ClipNegatives(SG::SpectralGrid) = ClipNegatives()

# function barrier
function hole_filling!(
        q::AbstractField,
        H::ClipNegatives,
        model::AbstractModel,
    )
    return hole_filling!(q, H)
end

function hole_filling!(q::AbstractField, ::ClipNegatives)
    return @. q = max(q, 0)
end
