abstract type AbstractLandGeometry <: AbstractModelComponent end

export LandGeometry
struct LandGeometry{VectorType} <: AbstractLandGeometry
    "[OPTION] thickness of each soil layer [m]"
    layer_thickness::VectorType
end

# default constructor
function LandGeometry(SG::SpectralGrid, nlayers::Int; layer_thickness=nothing, kwargs...)

    if isnothing(layer_thickness)
        (; NF) = SG
        nlayers = SG.nlayers_soil

        # for two layers use the default soil layer thickness of MITgcm's 2-layer model
        if nlayers == 2
            layer_thickness = NF[0.2, 2]
        else
            layer_thickness = ones(NF, nlayers)
        end
    end

    return LandGeometry(layer_thickness)
end

initialize!(::LandGeometry, model::PrimitiveEquation) = nothing

function Base.show(io::IO, geom::LandGeometry{V}) where {V}
    nlayers = length(geom.layer_thickness)
    println(io, "$nlayers-layer LandGeometry{$V}")
    print(io, "└ layer_thickness: $(geom.layer_thickness)")
end
