abstract type AbstractLandGeometry <: AbstractModelComponent end

export LandGeometry
struct LandGeometry{NF, VectorType} <: AbstractLandGeometry
    "[OPTION] number of soil layers"
    nlayers::Int

    "[OPTION] thickness of each soil layer [m]"
    layer_thickness::VectorType
end

# default constructor
function LandGeometry(SG::SpectralGrid; kwargs...)
    (; NF, VectorType) = SG
    nlayers = SG.nlayers_soil

    # for two layers use the default soil layer thickness of MITgcm's 2-layer model
    if nlayers == 2
        layer_thickness = NF[0.5, 4]
    else
        layer_thickness = ones(NF, nlayers)
    end

    return LandGeometry{NF, VectorType}(nlayers, layer_thickness)
end

initialize!(::LandGeometry, model::PrimitiveEquation) = nothing

function Base.show(io::IO, geom::LandGeometry{T, V}) where {T, V}
    (; nlayers) = geom
    println(io, "$nlayers-layer LandGeometry{$T, $V}")
    print(io, "└ layer_thickness: $(geom.layer_thickness)")
end
