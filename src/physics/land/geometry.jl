abstract type AbstractLandGeometry <: AbstractModelComponent end

export LandGeometry
mutable struct LandGeometry{NF} <: AbstractLandGeometry
    "[OPTION] number of soil layers"
    nlayers::Int

    "[OPTION] thickness of each soil layer [m]"
    layer_thickness::Vector{NF}
end

# default constructor
function LandGeometry(SG::SpectralGrid, nlayers::Int; kwargs...)
    (; NF) = SG

    # for two layers use the default soil layer thickness of MITgcm's 2-layer model
    if nlayers == 2
        layer_thickness = NF[0.2, 2]
    else
        layer_thickness = ones(NF, nlayers)
    end

    return LandGeometry{NF}(nlayers, layer_thickness)
end

initialize!(::LandGeometry, model::PrimitiveEquation) = nothing

function Base.show(io::IO, geom::LandGeometry{T}) where {T}
    (; nlayers) = geom
    println(io, "$nlayers-layer LandGeometry{$T}")
    print(io, "└ layer_thickness: $(geom.layer_thickness)")
end
