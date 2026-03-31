abstract type AbstractLandGeometry <: AbstractModelComponent end

export LandGeometry
struct LandGeometry{VectorType, IntType} <: AbstractLandGeometry
    "[OPTION] Number of soil layers"
    nlayers::IntType

    "[OPTION] thickness of each soil layer [m]"
    layer_thickness::VectorType
end

const LandGeometryOrNothing = Union{LandGeometry, Nothing}

# default constructor
function LandGeometry(SG::SpectralGrid; nlayers = DEFAULT_NLAYERS_SOIL, layer_thickness = nothing)

    if isnothing(layer_thickness)
        (; NF) = SG

        # for two layers use the default soil layer thickness of MITgcm's 2-layer model
        if nlayers == 2
            layer_thickness = NF[0.2, 2]
        else
            layer_thickness = ones(NF, nlayers)
        end
    end

    return LandGeometry{typeof(layer_thickness), typeof(nlayers)}(nlayers, layer_thickness)
end

initialize!(::LandGeometry, model::PrimitiveEquation) = nothing

function Base.show(io::IO, geom::LandGeometry{V}) where {V}
    println(io, "$(geom.nlayers)-layer LandGeometry{$V}")
    return print(io, "└ layer_thickness: $(geom.layer_thickness)")
end

# because model components can be `nothing`, their constructor being `Nothing()`
Base.Nothing(::SpectralGrid, ::LandGeometry) = Nothing()
