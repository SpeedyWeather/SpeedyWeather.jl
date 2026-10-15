abstract type AbstractBoundaryLayer <: AbstractParameterization end

export BoundaryLayer

"""Composite type, containing surface roughness computation
and drag coefficient computation. Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct BoundaryLayer{SC, NW, SR, D} <: AbstractBoundaryLayer
    @component surface_roughness::SR
    @component surface_condition::SC
    @component neutral_wind_speed::NW
    @component drag::D
end

Adapt.@adapt_structure BoundaryLayer
function BoundaryLayer(
        SG::SpectralGrid;
        surface_condition = SurfaceCondition(SG),
        neutral_wind_speed = NeutralWindSpeed(SG),
        surface_roughness = LearnedSurfaceRoughness(SG),
        drag = BulkRichardsonDrag(SG),
    )
    return BoundaryLayer(; surface_condition, neutral_wind_speed, surface_roughness, drag)
end

Base.show(io::IO, M::BoundaryLayer) = Base.show(io, M, values = false)

function initialize!(BL::BoundaryLayer, model::AbstractModel)
    initialize!(BL.surface_condition, model)
    initialize!(BL.neutral_wind_speed, model)
    initialize!(BL.surface_roughness, model)
    initialize!(BL.drag, model)
    return nothing
end

# variables of boundary layer are determined by its subcomponents
variables(BL::BoundaryLayer, model::AbstractModel) = (
    variables(BL.surface_condition, model)...,
    variables(BL.neutral_wind_speed, model)...,
    variables(BL.surface_roughness, model)...,
    variables(BL.drag, model)...,
)

# just call the sub-paramterizations one after another
@propagate_inbounds function parameterization!(ij, vars, BL::BoundaryLayer, model)
    parameterization!(ij, vars, BL.surface_condition, model)  # takes the surface roughness from previous time step
    parameterization!(ij, vars, BL.neutral_wind_speed, model)
    parameterization!(ij, vars, BL.surface_roughness, model)
    parameterization!(ij, vars, BL.drag, model)
    return nothing
end