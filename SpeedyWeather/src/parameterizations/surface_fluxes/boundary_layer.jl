abstract type AbstractBoundaryLayer <: AbstractParameterization end

export BoundaryLayer

"""Composite boundary layer parameterization computing the drag coefficients
for surface fluxes of momentum, heat and humidity. Its subcomponents are
executed in the following order:

1. `surface_condition`: extrapolates temperature, density and wind to the surface
2. `neutral_wind_speed`: neutral (stability-corrected) surface wind speed
3. `surface_roughness`: roughness lengths for momentum, heat and humidity over land and ocean
4. `drag`: drag coefficients for momentum, heat and humidity over land and ocean

The wind slowdown in the surface condition uses fixed roughness lengths by default, or optionally
the surface roughness from the previous time step (see `LearnedWindSlowdown`).
This is not the default to avoid hidden state.
Fields are $(TYPEDFIELDS)"""
@parameterized @kwdef struct BoundaryLayer{SC, NW, SR, D} <: AbstractBoundaryLayer
    @component surface_condition::SC
    @component neutral_wind_speed::NW
    @component surface_roughness::SR
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
    parameterization!(ij, vars, BL.surface_condition, model)  # optionally uses surface roughness from previous time step
    parameterization!(ij, vars, BL.neutral_wind_speed, model)
    parameterization!(ij, vars, BL.surface_roughness, model, BL.neutral_wind_speed)
    parameterization!(ij, vars, BL.drag, model)
    return nothing
end