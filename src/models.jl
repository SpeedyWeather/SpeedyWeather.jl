abstract type ModelSetup{D} end

"""
    M = BarotropicModel(::Parameters,
                        ::Constants,
                        ::GeoSpectral,
                        ::HorizontalDiffusion)

The BarotropicModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration. In contrast to
`ShallowWaterModel` or `PrimitiveEquationModel` it does not contain a `Boundaries` struct
as not needed."""
struct BarotropicModel{NF<:AbstractFloat, D<:AbstractDevice} <: ModelSetup{D}
    parameters::Parameters
    constants::Constants{NF}
    geospectral::GeoSpectral{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    device_setup::DeviceSetup{D}
end

"""
    M = ShallowWaterModel(  ::Parameters,
                            ::Constants,
                            ::GeoSpectral,
                            ::Boundaries,
                            ::HorizontalDiffusion)

The ShallowWaterModel struct holds all other structs that contain precalculated constants, whether scalars or
arrays that do not change throughout model integration."""
struct ShallowWaterModel{NF<:AbstractFloat, D<:AbstractDevice} <: ModelSetup{D}
    parameters::Parameters
    constants::Constants{NF}
    geospectral::GeoSpectral{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::Implicit{NF}
    device_setup::DeviceSetup{D}
end

"""
    M = PrimitiveEquationModel( ::Parameters,
                                ::Constants,
                                ::GeoSpectral,
                                ::Boundaries,
                                ::HorizontalDiffusion)

The PrimitiveEquationModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration."""
struct PrimitiveEquationModel{NF<:AbstractFloat,D<:AbstractDevice} <: ModelSetup{D}
    parameters::Parameters
    constants::Constants{NF}
    geospectral::GeoSpectral{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    device_setup::DeviceSetup{D}
end

# create a model hierarchy: BarotropicModel < ShallowWaterModel < PrimitiveEquationModel
# for instances
Base.:(<)(M1::BarotropicModel,M2::BarotropicModel) = false
Base.:(<)(M1::ShallowWaterModel,M2::ShallowWaterModel) = false
Base.:(<)(M1::PrimitiveEquationModel,M2::PrimitiveEquationModel) = false

Base.:(<)(M1::BarotropicModel,M2::ShallowWaterModel) = true
Base.:(<)(M1::BarotropicModel,M2::PrimitiveEquationModel) = true
Base.:(<)(M1::ShallowWaterModel,M2::PrimitiveEquationModel) = true

Base.:(<)(M1::ShallowWaterModel,M2::BarotropicModel) = false
Base.:(<)(M1::PrimitiveEquationModel,M2::BarotropicModel) = false
Base.:(<)(M1::PrimitiveEquationModel,M2::ShallowWaterModel) = false

# # and also for types
# Base.:(<)(::Type{ShallowWaterModel{NF1}},::Type{BarotropicModel{NF2}}) where {NF1,NF2} = false
# Base.:(<)(::Type{PrimitiveEquationModel{NF1}},::Type{BarotropicModel{NF2}}) where {NF1,NF2} = false
# Base.:(<)(::Type{PrimitiveEquationModel{NF1}},::Type{ShallowWaterModel{NF2}}) where {NF1,NF2} = false

# Base.:(<)(::Type{BarotropicModel{NF1}},::Type{ShallowWaterModel{NF2}}) where {NF1,NF2} = true
# Base.:(<)(::Type{BarotropicModel{NF1}},::Type{PrimitiveEquationModel{NF2}}) where {NF1,NF2} = true
# Base.:(<)(::Type{ShallowWaterModel{NF1}},::Type{PrimitiveEquationModel{NF2}}) where {NF1,NF2} = true

# Base.:(<)(::Type{BarotropicModel{NF1}},::Type{BarotropicModel{NF2}}) where {NF1,NF2} = false
# Base.:(<)(::Type{ShallowWaterModel{NF1}},::Type{ShallowWaterModel{NF2}}) where {NF1,NF2} = false
# Base.:(<)(::Type{PrimitiveEquationModel{NF1}},::Type{PrimitiveEquationModel{NF2}}) where {NF1,NF2} = false

# Base.:(<)(M::ModelSetup{NF},::Type{BarotropicModel}) where NF = typeof(M) < BarotropicModel{NF}
# Base.:(<)(M::ModelSetup{NF},::Type{ShallowWaterModel}) where NF = typeof(M) < ShallowWaterModel{NF}
# Base.:(<)(M::ModelSetup{NF},::Type{PrimitiveEquationModel}) where NF = typeof(M) < PrimitiveEquationModel{NF}

