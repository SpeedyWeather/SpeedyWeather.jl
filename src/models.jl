"""
    M = BarotropicModel(::Parameters,
                        ::Constants,
                        ::Geometry,
                        ::SpectralTransform,
                        ::HorizontalDiffusion)

The BarotropicModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration. In contrast to
`ShallowWaterModel` or `PrimitiveEquation` it does not contain a `Boundaries` struct
as not needed."""
struct BarotropicModel{NF<:AbstractFloat, D<:AbstractDevice} <: Barotropic
    parameters::Parameters
    constants::Constants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    device_setup::DeviceSetup{D}
end

has(::Type{<:Barotropic}, var_name::Symbol) = var_name in (:vor,)
default_model(::Type{Barotropic}) = BarotropicModel

"""
    M = ShallowWaterModel(  ::Parameters,
                            ::Constants,
                            ::Geometry,
                            ::SpectralTransform,
                            ::Boundaries,
                            ::HorizontalDiffusion)

The ShallowWaterModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration."""
struct ShallowWaterModel{NF<:AbstractFloat, D<:AbstractDevice} <: ShallowWater
    parameters::Parameters
    constants::Constants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::Implicit{NF}
    device_setup::DeviceSetup{D}
end

has(::Type{<:ShallowWater}, var_name::Symbol) = var_name in (:vor, :div, :pres)
default_model(::Type{ShallowWater}) = ShallowWaterModel

"""
    M = PrimitiveDryCoreModel( ::Parameters,
                                ::Constants,
                                ::Geometry,
                                ::SpectralTransform,
                                ::Boundaries,
                                ::HorizontalDiffusion
                                ::Implicit)

The PrimitiveDryCoreModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration."""
struct PrimitiveDryCoreModel{NF<:AbstractFloat,D<:AbstractDevice} <: PrimitiveDryCore
    parameters::Parameters
    constants::Constants{NF}
    parameterization_constants::ParameterizationConstants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::Implicit{NF}
    device_setup::DeviceSetup{D}
end

"""
    M = PrimitiveWetCoreModel( ::Parameters,
                                ::Constants,
                                ::Geometry,
                                ::SpectralTransform,
                                ::Boundaries,
                                ::HorizontalDiffusion
                                ::Implicit)

The PrimitiveWetCoreModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration."""
struct PrimitiveWetCoreModel{NF<:AbstractFloat,D<:AbstractDevice} <: PrimitiveWetCore
    parameters::Parameters
    constants::Constants{NF}
    parameterization_constants::ParameterizationConstants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::Implicit{NF}
    device_setup::DeviceSetup{D}
end

has(::Type{<:PrimitiveDryCore}, var_name::Symbol) = var_name in (:vor, :div, :temp, :pres)
has(::Type{<:PrimitiveWetCore}, var_name::Symbol) = var_name in (:vor, :div, :temp, :humid, :pres)
default_model(::Type{PrimitiveEquation}) = PrimitiveDryCoreModel
default_model(Model::Type{<:ModelSetup}) = Model

"""
    has(M::ModelSetup, var_name::Symbol)

Returns true if the model `M` has a prognostic variable `var_name`, false otherwise.
The default fallback is that all variables are included. 
"""
has(::Type{<:ModelSetup}, var_name::Symbol) = var_name in (:vor, :div, :temp, :humid, :pres)
has(M::ModelSetup, var_name) = has(typeof(M), var_name)