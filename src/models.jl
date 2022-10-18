abstract type ModelSetup{NF,D} end

"""
    M = BarotropicModel(::Parameters,
                        ::Constants,
                        ::Geometry,
                        ::SpectralTransform,
                        ::HorizontalDiffusion)

The BarotropicModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration. In contrast to
`ShallowWaterModel` or `PrimitiveEquationModel` it does not contain a `Boundaries` struct
as not needed."""
struct BarotropicModel{NF<:AbstractFloat, D<:AbstractDevice} <: ModelSetup{NF,D}
    parameters::Parameters
    constants::Constants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    device_setup::DeviceSetup{D}
end

has(::Type{BarotropicModel{NF,D}}, var_name::Symbol) where {NF,D} = var_name in [:vor]

"""
    M = ShallowWaterModel(  ::Parameters,
                            ::Constants,
                            ::Geometry,
                            ::SpectralTransform,
                            ::Boundaries,
                            ::HorizontalDiffusion)

The ShallowWaterModel struct holds all other structs that contain precalculated constants, whether scalars or
arrays that do not change throughout model integration."""
struct ShallowWaterModel{NF<:AbstractFloat, D<:AbstractDevice} <: ModelSetup{NF,D}
    parameters::Parameters
    constants::Constants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::Implicit{NF}
    device_setup::DeviceSetup{D}
end

has(::Type{ShallowWaterModel{NF,D}}, var_name::Symbol) where {NF,D} = var_name in [:vor, :div, :pres]


"""
    M = PrimitiveEquationModel( ::Parameters,
                                ::Constants,
                                ::Geometry,
                                ::SpectralTransform,
                                ::Boundaries,
                                ::HorizontalDiffusion
                                ::Implicit)

The PrimitiveEquationModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration."""
struct PrimitiveEquationModel{NF<:AbstractFloat,D<:AbstractDevice} <: ModelSetup{NF,D}
    parameters::Parameters
    constants::Constants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::Implicit{NF}
    device_setup::DeviceSetup{D}
end

has(::Type{PrimitiveEquationModel{NF,D}}, var_name::Symbol) where {NF,D} = var_name in [:vor, :div, :temp, :humid, :pres]

"""
    has(::ModelSetup, var_name::Symbol)

Returns true if the model `M` has a prognostic variable `var_name`, false otherwise. The default fallback is that all variables are included. 
"""
has(::Type{ModelSetup{NF,D}}, var_name::Symbol) where {NF,D} = var_name in [:vor, :div, :temp, :humid, :pres] 
has(M::ModelSetup, var_name) = has(typeof(M), var_name)