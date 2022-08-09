abstract type ModelSetup end

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
struct BarotropicModel{NF<:AbstractFloat} <: ModelSetup
    parameters::Parameters
    constants::Constants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
end

"""
    M = ShallowWaterModel(  ::Parameters,
                            ::Constants,
                            ::Geometry,
                            ::SpectralTransform,
                            ::Boundaries,
                            ::HorizontalDiffusion)

The ShallowWaterModel struct holds all other structs that contain precalculated constants, whether scalars or
arrays that do not change throughout model integration."""
struct ShallowWaterModel{NF<:AbstractFloat} <: ModelSetup
    parameters::Parameters
    constants::Constants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::Implicit{NF}
end

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
struct PrimitiveEquationModel{NF<:AbstractFloat} <: ModelSetup
    parameters::Parameters
    constants::Constants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::Implicit{NF}
end