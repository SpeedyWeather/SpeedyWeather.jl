abstract type ModelSetup end

"""
    M = BarotropicModel(::Parameters,
                        ::Constants,
                        ::GeoSpectral,
                        ::HorizontalDiffusion)

The BarotropicModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration. In contrast to
`ShallowWaterModel` or `PrimitiveEquationModel` it does not contain a `Boundaries` struct
as not needed."""
struct BarotropicModel{NF<:AbstractFloat} <: ModelSetup
    parameters::Parameters
    constants::Constants{NF}
    geospectral::GeoSpectral{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
end

"""
    M = ShallowWaterModel(  ::Parameters,
                            ::Constants,
                            ::GeoSpectral,
                            ::Boundaries,
                            ::HorizontalDiffusion)

The ShallowWaterModel struct holds all other structs that contain precalculated constants, whether scalars or
arrays that do not change throughout model integration."""
struct ShallowWaterModel{NF<:AbstractFloat} <: ModelSetup
    parameters::Parameters
    constants::Constants{NF}
    geospectral::GeoSpectral{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::Implicit{NF}
end

"""
    M = PrimitiveEquationModel( ::Parameters,
                                ::Constants,
                                ::GeoSpectral,
                                ::Boundaries,
                                ::HorizontalDiffusion)

The PrimitiveEquationModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration."""
struct PrimitiveEquationModel{NF<:AbstractFloat} <: ModelSetup
    parameters::Parameters
    constants::Constants{NF}
    geospectral::GeoSpectral{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
end