abstract type ModelSetup end

"""
    M = BarotropicModel(::Parameters,
                        ::Constants,
                        ::GeoSpectral,
                        ::Boundaries,
                        ::HorizontalDiffusion)

The BarotropicModel struct holds all other structs that contain precalculated constants, whether scalars or
arrays that do not change throughout model integration."""
struct BarotropicModel{NF<:AbstractFloat} <: ModelSetup
    parameters::Parameters
    constants::Constants{NF}
    geospectral::GeoSpectral{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
end

struct ShallowWaterModel{NF<:AbstractFloat} <: ModelSetup
    parameters::Parameters
    constants::Constants{NF}
    geospectral::GeoSpectral{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
end

struct PrimitiveEquationModel{NF<:AbstractFloat} <: ModelSetup
    parameters::Parameters
    constants::Constants{NF}
    geospectral::GeoSpectral{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
end