abstract type ModelSetup{D} end

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
struct BarotropicModel{NF<:AbstractFloat, D<:AbstractDevice} <: ModelSetup{D}
    parameters::Parameters
    constants::Constants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    device_setup::DeviceSetup{D}
end

# GPU methods BarotropicModel

function Adapt.adapt_structure(to, bm::BarotropicModel)
    BarotropicModel(Adapt.adapt(to, bm.parameters),
                    Adapt.adapt(to, bm.constants),
                    Adapt.adapt(to, bm.geometry),
                    Adapt.adapt(to, bm.spectral_transform),
                    Adapt.adapt(to, bm.horizontal_diffusion),
                    Adapt.adapt(to, bm.device_setup))
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
struct ShallowWaterModel{NF<:AbstractFloat, D<:AbstractDevice} <: ModelSetup{D}
    parameters::Parameters
    constants::Constants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::Implicit{NF}
    device_setup::DeviceSetup{D}
end

# GPU methods ShallowWaterModel

function Adapt.adapt_structure(to, swm::ShallowWaterModel)
    ShallowWaterModel(Adapt.adapt(to, swm.parameters),
                      Adapt.adapt(to, swm.constants),
                      Adapt.adapt(to, swm.geometry),
                      Adapt.adapt(to, swm.spectral_transform),
                      Adapt.adapt(to, swm.boundaries),
                      Adapt.adapt(to, swm.horizontal_diffusion),
                      Adapt.adapt(to, swm.implicit),
                      Adapt.adapt(to, swm.device_setup))
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
struct PrimitiveEquationModel{NF<:AbstractFloat,D<:AbstractDevice} <: ModelSetup{D}
    parameters::Parameters
    constants::Constants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::Implicit{NF}
    device_setup::DeviceSetup{D}
end

# GPU methods PrimitiveEquationModel

function Adapt.adapt_structure(to, pem::PrimitiveEquationModel)
    PrimitiveEquationModel(Adapt.adapt(to, pem.parameters),
                           Adapt.adapt(to, pem.constants),
                           Adapt.adapt(to, pem.geometry),
                           Adapt.adapt(to, pem.spectral_transform),
                           Adapt.adapt(to, pem.boundaries),
                           Adapt.adapt(to, pem.horizontal_diffusion),
                           Adapt.adapt(to, pem.implicit),
                           Adapt.adapt(to, pem.device_setup))
end