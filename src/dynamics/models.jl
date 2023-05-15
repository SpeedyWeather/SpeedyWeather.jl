"""
    model = BarotropicModel(::Parameters,
                        ::DynamicsConstants,
                        ::Geometry,
                        ::SpectralTransform,
                        ::HorizontalDiffusion)

The BarotropicModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration. In contrast to
`ShallowWaterModel` or `PrimitiveEquation` it does not contain a `Boundaries` struct
as not needed."""
Base.@kwdef struct BarotropicModel{NF<:AbstractFloat, D<:AbstractDevice} <: Barotropic
    
    # SpectralGrid object dictates resolution for many other components
    spectral_grid::SpectralGrid = SpectralGrid()

    # PHYSICS 
    planet::Planet = Earth()
    atmosphere::Atmosphere = EarthAthmosphere()

    # NUMERICS
    time_stepping::TimeIntegrator = LeapfrogSemiImplicit()
    spectral_transform::SpectralTransform{NF} = SpectralTransform(spectral_grid)
    horizontal_diffusion::HorizontalDiffusion{NF} = HyperDiffusion(spectral_grid)

    # INTERNALS
    geometry::Geometry{NF} = Geometry(spectral_grid)
    constants::DynamicsConstants{NF} = DynamicsConstants(spectral_grid,planet,atmosphere,time_stepping,geometry)
    device_setup::DeviceSetup{D} = DeviceSetup(CPUDevice())
end

has(::Type{<:Barotropic}, var_name::Symbol) = var_name in (:vor,)
default_concrete_model(::Type{Barotropic}) = BarotropicModel

"""
    M = ShallowWaterModel(  ::Parameters,
                            ::DynamicsConstants,
                            ::Geometry,
                            ::SpectralTransform,
                            ::Boundaries,
                            ::HorizontalDiffusion)

The ShallowWaterModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration."""
struct ShallowWaterModel{NF<:AbstractFloat, D<:AbstractDevice} <: ShallowWater
    parameters::Parameters
    constants::DynamicsConstants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::ImplicitShallowWater{NF}
    device_setup::DeviceSetup{D}
end

has(::Type{<:ShallowWater}, var_name::Symbol) = var_name in (:vor, :div, :pres)
default_concrete_model(::Type{ShallowWater}) = ShallowWaterModel

"""
    M = PrimitiveDryCoreModel(  ::Parameters,
                                ::DynamicsConstants,
                                ::Geometry,
                                ::SpectralTransform,
                                ::Boundaries,
                                ::HorizontalDiffusion
                                ::Implicit)

The PrimitiveDryCoreModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration."""
struct PrimitiveDryCoreModel{NF<:AbstractFloat,D<:AbstractDevice} <: PrimitiveDryCore
    parameters::Parameters
    constants::DynamicsConstants{NF}
    parameterization_constants::ParameterizationConstants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::ImplicitPrimitiveEq{NF}
    device_setup::DeviceSetup{D}
end

"""
    M = PrimitiveWetCoreModel(  ::Parameters,
                                ::DynamicsConstants,
                                ::Geometry,
                                ::SpectralTransform,
                                ::Boundaries,
                                ::HorizontalDiffusion
                                ::Implicit)

The PrimitiveWetCoreModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration."""
struct PrimitiveWetCoreModel{NF<:AbstractFloat,D<:AbstractDevice} <: PrimitiveWetCore
    parameters::Parameters
    constants::DynamicsConstants{NF}
    parameterization_constants::ParameterizationConstants{NF}
    geometry::Geometry{NF}
    spectral_transform::SpectralTransform{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
    implicit::ImplicitPrimitiveEq{NF}
    device_setup::DeviceSetup{D}
end

has(::Type{<:PrimitiveDryCore}, var_name::Symbol) = var_name in (:vor, :div, :temp, :pres)
has(::Type{<:PrimitiveWetCore}, var_name::Symbol) = var_name in (:vor, :div, :temp, :humid, :pres)
default_concrete_model(::Type{PrimitiveEquation}) = PrimitiveDryCoreModel
default_concrete_model(Model::Type{<:ModelSetup}) = Model

"""
    has(M::ModelSetup, var_name::Symbol)

Returns true if the model `M` has a prognostic variable `var_name`, false otherwise.
The default fallback is that all variables are included. 
"""
has(::Type{<:ModelSetup}, var_name::Symbol) = var_name in (:vor, :div, :temp, :humid, :pres)
has(M::ModelSetup, var_name) = has(typeof(M), var_name)