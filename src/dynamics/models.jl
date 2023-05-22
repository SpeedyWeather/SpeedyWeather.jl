"""
$(TYPEDSIGNATURES)
Simulation is a container struct to be used with `run!(::Simulation)`.
It contains
$(TYPEDFIELDS)"""
struct Simulation{Model<:ModelSetup}
    "define the current state of the model"
    prognostic_variables::PrognosticVariables

    "contain the tendencies and auxiliary arrays to compute them"
    diagnostic_variables::DiagnosticVariables

    "all parameters, constant at runtime"
    model::Model
end

"""
$(SIGNATURES)

The BarotropicModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration.

$(TYPEDFIELDS)"""
@kwdef struct BarotropicModel{NF<:AbstractFloat, D<:AbstractDevice} <: Barotropic
    "dictates resolution for many other components"
    spectral_grid::SpectralGrid = SpectralGrid()

    # PHYSICS 
    "contains physical and orbital characteristics"
    planet::AbstractPlanet = Earth()
    atmosphere::AbstractAtmosphere = EarthAtmosphere()
    forcing::AbstractForcing{NF} = NoForcing(spectral_grid)
    initial_conditions::InitialConditions = StartWithVorticity()

    # NUMERICS
    time_stepping::TimeStepper{NF} = Leapfrog(spectral_grid)
    spectral_transform::SpectralTransform{NF} = SpectralTransform(spectral_grid)
    horizontal_diffusion::HorizontalDiffusion{NF} = HyperDiffusion(spectral_grid)
    implicit::AbstractImplicit{NF} = NoImplicit(spectral_grid)

    # INTERNALS
    geometry::Geometry{NF} = Geometry(spectral_grid)
    constants::DynamicsConstants{NF} = DynamicsConstants(spectral_grid,planet,atmosphere,geometry)
    device_setup::DeviceSetup{D} = DeviceSetup(CPUDevice())

    # OUTPUT
    output::AbstractOutput = Output(spectral_grid,time_stepping)
    feedback::AbstractFeedback = Feedback(output,time_stepping)
end

has(::Type{<:Barotropic}, var_name::Symbol) = var_name in (:vor,)
default_concrete_model(::Type{Barotropic}) = BarotropicModel

function initialize!(model::BarotropicModel)
    initialize!(model.forcing,model)
    initialize!(model.horizontal_diffusion,model.time_stepping)

    prognostic_variables = initial_conditions(model)
    diagnostic_variables = DiagnosticVariables(model.spectral_grid)
    return Simulation(prognostic_variables,diagnostic_variables,model)
end

# """
#     M = ShallowWaterModel(  ::Parameters,
#                             ::DynamicsConstants,
#                             ::Geometry,
#                             ::SpectralTransform,
#                             ::Boundaries,
#                             ::HorizontalDiffusion)

# The ShallowWaterModel struct holds all other structs that contain precalculated constants,
# whether scalars or arrays that do not change throughout model integration."""
# struct ShallowWaterModel{NF<:AbstractFloat, D<:AbstractDevice} <: ShallowWater
#     parameters::Parameters
#     constants::DynamicsConstants{NF}
#     geometry::Geometry{NF}
#     spectral_transform::SpectralTransform{NF}
#     boundaries::Boundaries{NF}
#     horizontal_diffusion::HorizontalDiffusion{NF}
#     implicit::ImplicitShallowWater{NF}
#     device_setup::DeviceSetup{D}
# end

# has(::Type{<:ShallowWater}, var_name::Symbol) = var_name in (:vor, :div, :pres)
# default_concrete_model(::Type{ShallowWater}) = ShallowWaterModel

# """
#     M = PrimitiveDryModel(  ::Parameters,
#                                 ::DynamicsConstants,
#                                 ::Geometry,
#                                 ::SpectralTransform,
#                                 ::Boundaries,
#                                 ::HorizontalDiffusion
#                                 ::Implicit)

# The PrimitiveDryModel struct holds all other structs that contain precalculated constants,
# whether scalars or arrays that do not change throughout model integration."""
# struct PrimitiveDryModel{NF<:AbstractFloat,D<:AbstractDevice} <: PrimitiveDry
#     parameters::Parameters
#     constants::DynamicsConstants{NF}
#     parameterization_constants::ParameterizationConstants{NF}
#     geometry::Geometry{NF}
#     spectral_transform::SpectralTransform{NF}
#     boundaries::Boundaries{NF}
#     horizontal_diffusion::HorizontalDiffusion{NF}
#     implicit::ImplicitPrimitiveEq{NF}
#     device_setup::DeviceSetup{D}
# end

# """
#     M = PrimitiveWetModel(  ::Parameters,
#                                 ::DynamicsConstants,
#                                 ::Geometry,
#                                 ::SpectralTransform,
#                                 ::Boundaries,
#                                 ::HorizontalDiffusion
#                                 ::Implicit)

# The PrimitiveWetModel struct holds all other structs that contain precalculated constants,
# whether scalars or arrays that do not change throughout model integration."""
# struct PrimitiveWetModel{NF<:AbstractFloat,D<:AbstractDevice} <: PrimitiveWet
#     parameters::Parameters
#     constants::DynamicsConstants{NF}
#     parameterization_constants::ParameterizationConstants{NF}
#     geometry::Geometry{NF}
#     spectral_transform::SpectralTransform{NF}
#     boundaries::Boundaries{NF}
#     horizontal_diffusion::HorizontalDiffusion{NF}
#     implicit::ImplicitPrimitiveEq{NF}
#     device_setup::DeviceSetup{D}
# end

# has(::Type{<:PrimitiveDry}, var_name::Symbol) = var_name in (:vor, :div, :temp, :pres)
# has(::Type{<:PrimitiveWet}, var_name::Symbol) = var_name in (:vor, :div, :temp, :humid, :pres)
# default_concrete_model(::Type{PrimitiveEquation}) = PrimitiveDryModel
# default_concrete_model(Model::Type{<:ModelSetup}) = Model

# """
#     has(M::ModelSetup, var_name::Symbol)

# Returns true if the model `M` has a prognostic variable `var_name`, false otherwise.
# The default fallback is that all variables are included. 
# """
# has(::Type{<:ModelSetup}, var_name::Symbol) = var_name in (:vor, :div, :temp, :humid, :pres)
# has(M::ModelSetup, var_name) = has(typeof(M), var_name)

function Model(spectral_grid::SpectralGrid{WhichModel};kwargs...) where WhichModel
    return default_concrete_model(WhichModel)(;kwargs...)
end