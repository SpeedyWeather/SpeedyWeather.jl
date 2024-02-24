export ShallowWaterModel

"""
$(SIGNATURES)
The ShallowWaterModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration.
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct ShallowWaterModel{NF<:AbstractFloat, D<:AbstractDevice} <: ShallowWater
    spectral_grid::SpectralGrid = SpectralGrid(nlev=1)

    # DYNAMICS
    planet::AbstractPlanet = Earth()
    atmosphere::AbstractAtmosphere = EarthAtmosphere()
    coriolis::AbstractCoriolis = Coriolis(spectral_grid)
    forcing::AbstractForcing = NoForcing()
    drag::AbstractDrag = NoDrag()
    initial_conditions::InitialConditions = ZonalJet()
    orography::AbstractOrography{NF} = EarthOrography(spectral_grid)

    # NUMERICS
    time_stepping::TimeStepper = Leapfrog(spectral_grid)
    spectral_transform::SpectralTransform{NF} = SpectralTransform(spectral_grid)
    horizontal_diffusion::HorizontalDiffusion = HyperDiffusion(spectral_grid)
    implicit::AbstractImplicit = ImplicitShallowWater(spectral_grid)

    # INTERNALS
    geometry::AbstractGeometry = Geometry(spectral_grid)
    constants::DynamicsConstants{NF} = DynamicsConstants(spectral_grid,planet,atmosphere,geometry)
    device_setup::DeviceSetup{D} = DeviceSetup(CPUDevice())

    # OUTPUT
    output::AbstractOutputWriter = OutputWriter(spectral_grid,ShallowWater)
    feedback::AbstractFeedback = Feedback()
end

has(::Type{<:ShallowWater}, var_name::Symbol) = var_name in (:vor, :div, :pres)
default_concrete_model(::Type{ShallowWater}) = ShallowWaterModel

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for components of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!` and `model.implicit` which is done in `first_timesteps!`."""
function initialize!(model::ShallowWater;time::DateTime = DEFAULT_DATE)
    (;spectral_grid) = model

    spectral_grid.nlev > 1 && @warn "Only nlev=1 supported for ShallowWaterModel, \
                                SpectralGrid with nlev=$(spectral_grid.nlev) provided."

    # slightly adjust model time step to be a convenient divisor of output timestep
    initialize!(model.time_stepping,model)

    # initialize components
    initialize!(model.forcing,model)
    initialize!(model.drag,model)
    initialize!(model.horizontal_diffusion,model)
    initialize!(model.orography,model)

    # initial conditions
    prognostic_variables = PrognosticVariables(spectral_grid,model)
    initialize!(prognostic_variables,model.initial_conditions,model)
    prognostic_variables.clock.time = time       # set the time

    diagnostic_variables = DiagnosticVariables(spectral_grid,model)
    return Simulation(prognostic_variables,diagnostic_variables,model)
end