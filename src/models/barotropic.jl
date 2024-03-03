export BarotropicModel, initialize!

"""
The BarotropicModel contains all model components needed for the simulation of the
barotropic vorticity equations. To be constructed like

    model = BarotropicModel(;spectral_grid, kwargs...)

with `spectral_grid::SpectralGrid` used to initalize all non-default components
passed on as keyword arguments, e.g. `planet=Earth(spectral_grid)`. Fields, representing
model components, are
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct BarotropicModel{
    NF<:AbstractFloat,
    DS<:DeviceSetup,
    PL<:AbstractPlanet,
    AT<:AbstractAtmosphere,
    CO<:AbstractCoriolis,
    FR<:AbstractForcing,
    DR<:AbstractDrag,
    PA<:AbstractParticleAdvection,
    IC<:InitialConditions,
    TS<:AbstractTimeStepper,
    ST<:SpectralTransform{NF},
    IM<:AbstractImplicit,
    HD<:AbstractHorizontalDiffusion,
    GE<:AbstractGeometry,
    OW<:AbstractOutputWriter,
    FB<:AbstractFeedback,
} <: Barotropic
    
    # GRID
    spectral_grid::SpectralGrid = SpectralGrid(nlev=1)
    geometry::GE = Geometry(spectral_grid)
    
    # DYNAMICS
    planet::PL = Earth(spectral_grid)
    atmosphere::AT = EarthAtmosphere(spectral_grid)
    coriolis::CO = Coriolis(spectral_grid)
    forcing::FR = NoForcing()
    drag::DR = NoDrag()
    particle_advection::PA = ParticleAdvection(spectral_grid)
    initial_conditions::IC = StartWithRandomVorticity()
    
    # NUMERICS
    device_setup::DS = DeviceSetup(CPUDevice())
    time_stepping::TS = Leapfrog(spectral_grid)
    spectral_transform::ST = SpectralTransform(spectral_grid)
    implicit::IM = NoImplicit(spectral_grid)
    horizontal_diffusion::HD = HyperDiffusion(spectral_grid)

    # OUTPUT
    output::OW = OutputWriter(spectral_grid,Barotropic)
    callbacks::Dict{Symbol,AbstractCallback} = Dict{Symbol,AbstractCallback}()
    feedback::FB = Feedback()
end

has(::Type{<:Barotropic}, var_name::Symbol) = var_name in (:vor,)
default_concrete_model(::Type{Barotropic}) = BarotropicModel

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for most fields, representing components, of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!`."""
function initialize!(model::Barotropic; time::DateTime = DEFAULT_DATE)
    (;spectral_grid) = model

    spectral_grid.nlev > 1 && @warn "Only nlev=1 supported for BarotropicModel, \
        SpectralGrid with nlev=$(spectral_grid.nlev) provided."

    # initialize components
    initialize!(model.time_stepping, model)
    initialize!(model.coriolis, model)
    initialize!(model.forcing, model)
    initialize!(model.drag, model)
    initialize!(model.horizontal_diffusion, model)

    # initial conditions
    prognostic_variables = PrognosticVariables(spectral_grid, model)
    initialize!(prognostic_variables, model.initial_conditions, model)
    prognostic_variables.clock.time = time       # set the time

    # particle advection
    initialize!(prognostic_variables.particles, model)

    diagnostic_variables = DiagnosticVariables(spectral_grid, model)
    return Simulation(prognostic_variables, diagnostic_variables, model)
end