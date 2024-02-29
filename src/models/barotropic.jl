export BarotropicModel, initialize!

"""
$(SIGNATURES)
The BarotropicModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration.
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct BarotropicModel{
    NF<:AbstractFloat,
    DS<:DeviceSetup,
    PL<:AbstractPlanet,
    AT<:AbstractAtmosphere,
    CO<:AbstractCoriolis,
    FR<:AbstractForcing,
    DR<:AbstractDrag,
    IC<:InitialConditions,
    TS<:AbstractTimeStepper,
    ST<:SpectralTransform{NF},
    IM<:AbstractImplicit,
    HD<:AbstractHorizontalDiffusion,
    GE<:AbstractGeometry,
    OW<:AbstractOutputWriter,
    FB<:AbstractFeedback
} <: Barotropic
    
    spectral_grid::SpectralGrid = SpectralGrid(nlev=1)
    device_setup::DS = DeviceSetup(CPUDevice())

    # DYNAMICS
    planet::PL = Earth(spectral_grid)
    atmosphere::AT = EarthAtmosphere(spectral_grid)
    coriolis::CO = Coriolis(spectral_grid)
    forcing::FR = NoForcing()
    drag::DR = NoDrag()
    initial_conditions::IC = StartWithRandomVorticity()

    # NUMERICS
    time_stepping::TS = Leapfrog(spectral_grid)
    spectral_transform::ST = SpectralTransform(spectral_grid)
    implicit::IM = NoImplicit(spectral_grid)
    horizontal_diffusion::HD = HyperDiffusion(spectral_grid)
    geometry::GE = Geometry(spectral_grid)

    # OUTPUT
    output::OW = OutputWriter(spectral_grid,Barotropic)
    feedback::FB = Feedback()
end

has(::Type{<:Barotropic}, var_name::Symbol) = var_name in (:vor,)
default_concrete_model(::Type{Barotropic}) = BarotropicModel

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for components of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!`."""
function initialize!(model::Barotropic; time::DateTime = DEFAULT_DATE)
    (;spectral_grid) = model

    spectral_grid.nlev > 1 && @warn "Only nlev=1 supported for BarotropicModel, \
        SpectralGrid with nlev=$(spectral_grid.nlev) provided."

    # slightly adjust model time step to be a convenient divisor of output timestep
    initialize!(model.time_stepping, model)

    # initialize components
    initialize!(model.coriolis, model)
    initialize!(model.forcing, model)
    initialize!(model.drag, model)
    initialize!(model.horizontal_diffusion, model)

    # initial conditions
    prognostic_variables = PrognosticVariables(spectral_grid, model)
    initialize!(prognostic_variables, model.initial_conditions, model)
    prognostic_variables.clock.time = time       # set the time

    diagnostic_variables = DiagnosticVariables(spectral_grid, model)
    return Simulation(prognostic_variables, diagnostic_variables, model)
end