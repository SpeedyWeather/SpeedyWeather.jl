module SpeedyWeatherReactantExt

using SpeedyWeather
using Reactant
using SpeedyWeather.DocStringExtensions

using SpeedyWeather: ReactantDevice

# initialize a BarotropicModel with ReactantDevice, only use @jit on some parts 
# TODO: might define a custom @jit_or_not macro for this to use in main code
function initialize!(model::BarotropicModel{SG,<:ReactantDevice}; time::DateTime = DEFAULT_DATE) where {SG <: SpectralGrid}
    (; spectral_grid) = model

    spectral_grid.nlayers > 1 && @error "Only nlayers=1 supported for BarotropicModel, \
        SpectralGrid with nlayers=$(spectral_grid.nlayers) provided."

    # initialize components
    @jit initialize!(model.geometry, model)
    @jit initialize!(model.time_stepping, model)
    @jit initialize!(model.coriolis, model)
    @jit initialize!(model.forcing, model)
    @jit initialize!(model.drag, model)
    @jit initialize!(model.horizontal_diffusion, model)
    @jit initialize!(model.random_process, model)
    @jit initialize!(model.particle_advection, model)

    # allocate prognostic and diagnostic variables
    prognostic_variables = PrognosticVariables(model)
    diagnostic_variables = DiagnosticVariables(model)
    # initialize particles (or other non-atmosphere prognostic variables)
    initialize!(prognostic_variables.particles, prognostic_variables, diagnostic_variables, model)

    # set the initial conditions
    @jit initialize!(prognostic_variables, model.initial_conditions, model)
    (; clock) = prognostic_variables
    clock.time = time       # set the current time
    clock.start = time      # and store the start time

    return Simulation(prognostic_variables, diagnostic_variables, model)
end

end
