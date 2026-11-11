# The implicit solver's reference temperature profile is recomputed when the time step changed or
# when its schedule executes. With Reactant `reinitialize!` (SpeedyWeatherReactantExt) decides this
# without branching on the traced clock, compared here against the CPU.

"""Initialized PrimitiveWetModel simulation on `arch` whose implicit solver uses `schedule`."""
function implicit_test_simulation(arch, schedule)
    spectral_grid = SpectralGrid(; architecture = arch, nlayers = nlayers_for_model(PrimitiveWetModel), truncation = TRUNCATION)
    reference_temperature = GlobalMeanReferenceTemperature{spectral_grid.NF}()     # recomputed on the schedule
    implicit = ImplicitPrimitiveEquation(spectral_grid; schedule = deepcopy(schedule), reference_temperature)
    return initialize!(create_model(PrimitiveWetModel, arch; implicit))
end

function set_temperature_profiles!(implicit, vars, reference, average)
    implicit.temp_profile .= reference
    vars.dynamics.average_temperature_profile .= average
    return nothing
end

"""Reference temperature profile of the implicit solver after `reinitialize!` on time step
`time_step_counter`, starting from a reference profile of 250 K and a global-mean profile of 300 K.
A leapfrog step (time step 2Δt) with the time step the implicit solver was initialized with, so that
only the schedule decides whether the reference profile is updated."""
function reference_temperature_after_reinitialize!(simulation, time_step_counter)
    vars, model = SpeedyWeather.unpack(simulation)
    (; implicit, time_stepping, architecture) = model
    NF = eltype(implicit.temp_profile)

    copy!(vars.prognostic.clock, Clock(step_counter = 2, time_step_counter = time_step_counter, Δt = time_stepping.Δt_millisec))
    implicit.Δt[] = SpeedyWeather.scalar_ref(architecture, 2 * time_stepping.Δt)[]
    @maybe_jit architecture set_temperature_profiles!(implicit, vars, NF(250), NF(300))
    @maybe_jit architecture SpeedyWeather.reinitialize!(implicit, model, vars)
    return Array(on_architecture(SpeedyWeather.CPU(), implicit.temp_profile))
end

@testset "Implicit reference temperature on schedule: $(nameof(typeof(schedule)))" for (schedule, cases) in (
        (StepSchedule(2), ((4, true), (3, false))),
        (PeriodicSchedule(every = Hour(2)), ((3, true), (4, false))),       # every 3 time steps of 40 min
        (PeriodicSchedule(every = Second(0)), ((3, false), (4, false))),    # never
    )
    simulation_cpu = implicit_test_simulation(SpeedyWeather.CPU(), schedule)
    simulation_reactant = implicit_test_simulation(ReactantDevice(), schedule)

    @testset "time step $time_step_counter" for (time_step_counter, scheduled) in cases
        reference_cpu = reference_temperature_after_reinitialize!(simulation_cpu, time_step_counter)
        reference_reactant = reference_temperature_after_reinitialize!(simulation_reactant, time_step_counter)
        @test all(==(scheduled ? 300 : 250), reference_cpu)
        @test reference_reactant == reference_cpu
    end
end
