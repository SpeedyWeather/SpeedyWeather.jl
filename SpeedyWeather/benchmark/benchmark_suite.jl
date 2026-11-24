using BenchmarkTools
using Statistics: median
using SpeedyWeather: synchronize, ReactantDevice, MatrixSpectralTransform, SpectralTransform

# The dense MatrixSpectralTransform is too large (memory and speed) above this truncation
const MAX_MATRIX_TRUNCATION = 150

const SPEEDY_CONFIGURATION = "The PrimitiveWet model in a configuration closer to SPEEDY.f90, at the resolutions " *
    "of the PrimitiveWet resolution sweep with the default transform. " *
    "As in SPEEDY.f90: prescribed seasonal climatologies of sea surface temperature (`SeasonalOceanClimatology`), " *
    "land surface temperature (`SeasonalLandTemperature`) and soil moisture (`SeasonalSoilMoisture`); daily-mean " *
    "insolation with a seasonal cycle on an eccentric orbit (`SolarZenithSeason`, `KeplerOrbit`); two-band shortwave " *
    "radiation with diagnostic clouds and seasonal ozone (`TwoBandShortwave`); albedos of 0.07 for the ocean and 0.6 for " *
    "sea ice, blended by sea-ice concentration (`OceanSeaIceAlbedo`), and a land albedo of bare soil and vegetation " *
    "with a snow cover linear in snow depth up to 60 mm water equivalent (`LandSnowAlbedo`); a near-surface wind of " *
    "0.95 times the lowermost-layer wind (`ConstantWindSlowdown`) with 5 m/s gusts (`ConstantWindGust`) and constant " *
    "roughness lengths (`ConstantSurfaceRoughness`); ∇⁸ hyperdiffusion with a 2.4-hour time scale at T32; a " *
    "40-minute time step at T32 with Robert–Williams filter. Without an equivalent of SPEEDY.f90's scheme in " *
    "SpeedyWeather: Betts–Miller convection (SPEEDY.f90: mass flux), one-band longwave radiation (SPEEDY.f90: " *
    "four bands), Richardson-number-dependent drag coefficients (SPEEDY.f90: constant coefficients with a " *
    "stability correction), `BulkRichardsonDiffusion` vertical diffusion, thermodynamic sea ice (SPEEDY.f90: " *
    "sea-ice climatology), a snow model (SPEEDY.f90: snow climatology) and the speed-limit drag. SpeedyWeather's " *
    "σ levels (Frierson et al. 2006) instead of SPEEDY.f90's, with which SpeedyWeather's shortwave absorption " *
    "moistens the stratosphere."

# Model components of the SPEEDY.f90-like configuration (`--SPEEDY` suite), see `SPEEDY_CONFIGURATION`,
# as `model_kwargs` of a `BenchmarkSuite`: types or functions of the spectral grid (`resolve_model_kwargs`).
# Every component of the PrimitiveWetModel and all their subcomponents are set explicitly, so that a
# changed default of the model does not change this configuration. Left to the benchmark as for all
# other suites: the spectral transform, output, feedback, diagnostics and callbacks. Parameters are the
# components' defaults unless set here, mostly to the values of SPEEDY.f90 (github.com/samhatfield/speedy.f90).
speedy_components(::Type{<:PrimitiveWet}) = (;
    # DYNAMICS
    geometry = speedy_geometry,
    planet = speedy_planet,
    atmosphere = EarthAtmosphere,
    coriolis = Coriolis,
    geopotential = Geopotential,
    adiabatic_conversion = AdiabaticConversion,
    particle_advection = nothing,
    initial_conditions = spectral_grid -> (;
        zonal = BalancedZonalState(spectral_grid),
        humid = ConstantRelativeHumidity(spectral_grid),
    ),
    forcing = nothing,
    drag = SpeedLimitDrag,                  # SPEEDY.f90 instead: 30-day drag on the zonal-mean wind in the top layer
    random_process = nothing,

    # BOUNDARY CONDITIONS
    orography = EarthOrography,
    land_sea_mask = EarthLandSeaMask,
    ocean = SeasonalOceanClimatology,       # SPEEDY.f90: SST climatology (+ observed anomalies)
    sea_ice = ThermodynamicSeaIce,          # SPEEDY.f90: sea-ice climatology, not available in SpeedyWeather
    land = speedy_land,

    # PHYSICS
    solar_zenith = speedy_solar_zenith,
    albedo = speedy_albedo,
    boundary_layer = speedy_boundary_layer,
    vertical_diffusion = BulkRichardsonDiffusion,
    surface_momentum_flux = spectral_grid -> SurfaceMomentumFlux(spectral_grid; time_integration = RationalFluxes()),
    surface_heat_flux = spectral_grid -> SurfaceHeatFlux(
        spectral_grid;
        ocean = SurfaceOceanHeatFlux(spectral_grid; time_integration = RationalFluxes()),
        land = SurfaceLandHeatFlux(spectral_grid; time_integration = RationalFluxes()),
    ),
    surface_humidity_flux = spectral_grid -> SurfaceHumidityFlux(
        spectral_grid;
        ocean = SurfaceOceanHumidityFlux(spectral_grid; time_integration = RationalFluxes()),
        land = SurfaceLandHumidityFlux(spectral_grid; time_integration = RationalFluxes()),
    ),
    large_scale_condensation = ImplicitCondensation,
    convection = spectral_grid -> BettsMillerConvection(spectral_grid; entrainment = LinearEntrainment(spectral_grid)),
    radiation = speedy_radiation,
    greenhouse_gases = spectral_grid -> (; co2 = CO2(spectral_grid, 360)),     # [ppm], constant as in SPEEDY.f90
    stochastic_physics = nothing,
    custom_parameterization = nothing,

    # NUMERICS
    time_stepping = speedy_time_stepping,
    implicit = spectral_grid -> ImplicitPrimitiveEquation(
        spectral_grid;
        centering = 0.5,                    # SPEEDY.f90's alph = 0.5
        reference_temperature = IsothermalReferenceTemperature{spectral_grid.NF}(),
        schedule = PeriodicSchedule(every = Day(1)),
    ),
    horizontal_diffusion = spectral_grid -> HyperDiffusion(
        spectral_grid;
        power = 4,                          # ∇⁸ as SPEEDY.f90's npowhd = 4
        time_scale = Minute(144),           # SPEEDY.f90's thd = 2.4 hours, scaled with resolution from T32
        resolution_scaling = 1,
        # divergence: SpeedyWeather's stronger damping, scale-selective in the troposphere
        # and ∇⁴ in the stratosphere as a sponge, SPEEDY.f90 instead: thdd = 2.4 hours
        # and an additional ∇² diffusion (12 hours) of vorticity, divergence and temperature
        time_scale_div = Hour(1),
        power_stratosphere = 2,
        tapering_σ = 0.2,
    ),
    vertical_advection = CenteredVerticalAdvection,
    hole_filling = ClipNegatives,
)

# Fields of the model that `speedy_components` leaves to the benchmark, as for all other suites
const SPEEDY_UNSET_FIELDS = (
    :spectral_grid, :architecture, :dynamics, :dynamics_only, :tracers, :spectral_transform, :output,
    :diagnostics, :callbacks, :feedback, :core_components, :parameterizations, :params,
)

# Error if `Model` has a field that `speedy_components` neither sets nor leaves to the benchmark, e.g.
# a new component, which would otherwise silently take its default in the SPEEDY.f90-like configuration
function check_speedy_components(Model)
    unset = setdiff(fieldnames(Model), keys(speedy_components(Model)), SPEEDY_UNSET_FIELDS)
    isempty(unset) || error(
        "speedy_components($Model) does not set $(Tuple(unset)): set them explicitly (benchmark_suite.jl) " *
            "so that the SPEEDY.f90-like configuration does not change with model defaults."
    )
    return nothing
end

# σ levels of Frierson et al. 2006, not SPEEDY.f90's (half levels 0, 0.05, 0.14, 0.26, 0.42, 0.6, 0.77,
# 0.9, 1 with 8 layers): with those, the shortwave water vapour absorption of SpeedyWeather (calibrated on
# these levels in #1414, SPEEDY.f90 has none in its top layer) moistens the top layer, 0.002 to 0.05 g/kg
# in the first year at T32L8 and still rising, while it stays below 0.02 g/kg with these levels
speedy_geometry(spectral_grid) = Geometry(spectral_grid; vertical_coordinates = FriersonSigmaCoordinates(spectral_grid))

# Earth with daily-mean insolation (no daily cycle) as in SPEEDY.f90, seasonal cycle on an eccentric
# orbit, also SPEEDY.f90 computes the Sun-Earth distance factor (via Spencer 1971's Fourier series)
speedy_planet(spectral_grid) = Earth(
    spectral_grid; daily_cycle = false, seasonal_cycle = true, orbit = KeplerOrbit(spectral_grid)
)

# daily-mean cos of the zenith angle as in `WhichZenith` for `daily_cycle = false`
# (`KeplerOrbit` computes the declination itself, the declination function is then not used)
speedy_solar_zenith(spectral_grid) = SolarZenithSeason(
    spectral_grid; seasonal_cycle = true,
    solar_declination = SpeedyWeather.SinSolarDeclination(speedy_planet(spectral_grid))
)

# land with SPEEDY.f90's prescribed climatologies of land surface temperature and soil moisture
# (SPEEDY.f90 adds anomalies of a slab land model relaxed to its climatology), the snow model as
# SpeedyWeather has no snow climatology (SPEEDY.f90's snow is prescribed)
function speedy_land(spectral_grid)
    geometry = LandGeometry(spectral_grid, nlayers = 2)
    return LandModel(;
        spectral_grid,
        geometry,
        thermodynamics = LandThermodynamics(spectral_grid, geometry),
        temperature = SeasonalLandTemperature(spectral_grid, geometry),
        soil_moisture = SeasonalSoilMoisture(spectral_grid, geometry),
        snow = SnowModel(spectral_grid, geometry),
        vegetation = VegetationClimatology(spectral_grid, geometry),
        rivers = nothing,
    )
end

# SPEEDY.f90's albedos: ocean 0.07 and sea ice 0.6 blended by sea-ice concentration, and over land
# snow cover = min(1, snow depth / 60 mm water equivalent). SPEEDY.f90 reads the snow-free land
# albedo from a climatology and blends it with snow (0.6), here vegetation and bare soil with
# additive snow albedo (0.4, so ≈0.6 under full snow cover over vegetation as in SPEEDY.f90)
speedy_albedo(spectral_grid) = OceanLandAlbedo(
    spectral_grid;
    ocean = OceanSeaIceAlbedo(spectral_grid; albedo_ocean = 0.07, albedo_ice = 0.6),
    land = LandSnowAlbedo(
        spectral_grid; snow_cover = LinearSnowCover(), snow_depth_scale = 0.06, albedo_snow = 0.4
    ),
)

# SPEEDY.f90's constant wind slowdown (0.95) and gusts (5 m/s), Kucharski et al. 2006, constant
# surface roughness lengths (SPEEDY.f90 has constant drag coefficients instead, `ConstantDrag` has
# no stability correction and one coefficient for momentum, heat and humidity), no neutral wind
# speed (only used by the learned ocean roughness)
speedy_boundary_layer(spectral_grid) = BoundaryLayer(
    spectral_grid;
    surface_condition = SurfaceCondition(
        spectral_grid;
        wind_slowdown = ConstantWindSlowdown(spectral_grid, wind_slowdown = 0.95),
        wind_gust = ConstantWindGust(spectral_grid, wind_gust = 5),
    ),
    neutral_wind_speed = nothing,
    surface_roughness = ConstantSurfaceRoughness(spectral_grid),
    drag = BulkRichardsonDrag(spectral_grid),
)

# SPEEDY.f90's two-band shortwave scheme with its diagnostic clouds and ozone (SpeedyWeather's
# absorptivities, calibrated in #1414) and one-band longwave (SPEEDY.f90: four bands)
speedy_radiation(spectral_grid) = Radiation(
    spectral_grid;
    shortwave = TwoBandShortwave(
        spectral_grid;
        clouds = DiagnosticClouds(spectral_grid),
        transmissivity = TwoBandShortwaveTransmissivity(spectral_grid),
        ozone = SeasonalOzone(spectral_grid),
        radiative_transfer = TwoBandShortwaveRadiativeTransfer(spectral_grid),
    ),
    longwave = OneBandLongwave(
        spectral_grid;
        transmissivity = ByrneOGormanLongwaveTransmissivity(spectral_grid),
        radiative_transfer = OneBandLongwaveRadiativeTransfer(spectral_grid, stratospheric_emissivity = 0.05),
    ),
)

# SPEEDY.f90's leapfrog: 40 min time step at T32 (36 steps/day), Robert filter rob = 0.05 and
# Williams filter wil = 0.53; SpeedyWeather's `robert_filter` is ν in Williams (2009) with the
# displacement ν/2, so 0.1 here is SPEEDY.f90's rob = 0.05
function speedy_time_stepping(spectral_grid)
    Δt_at_T32 = Minute(40)
    adjust_with_output = true
    return Leapfrog(
        spectral_grid;
        Δt_at_T32,
        adjust_with_output,
        robert_filter = 0.1,
        williams_filter = 0.53,
        ocean = EulerForward(spectral_grid; Δt_at_T32, adjust_with_output),
        land = EulerForward(spectral_grid; Δt_at_T32, adjust_with_output),
        noise = nothing,
    )
end

# Build the model with the spectral transform of `transform_kind` for `arch`
# (`:default` for FFT+Legendre, `:matrix` for MatrixSpectralTransform).
function build_model(Model, spectral_grid, arch; transform_kind::Symbol = :default, kwargs...)
    if arch isa ReactantDevice
        transform_kind in (:default, :matrix) || error("Unknown transform_kind = $transform_kind. Use :default or :matrix.")
        # `:default` explicitly the FFT-based SpectralTransform: the models' default on Reactant
        # (`WhichTransform`) is the MatrixSpectralTransform at every resolution
        S = transform_kind == :matrix ? MatrixSpectralTransform(spectral_grid) : SpectralTransform(spectral_grid)
        # no progress meter or output: they would read the clock from the device every time step,
        # synchronizing it. Passed last as repeated keyword arguments silently take the last value
        return Model(spectral_grid; kwargs..., spectral_transform = S, feedback = nothing, output = nothing)
    elseif transform_kind == :matrix
        M = MatrixSpectralTransform(spectral_grid)
        return Model(spectral_grid; spectral_transform = M, kwargs...)
    elseif transform_kind == :default
        return Model(spectral_grid; kwargs...)
    else
        error("Unknown transform_kind = $transform_kind. Use :default or :matrix.")
    end
end

# Format SYPD: one decimal digit if < 10, integer otherwise. Used by both the
# per-arch tables and the cross-arch overview so the column reads consistently.
function format_sypd(s)
    (s isa Number && isfinite(s)) || return "0"
    return s < 10 ? string(round(s; digits = 1)) : string(round(Int, s))
end

# Format a duration in seconds, "—" if not available (e.g. compile time without Reactant)
format_seconds(t) = isfinite(t) ? @sprintf("%.1f s", t) : "—"

# Spread of the timed runs, (slowest - fastest) / median wallclock time, from their SYPD
# (∝ 1/time). "—" for a single run.
function format_spread(SYPD, SYPD_min, SYPD_max, nrepeats)
    nrepeats > 1 || return "—"
    return @sprintf("%.0f%%", 100 * (SYPD / SYPD_min - SYPD / SYPD_max))
end

# Resolve a NamedTuple of per-run model kwargs against `spectral_grid`.
# Convention: a `Type` value is instantiated with `spectral_grid` (every
# component has a `T(spectral_grid)` constructor); a
# `Function` value is called with `spectral_grid` (for parameterised cases
# like `sg -> WhichZenith(sg, my_planet)`); everything else is passed through.
_resolve_kwarg_value(v::Type, sg) = v(sg)
_resolve_kwarg_value(v::Function, sg) = v(sg)
_resolve_kwarg_value(v, _) = v
function resolve_model_kwargs(nt::NamedTuple, spectral_grid)
    return NamedTuple{keys(nt)}(map(v -> _resolve_kwarg_value(v, spectral_grid), values(nt)))
end

abstract type AbstractBenchmarkSuite end
@kwdef mutable struct BenchmarkSuite <: AbstractBenchmarkSuite
    title::String
    # written below the title in the README, e.g. the configuration of the models
    description::String = ""
    nruns::Int = 1
    model::Vector = fill(PrimitiveWetModel, nruns)
    NF::Vector = fill(SpeedyWeather.DEFAULT_NF, nruns)
    truncation::Vector{Int} = fill(SpeedyWeather.DEFAULT_TRUNCATION, nruns)
    nlayers::Vector{Int} = default_nlayers(model)
    Grid::Vector = fill(SpeedyWeather.DEFAULT_GRID, nruns)
    nlat::Vector{Int} = fill(0, nruns)
    dynamics::Vector{Bool} = fill(true, nruns)
    physics::Vector{Bool} = fill(true, nruns)
    spectral_transform::Vector{Symbol} = fill(:default, nruns)
    # Per-run extra kwargs to splat into the model constructor. Values that
    # are `Type` or `Function` are resolved against `spectral_grid` first;
    # Use this to vary arbitrary model components, e.g.
    # `model_kwargs = [(convection = NoConvection,), (convection = BettsMillerConvection,)]`.
    model_kwargs::Vector = fill(NamedTuple(), nruns)
    # median SYPD of the timed runs, and of the slowest and the fastest run
    SYPD::Vector{Float64} = fill(0.0, nruns)
    SYPD_min::Vector{Float64} = fill(0.0, nruns)
    SYPD_max::Vector{Float64} = fill(0.0, nruns)
    # number of timed runs, see `repeats` and `max_repeat_time`
    nrepeats::Vector{Int} = fill(0, nruns)
    Δt::Vector{Float64} = fill(0.0, nruns)
    memory::Vector{Int} = fill(0, nruns)
    # one-off compilation of the time step (Reactant only, NaN otherwise), excluded from SYPD
    compile_time::Vector{Float64} = fill(NaN, nruns)
    # "ok" or "failed: <error type>"
    status::Vector{String} = fill("", nruns)
    architecture::Any = SpeedyWeather.CPU()
    # Scale factor on `n_timesteps` for the timed run; >1 lengthens every run for
    # more robust (e.g. publication-ready) timings. Set by the driver for :benchmark201 only.
    timestep_multiplier::Real = 1
    # Time every configuration up to `repeats` times from the same initial state, but no further
    # run once all timed runs together would take longer than `max_repeat_time` [s]
    repeats::Int = 5
    max_repeat_time::Float64 = 30
end

default_nlayers(::Type{<:Barotropic}) = 1
default_nlayers(::Type{<:ShallowWater}) = 1
default_nlayers(::Type{<:PrimitiveEquation}) = 8
default_nlayers(models) = [default_nlayers(model) for model in models]

# Number of timesteps for a timed run, chosen so every simulation takes a few
# seconds. The step count is clamped between a floor and a ceiling:
# `multiplier` scales the (clamped) result — pass e.g. `10` for longer,
# publication-ready runs with more robust timings; applied after the clamp so it
# lengthens even the floored/ceilinged configs.
n_timesteps(truncation, nlayers, multiplier = 1) =
    round(Int, multiplier * clamp(round(Int, 4.0e9 / truncation^3 / nlayers^2), 50, 1200))

# Number of steps the time loop actually takes for `period` (`clock.n_steps`, as set by
# `initialize!(simulation; period)`), incl. the spin-up steps of the time stepper. With the timed
# `period = Δt * (n_timesteps + 1)` usually `n_timesteps + 2` (1 Euler spin-up step of Leapfrog), one
# more if Δt isn't whole seconds and rounding the period to seconds rounds up (on Reactant, without
# output Δt isn't adjusted to the output interval, e.g. at T43, T342). All are timed, the SYPD is
# computed from these.
function n_steps_run(model, period)
    clock = SpeedyWeather.Clock()
    initialize!(clock, model.time_stepping, model.planet, period)
    return Int(clock.n_steps)
end

# Architecture-specific hooks. The methods here are for CPU and GPU, `reactant_runner.jl`
# adds the ones for `ReactantDevice` (it needs Reactant's `@compile` macro so it is only
# included by the driver once Reactant is loaded).

"""Return `(; time_stepping, compile_time)` where `time_stepping(simulation)` runs the main
time loop of an initialized `simulation` for the timed `period`. `compile_time` is the time it took
to compile the time loop (Reactant only, NaN otherwise), which is not included in the timed run."""
build_runner(arch, simulation, period) = (time_stepping = SpeedyWeather.time_stepping!, compile_time = NaN)

"""Wait for the device to have finished all work on `simulation`."""
synchronize_device(arch, simulation) = synchronize(arch)

"""Memory footprint of `simulation` in bytes."""
memory_footprint(arch, simulation) = Base.summarysize(simulation)

"""Return a function without arguments that calls `f(args...)` and waits for its result."""
compiled_call(arch, f, args...) = () -> (f(args...); synchronize(arch))

"""Adapt `suite` to what `arch` supports, by default unchanged."""
adapt_to_architecture!(suite, arch) = suite

"""Reset `simulation` to the initial state `initial_variables` for a timed run and return the
simulation to run. By default copied into the variables of `simulation`, so that every timed run
steps the same buffers as the warm-up and replays the GPU graphs captured there (they are per
buffer, see `graph_key` in SpeedyTransforms) instead of capturing new ones in the timed region."""
function reset_simulation!(arch, simulation, initial_variables, model)
    copy!(simulation.variables, initial_variables)
    return simulation
end

# Run every configuration of a suite. Each run lives in its own `run_benchmark!`
# call so that its model and simulation go out of scope on return; `free_memory!`
# (defined by the driver script) then collects them and hands device memory back
# before the next run. This happens outside the timed region, so it does not
# affect the measured SYPD but prevents memory piling up across many runs.
# A configuration that errors is recorded as failed without losing the rest of the suite.
# One with a `status` already, set by the driver for what the architecture doesn't support
# (e.g. Float64 on Metal), is skipped and keeps that status.
function run_benchmark_suite!(suite::AbstractBenchmarkSuite)
    for i in 1:suite.nruns
        if hasproperty(suite, :status) && !isempty(suite.status[i])
            @info "Skipping benchmark $i of $(suite.title): $(suite.status[i])"
            continue
        end
        try
            run_benchmark!(suite, i)
        catch err
            @warn "Benchmark $i of $(suite.title) failed on $(typeof(suite.architecture)); recording as failed" exception = (err, catch_backtrace())
            hasproperty(suite, :status) && (suite.status[i] = "failed: $(nameof(typeof(err)))")
        end
        free_memory!()
    end
    return suite
end

function run_benchmark!(suite::BenchmarkSuite, i::Integer)
    # unpack
    Model = suite.model[i]
    NF = suite.NF[i]
    truncation = suite.truncation[i]
    nlayers = suite.nlayers[i]
    Grid = suite.Grid[i]
    dynamics = suite.dynamics[i]
    physics = suite.physics[i]
    transform_kind = suite.spectral_transform[i]
    architecture = suite.architecture

    spectral_grid = SpectralGrid(; NF, truncation, Grid, nlayers, architecture)
    suite.nlat[i] = spectral_grid.nlat

    extra_components = resolve_model_kwargs(suite.model_kwargs[i], spectral_grid)
    model = build_model(Model, spectral_grid, architecture; transform_kind, extra_components..., feedback = Feedback(verbose = true))
    if Model <: PrimitiveEquation
        model.dynamics = dynamics
        model.dynamics_only = !physics
    else
        suite.dynamics[i] = true
        suite.physics[i] = false
    end

    simulation = initialize!(model)
    suite.memory[i] = memory_footprint(architecture, simulation)

    nsteps = n_timesteps(truncation, nlayers, suite.timestep_multiplier)
    period = Second(round(Int, model.time_stepping.Δt * (nsteps + 1)))
    nsteps_run = n_steps_run(model, period)     # steps of every timed run, incl. spin-up, for the SYPD

    # Warm up before timing: run a few steps so JIT compilation of the time
    # loop happens here rather than inside the timed run below. The timed run
    # is short (≈nsteps steps, only a few seconds), so without this its
    # wallclock is dominated by one-off compilation, which massively
    # understates SYPD and can even invert the ordering between models
    # (e.g. the model compiled first appears slower than a later one that
    # reuses the compiled code). A handful of steps covers the Euler spin-up
    # plus regular steps. We then re-initialize so the timed run starts from
    # the same initial state as before. With Reactant the whole time loop of the
    # timed run is compiled once here (timed separately), warmed up and reused.
    # (`initialize!`, `time_stepping!`, `finalize!` is what `run!` does)
    warmup_period = Second(round(Int, model.time_stepping.Δt * 3))
    initialize!(simulation; period = warmup_period)
    runner = build_runner(architecture, simulation, period)
    runner.time_stepping(simulation)
    finalize!(simulation)
    synchronize_device(architecture, simulation)

    # Time only the main time loop with a real wallclock timer rather than the
    # progress meter's `tlast - tinit`. The latter measures only up to the *last*
    # progress-bar render (every `feedback_dt = 0.1 s`). Repeat from the same initial
    # state, see `repeats` and `max_repeat_time`: initialize once (with Reactant every
    # `initialize!` compiles again, `@maybe_jit`) and reset the warm-up's simulation to it
    # before each run, see `reset_simulation!`. Every timed run thereby steps the warm-up's
    # buffers and replays the GPU graphs of the Fourier transform captured in the warm-up;
    # none are captured in the timed region. (Graphs are per buffer: new variables would
    # capture new ones. Hence also no `free_memory!` here, it would clear the graphs.)
    # After the reset the GC frees what is left over, e.g. Reactant's previous copy, outside
    # the timed region.
    initial_simulation = initialize!(model)
    initialize!(initial_simulation; period)
    initial_variables = initial_simulation.variables    # initial state of every timed run, never run itself
    initial_simulation = nothing
    times = Float64[]
    while length(times) < suite.repeats && (isempty(times) || sum(times) + times[end] <= suite.max_repeat_time)
        simulation = reset_simulation!(architecture, simulation, initial_variables, model)
        GC.gc()
        initialize!(model.feedback, simulation.variables, model)    # a new progress meter
        synchronize_device(architecture, simulation)
        t0 = time()
        runner.time_stepping(simulation)
        synchronize_device(architecture, simulation)
        push!(times, time() - t0)
        finalize!(simulation)   # but keep `simulation`, the next run resets it
    end
    @info @sprintf("     %d timed run(s) of %d steps: median %.3g s, fastest %.3g s, slowest %.3g s", length(times), nsteps_run, median(times), minimum(times), maximum(times))

    # all steps of a timed run take a full step's work, incl. the spin-up step
    sypd(time_elapsed) = model.time_stepping.Δt * nsteps_run / (time_elapsed * 365.25)

    suite.Δt[i] = model.time_stepping.Δt
    suite.SYPD[i] = sypd(median(times))
    suite.SYPD_min[i] = sypd(maximum(times))
    suite.SYPD_max[i] = sypd(minimum(times))
    suite.nrepeats[i] = length(times)
    suite.compile_time[i] = runner.compile_time
    suite.status[i] = "ok"
    return
end

function write_results(md, suite::BenchmarkSuite)

    write(md, "\n### $(suite.title)\n\n")
    isempty(suite.description) || write(md, "$(suite.description)\n\n")

    print_NF = any(suite.NF .!= suite.NF[1])
    print_Grid = any(suite.Grid .!= suite.Grid[1])
    print_nlat = any(suite.nlat .!= suite.nlat[1])
    print_dynamics = any(suite.dynamics .!= suite.dynamics[1])
    print_physics = any(suite.physics .!= suite.physics[1])
    print_transform = any(suite.spectral_transform .!= suite.spectral_transform[1])
    print_compile = any(isfinite, suite.compile_time)

    column_header = "| Model "
    column_header *= print_NF ? "| NF " : ""
    column_header *= "| truncation "
    column_header *= "| L "
    column_header *= print_Grid ? "| Grid " : ""
    column_header *= print_nlat ? "| Rings " : ""
    column_header *= print_dynamics ? "| Dynamics " : ""
    column_header *= print_physics ? "| Physics " : ""
    column_header *= print_transform ? "| Transform " : ""
    column_header *= "| Δt | SYPD | Runs | Spread | Memory"
    column_header *= print_compile ? "| Compile |" : "|"

    ncolumns = length(findall('|', column_header)) - 1
    second_row = repeat("| --- ", ncolumns) * "|"

    write(md, "$column_header\n")
    write(md, "$second_row\n")

    for i in 1:suite.nruns

        row = "| $(suite.model[i]) "
        row *= print_NF ? "| $(suite.NF[i]) " : ""
        row *= "| $(suite.truncation[i]) "
        row *= "| $(suite.nlayers[i]) "
        row *= print_Grid ? "| $(suite.Grid[i]) " : ""
        row *= print_nlat ? "| $(suite.nlat[i]) " : ""
        row *= print_dynamics ? "| $(suite.dynamics[i]) " : ""
        row *= print_physics ? "| $(suite.physics[i]) " : ""
        row *= print_transform ? "| $(suite.spectral_transform[i]) " : ""

        Δt = round(Int, suite.Δt[i])
        ok = suite.status[i] == "ok"
        SYPD = ok ? format_sypd(suite.SYPD[i]) : suite.status[i]
        runs = ok ? suite.nrepeats[i] : "—"
        spread = ok ? format_spread(suite.SYPD[i], suite.SYPD_min[i], suite.SYPD_max[i], suite.nrepeats[i]) : "—"
        memory = prettymemory(suite.memory[i])
        row *= "| $Δt | $SYPD | $runs | $spread | $memory "
        row *= print_compile ? "| $(format_seconds(suite.compile_time[i])) |" : "|"

        write(md, "$row\n")
    end
    return
end

abstract type AbstractBenchmarkSuiteTimed <: AbstractBenchmarkSuite end

@kwdef mutable struct BenchmarkSuiteDynamics <: AbstractBenchmarkSuiteTimed
    title::String
    nruns::Int = 1
    model::Vector = fill(PrimitiveWetModel, nruns)
    NF::Vector = fill(SpeedyWeather.DEFAULT_NF, nruns)
    truncation::Vector{Int} = fill(SpeedyWeather.DEFAULT_TRUNCATION, nruns)
    nlayers::Vector{Int} = default_nlayers(model)
    Grid::Vector = fill(SpeedyWeather.DEFAULT_GRID, nruns)
    nlat::Vector{Int} = fill(0, nruns)
    function_names::Vector{String} = default_function_names()
    time::Vector{Vector{Float64}} = [fill(0.0, length(function_names)) for i in 1:nruns]
    memory::Vector{Vector{Int}} = [fill(0, length(function_names)) for i in 1:nruns]
    allocs::Vector{Vector{Int}} = [fill(0, length(function_names)) for i in 1:nruns]
    architecture::Any = SpeedyWeather.CPU()
end

default_function_names() = ["pressure_gradient_flux!", "linear_virtual_temperature!", "geopotential!", "vertical_integration!", "surface_pressure_tendency!", "vertical_velocity!", "linear_pressure_gradient!", "vertical_advection!", "vordiv_tendencies!", "temperature_tendency!", "humidity_tendency!", "bernoulli_potential!"]

function add_results!(suite::AbstractBenchmarkSuiteTimed, trial::BenchmarkTools.Trial, i_run::Integer, i_func::Integer)

    t = minimum(trial)

    suite.time[i_run][i_func] = t.time
    suite.memory[i_run][i_func] = t.memory
    suite.allocs[i_run][i_func] = t.allocs

    return suite
end

# Individual dynamics functions that are not GPU-compatible: each performs a
# per-variable grid→spectral `transform!` on a non-contiguous view (a layer slot of
# the fused tendency buffer). `reinterpret` of such a strided view is not GPU-safe,
# so cuFFT falls back to a CPU path (StackOverflowError) and the Legendre kernel
# fails to compile (InvalidIRError). The full model avoids this by transforming the
# contiguous fused parent in one batched call, so these are only exercised in
# isolation here. Skip them on GPU rather than triggering the (slow) failures.
const GPU_INCOMPATIBLE_FUNCTIONS = Set(
    [
        "surface_pressure_tendency!",
        "vordiv_tendencies!",
        "temperature_tendency!",
        "humidity_tendency!",
        "bernoulli_potential!",
    ]
)

# Run one @benchmark and store the result. On failure, record NaN/0 and warn
# (used for individual function benchmarks that may not be GPU-compatible).
function safe_benchmark!(f, suite::AbstractBenchmarkSuiteTimed, i_run::Integer, i_func::Integer)
    name = suite.function_names[i_func]
    if suite.architecture isa SpeedyWeather.GPU && name in GPU_INCOMPATIBLE_FUNCTIONS
        @info "Skipping $name on $(suite.architecture): not GPU-compatible in isolation; recording N/A"
        suite.time[i_run][i_func] = NaN
        suite.memory[i_run][i_func] = 0
        suite.allocs[i_run][i_func] = 0
        return suite
    end
    try
        trial = f()
        add_results!(suite, trial, i_run, i_func)
    catch err
        @warn "Benchmark for $name failed on $(suite.architecture); recording N/A" exception = (err, catch_backtrace())
        suite.time[i_run][i_func] = NaN
        suite.memory[i_run][i_func] = 0
        suite.allocs[i_run][i_func] = 0
    end
    return suite
end

function run_benchmark!(suite::BenchmarkSuiteDynamics, i::Integer)
    arch = suite.architecture

    Model = suite.model[i]
    NF = suite.NF[i]
    truncation = suite.truncation[i]
    nlayers = suite.nlayers[i]
    Grid = suite.Grid[i]

    spectral_grid = SpectralGrid(; NF, truncation, Grid, nlayers, architecture = arch)
    suite.nlat[i] = spectral_grid.nlat

    model = build_model(Model, spectral_grid, arch)

    simulation = initialize!(model)

    vars, model = SpeedyWeather.unpack(simulation)
    (; orography, geometry, spectral_transform, geopotential, atmosphere, implicit, time_stepping) = model

    # Each function is called through `compiled_call` (compiled once with Reactant), each
    # benchmark sample also waits for the device. Building the call inside `safe_benchmark!`
    # records a function that fails to compile as N/A too.
    functions_and_arguments = (
        (SpeedyWeather.pressure_gradient_flux!, (vars, time_stepping)),
        (SpeedyWeather.linear_virtual_temperature!, (vars, model)),
        (SpeedyWeather.geopotential!, (vars, geopotential, orography)),
        (SpeedyWeather.vertical_integration!, (vars, geometry, time_stepping)),
        (SpeedyWeather.surface_pressure_tendency!, (vars, spectral_transform, time_stepping)),
        (SpeedyWeather.vertical_velocity!, (vars, geometry, time_stepping)),
        (SpeedyWeather.linear_pressure_gradient!, (vars, atmosphere, implicit, time_stepping)),
        (SpeedyWeather.vertical_advection!, (vars, model)),
        (SpeedyWeather.vordiv_tendencies!, (vars, model)),
        (SpeedyWeather.temperature_tendency!, (vars, model)),
        (SpeedyWeather.humidity_tendency!, (vars, model)),
        (SpeedyWeather.bernoulli_potential!, (vars, spectral_transform, time_stepping)),
    )
    @assert [string(nameof(f)) for (f, _) in functions_and_arguments] == suite.function_names

    for (i_func, (f, args)) in enumerate(functions_and_arguments)
        safe_benchmark!(suite, i, i_func) do
            call = compiled_call(arch, f, args...)
            @benchmark $call()
        end
    end

    return suite
end

function write_results(md, suite::AbstractBenchmarkSuiteTimed)

    write(md, "\n### $(suite.title)\n\n")

    for i_run in 1:suite.nruns

        title_row = "$(suite.model[i_run]) "
        title_row *= "| $(suite.NF[i_run]) "
        title_row *= "| T$(suite.truncation[i_run]) "
        title_row *= "L$(suite.nlayers[i_run]) "
        title_row *= "| $(suite.Grid[i_run]) "
        title_row *= "| $(suite.nlat[i_run]) Rings"

        write(md, "\n#### $title_row\n\n")

        column_header = "| Function "
        column_header *= "| Time "
        column_header *= "| Memory "
        column_header *= "| Allocations |"

        ncolumns = length(findall('|', column_header)) - 1
        second_row = repeat("| --- ", ncolumns) * "|"

        write(md, "$column_header\n")
        write(md, "$second_row\n")

        for i_func in 1:length(suite.function_names)

            t = suite.time[i_run][i_func]
            if isnan(t)
                time = "N/A"
                memory = "N/A"
                allocs = "N/A"
            else
                time = BenchmarkTools.prettytime(t)
                memory = BenchmarkTools.prettymemory(suite.memory[i_run][i_func])
                allocs = string(suite.allocs[i_run][i_func])
            end

            row = "| $(suite.function_names[i_func]) "
            row *= "| $time"
            row *= "| $memory"
            row *= "| $allocs |"

            write(md, "$row\n")
        end
    end
    return
end
