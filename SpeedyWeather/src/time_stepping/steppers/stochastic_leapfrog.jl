export StochasticLeapfrog, AdditiveWienerNoise

"""Supertype of the noise components of stochastic time steppers, see `StochasticLeapfrog`."""
abstract type AbstractWienerNoise <: AbstractModelComponent end

"""Additive Wiener noise for the prognostic variables, used by `StochasticLeapfrog`.
Every prognostic variable `x` and every layer is forced by its own independent white noise
pattern `ξ` (independent in time, between variables and between layers) with a fixed
amplitude `a`, so the noise is additive: its covariance depends neither on the state nor on
any parameter of the model. With time step `Δt` the tendency receives `a/√Δt ⋅ taper ⋅ ξ`,
so the state receives an increment ∝ `a√Δt ⋅ ξ` per step, a Wiener increment whose variance
is independent of the time step.

- vorticity, divergence (and η in the ShallowWaterModel) are forced in spectral space by a
  `SpectralWhiteNoiseProcess` (spectrally flat by default), so the noise covariance is diagonal
  in the spectral coefficients.
- temperature, humidity and surface pressure are forced in grid space with the grid
  transform of the same kind of spectral pattern. The pattern is band-limited, so the forcing
  arrives in the spectral tendencies with the spectrum of the pattern.
- the sea surface temperature (if the ocean has a tendency for it) is forced by a
  `GridWhiteNoiseProcess`, independently at every grid point.

The amplitudes multiply a pattern of standard deviation `standard_deviation` (1/3 by default,
as for `SpectralAR1Process`). The defaults are the amplitudes used for T32 L8 PrimitiveWet
models in a 1500-year stochastic control run: climate-neutral, stable, humidity bounded by
clipping rather than by the climate (see the docstring of `StochasticLeapfrog`). They are not
universal, check them for other resolutions, time steps or models.
Every variable of the model gets a pattern slot (also with amplitude zero), so switching one
amplitude on or off does not change the noise realisation of the others. The noise is drawn on
the model's architecture (see `architecture_rng`); `draw_on_cpu = true` draws it on the CPU also
on GPU (slower, but the same numbers as on CPU).

Fields are $(TYPEDFIELDS)"""
struct AdditiveWienerNoise{NF, F, VectorType, SP, GP} <: AbstractWienerNoise
    "[OPTION] Noise amplitude for vorticity [1/s^(3/2)]"
    vorticity::NF

    "[OPTION] Noise amplitude for divergence [1/s^(3/2)]"
    divergence::NF

    "[OPTION] Noise amplitude for temperature [K/s^(1/2)]"
    temperature::NF

    "[OPTION] Noise amplitude for specific humidity [kg/kg/s^(1/2)]"
    humidity::NF

    "[OPTION] Noise amplitude for the logarithm of surface pressure [1/s^(1/2)], or for the interface displacement η [m/s^(1/2)] in the ShallowWaterModel"
    pressure::NF

    "[OPTION] Noise amplitude for the sea surface temperature [K/s^(1/2)]"
    sea_surface_temperature::NF

    "[OPTION] Vertical tapering as a function of sigma, for vorticity, divergence, temperature and humidity"
    tapering::F

    "[DERIVED] Vertical tapering per layer"
    taper::VectorType

    "[OPTION] White noise process in spectral space for vorticity, divergence, temperature, humidity, pressure"
    spectral::SP

    "[OPTION] White noise process in grid space for the sea surface temperature"
    grid::GP
end

"""$(TypedMethodSignatures(defaults = true))
Additive Wiener noise for the prognostic variables, amplitudes per variable. `seed` seeds
both white noise processes (they draw from independent streams), `wavenumber` and
`standard_deviation` are passed on to the `SpectralWhiteNoiseProcess`, `draw_on_cpu` to both."""
function AdditiveWienerNoise(
        SG::SpectralGrid;
        vorticity = 3.0e-8,
        divergence = 3.0e-8,
        temperature = 3.0e-3,
        humidity = 3.0e-8,
        pressure = 0,
        sea_surface_temperature = 1.0e-3,
        tapering = σ -> 1,
        seed = 0,
        wavenumber = 0,
        standard_deviation = 1 / 3,
        draw_on_cpu = false,
    )
    NF = SG.NF
    taper = on_architecture(SG.architecture, zeros(NF, SG.nlayers))
    spectral = SpectralWhiteNoiseProcess(SG; seed, wavenumber, standard_deviation, draw_on_cpu)
    grid = GridWhiteNoiseProcess(SG; seed, standard_deviation, draw_on_cpu)
    return AdditiveWienerNoise{NF, typeof(tapering), typeof(taper), typeof(spectral), typeof(grid)}(
        NF(vorticity), NF(divergence), NF(temperature), NF(humidity), NF(pressure),
        NF(sea_surface_temperature), tapering, taper, spectral, grid,
    )
end

"""$(TYPEDSIGNATURES)
Slots of the spectral noise pattern per variable: `vorticity`, `divergence`, `temperature`,
`humidity` are offsets (slot = offset + layer), `pressure` is the slot itself, `nothing` if the
model does not have that variable, `n` is the total number of slots."""
function wiener_slots(model::AbstractModel)
    nlayers = model.spectral_grid.nlayers
    return wiener_slots(model, nlayers)
end

wiener_slots(::Barotropic, nlayers) =
    (vorticity = 0, divergence = nothing, temperature = nothing, humidity = nothing, pressure = nothing, n = nlayers)
wiener_slots(::ShallowWater, nlayers) =
    (vorticity = 0, divergence = nlayers, temperature = nothing, humidity = nothing, pressure = 2nlayers + 1, n = 2nlayers + 1)
wiener_slots(::PrimitiveDry, nlayers) =
    (vorticity = 0, divergence = nlayers, temperature = 2nlayers, humidity = nothing, pressure = 3nlayers + 1, n = 3nlayers + 1)
wiener_slots(::PrimitiveWet, nlayers) =
    (vorticity = 0, divergence = nlayers, temperature = 2nlayers, humidity = 3nlayers, pressure = 4nlayers + 1, n = 4nlayers + 1)

"""$(TYPEDSIGNATURES) Does the ocean of `model` have a tendency for the sea surface temperature?"""
function has_sea_surface_temperature_tendency(model::AbstractModel)
    hasproperty(model, :ocean) || return false
    return any(variables(model.ocean, model)) do var
        var isa TendencyVariable && var.name == :sea_surface_temperature && var.namespace == :ocean
    end
end

function variables(noise::AdditiveWienerNoise, model::AbstractModel)
    n = wiener_slots(model).n
    spectral = (
        PrognosticVariable(
            :wiener_pattern, Spectral3D(; n), units = "1",
            desc = "Pattern of the Wiener increment of the next time step, one slot per variable and layer",
        ),
    )
    # temperature, humidity and pressure are forced in grid space, with the grid transform of the pattern
    grid = model isa PrimitiveEquation ? (
            GridVariable(
                :wiener_pattern, Grid3D(; n), units = "1",
                desc = "Pattern of the Wiener increment of the next time step on the grid",
            ),
        ) : ()
    sst = has_sea_surface_temperature_tendency(model) ? (
            PrognosticVariable(
                :wiener_pattern_grid, Grid2D(), units = "1",
                desc = "Grid-point pattern of the Wiener increment of the sea surface temperature of the next time step",
            ),
        ) : ()
    return (spectral..., grid..., sst...)
end

function initialize!(noise::AdditiveWienerNoise, model::AbstractModel)
    initialize!(noise.spectral, model)
    initialize!(noise.grid, model)

    # with the same seed both processes would draw the same normals, derive a distinct
    # (but reproducible) stream for the grid process
    if noise.spectral.seed == noise.grid.seed && noise.grid.seed != 0
        Random.seed!(noise.grid.random_number_generator, hash(noise.grid.seed, UInt(0x9e3779b97f4a7c15)))
    end

    # vertical tapering, in a kernel as a comprehension would scalar-index the vertical coordinates on GPU
    coord = model.geometry.vertical_coordinates
    launch!(
        architecture(noise.taper), LinearWorkOrder, (get_nlayers(coord),),
        wiener_taper_kernel!, noise.taper, noise.tapering, coord,
    )

    # models with an ocean that has no sea surface temperature tendency (e.g. prescribed), models
    # without ocean (Barotropic, ShallowWater) ignore the sea surface temperature amplitude
    if noise.sea_surface_temperature != 0 && hasproperty(model, :ocean) && !has_sea_surface_temperature_tendency(model)
        @warn "AdditiveWienerNoise: sea_surface_temperature = $(noise.sea_surface_temperature) but " *
            "the ocean $(typeof(model.ocean).name.name) has no sea surface temperature tendency, " *
            "the sea surface temperature is not forced. Set sea_surface_temperature = 0 to silence."
    end
    return nothing
end

@kernel inbounds = true function wiener_taper_kernel!(taper, tapering, coordinate)
    k = @index(Global, Linear)
    taper[k] = tapering(sigma(k, coordinate))
end

"""Leapfrog time stepping with additive Wiener noise, i.e. the stochastic differential equation

    dx = F(x) dt + a dW

with the deterministic tendencies F of the model, white noise patterns of fixed amplitude `a`
(see `AdditiveWienerNoise`) and the semi-implicit leapfrog scheme of `Leapfrog` (Euler start,
Robert-Williams filter, implicit diffusion). Every time step adds `a/√Δt ⋅ ξ` to the tendencies
after the parameterizations and before the dynamical core, i.e. the noise enters exactly like a
forcing; it is radius-scaled like the tendencies, and is part of the tendency that the
semi-implicit correction, the horizontal diffusion and the leapfrog step (2Δt) act on. The
leapfrog step from `n-1` to `n+1` therefore adds `2a√Δt ⋅ ξⁿ`. Because of the two interleaved
leapfrog sequences, the time filter and the implicit terms, the diffusion this scheme
effectively injects is not the nominal `a²σ²/2` and should be measured where it matters
(e.g. 1.33× nominal for a T8 barotropic model).

The pattern `ξ` applied in a time step is drawn at the end of the previous time step and stored
in the variables (`vars.prognostic.wiener_pattern`, `vars.grid.wiener_pattern`,
`vars.prognostic.wiener_pattern_grid`), so the first time step of a new simulation is
deterministic and a simulation continued with another `run!` continues with the stored pattern.
The noise realisation is reproducible through `seed`, independent of other random processes.

`StochasticLeapfrog` is a `Leapfrog` whose `noise` field holds an `AbstractWienerNoise` (a type
alias, not a new time stepper): every method of `Leapfrog` applies, only the Wiener increment is
added. Use as `time_stepping = StochasticLeapfrog(spectral_grid; seed = 123)`, keyword arguments
other than those of `Leapfrog` are passed on to `AdditiveWienerNoise`, or equivalently
`Leapfrog(spectral_grid; noise = AdditiveWienerNoise(spectral_grid; seed = 123))`.
All amplitudes zero gives the same trajectory as `Leapfrog` without noise."""
const StochasticLeapfrog{NF, S, B, MS, O, L, N <: AbstractWienerNoise} = Leapfrog{NF, S, B, MS, O, L, N}

"""$(TypedMethodSignatures(defaults = true))
Leapfrog with additive Wiener noise. Keyword arguments of `Leapfrog` are passed on to it, all
others to `AdditiveWienerNoise` unless `noise` is provided."""
function StochasticLeapfrog(
        spectral_grid::SpectralGrid;
        Δt_at_T32 = Minute(40),
        adjust_with_output = true,
        robert_filter = 0.1,
        williams_filter = 0.53,
        ocean = EulerForward(spectral_grid; Δt_at_T32, adjust_with_output),
        land = EulerForward(spectral_grid; Δt_at_T32, adjust_with_output),
        noise = nothing,
        kwargs...,
    )
    if isnothing(noise)
        noise = AdditiveWienerNoise(spectral_grid; kwargs...)
    elseif !isempty(kwargs)
        throw(ArgumentError("Provide either `noise` or its keyword arguments $(keys(kwargs)), not both."))
    end
    return Leapfrog(spectral_grid; Δt_at_T32, adjust_with_output, robert_filter, williams_filter, ocean, land, noise)
end

"""$(TYPEDSIGNATURES)
Add the Wiener increment of this time step (drawn at the end of the previous one) to the
tendencies, then draw the pattern for the next time step."""
function wiener_increment!(vars::Variables, L::StochasticLeapfrog, model::AbstractModel)
    wiener_tendencies!(vars, L.noise, L, model)
    wiener_draw!(vars, L.noise, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Add the stored noise pattern times `amplitude/√Δt` to the tendencies."""
function wiener_tendencies!(
        vars::Variables,
        noise::AdditiveWienerNoise,
        L::AbstractTimeStepper,
        model::AbstractModel,
    )
    slots = wiener_slots(model)
    (; taper) = noise
    isdt = 1 / sqrt(L.Δt)               # Δt in [s], the Wiener increment scaling
    R = vars.prognostic.scale[]         # radius scaling of vorticity and divergence

    # VORTICITY, DIVERGENCE in spectral space. scale_tendencies! multiplies every tendency
    # by the radius R and the time step is Δt/R, but vorticity and divergence are stored
    # scaled by R, so their tendency needs R once more here
    spec = vars.prognostic.wiener_pattern
    wiener_spectral!(vars.tendencies.vorticity, spec, noise.vorticity * isdt * R, taper, slots.vorticity, L, noise)
    isnothing(slots.divergence) ||
        wiener_spectral!(vars.tendencies.divergence, spec, noise.divergence * isdt * R, taper, slots.divergence, L, noise)

    # TEMPERATURE, HUMIDITY in grid space, their spectral tendencies are overwritten by the
    # batched transform of the grid tendencies in the dynamical core
    if !isnothing(slots.temperature)
        grid = vars.grid.wiener_pattern
        wiener_grid!(vars.tendencies.grid.temperature, grid, noise.temperature * isdt, taper, slots.temperature, L, noise)
        isnothing(slots.humidity) ||
            wiener_grid!(vars.tendencies.grid.humidity, grid, noise.humidity * isdt, taper, slots.humidity, L, noise)
    end

    # PRESSURE (log surface pressure on the grid, or η in spectral space for ShallowWater), no taper
    if !isnothing(slots.pressure) && noise.pressure != 0
        wiener_pressure!(vars, noise.pressure * isdt, slots.pressure, L, noise, model)
    end

    # SEA SURFACE TEMPERATURE, unscaled: ocean variables are stepped with the unscaled time step
    if noise.sea_surface_temperature != 0 && haskey(vars.prognostic, :wiener_pattern_grid)
        dsst = get_tendency_step(vars.tendencies.ocean.sea_surface_temperature, L, noise)
        launch!(
            architecture(dsst), RingGridWorkOrder, size(dsst), wiener_surface_kernel!,
            dsst, vars.prognostic.wiener_pattern_grid, noise.sea_surface_temperature * isdt,
        )
    end
    return nothing
end

function wiener_spectral!(tendency, pattern, amplitude, taper, offset, L, noise)
    iszero(amplitude) && return nothing
    tend = get_tendency_step(tendency, L, noise)
    launch!(
        architecture(tend), SpectralWorkOrder, size(tend), wiener_layered_kernel!,
        tend, pattern, amplitude, taper, offset,
    )
    return nothing
end

function wiener_grid!(tendency, pattern, amplitude, taper, offset, L, noise)
    iszero(amplitude) && return nothing
    tend = get_tendency_step(tendency, L, noise)
    launch!(
        architecture(tend), RingGridWorkOrder, size(tend), wiener_layered_kernel!,
        tend, pattern, amplitude, taper, offset,
    )
    return nothing
end

# log surface pressure, forced in grid space with the transformed pattern
function wiener_pressure!(vars, amplitude, slot, L, noise, ::PrimitiveEquation)
    tend = get_tendency_step(vars.tendencies.grid.pressure, L, noise)
    launch!(
        architecture(tend), RingGridWorkOrder, size(tend), wiener_slot_kernel!,
        tend, vars.grid.wiener_pattern, amplitude, slot,
    )
    return nothing
end

# interface displacement η, forced in spectral space (its spectral tendency accumulates)
function wiener_pressure!(vars, amplitude, slot, L, noise, ::ShallowWater)
    tend = get_tendency_step(vars.tendencies.η, L, noise)
    launch!(
        architecture(tend), SpectralWorkOrder, size(tend), wiener_slot_kernel!,
        tend, vars.prognostic.wiener_pattern, amplitude, slot,
    )
    return nothing
end

# tendency[i, k] += a * taper[k] * pattern[i, offset + k] for every (spectral or grid) index i and layer k
@kernel inbounds = true function wiener_layered_kernel!(tend, pattern, a, taper, offset)
    I = @index(Global, NTuple)
    i = I[1]
    k = length(I) > 1 ? I[2] : 1
    t = taper[k]
    tend[i, k] += a * t * pattern[i, offset + k]
end

# single-layer variables, one slot of the pattern
@kernel inbounds = true function wiener_slot_kernel!(tend, pattern, a, slot)
    i = @index(Global, Linear)
    tend[i] += a * pattern[i, slot]
end

# grid-point pattern for surface variables
@kernel inbounds = true function wiener_surface_kernel!(tend, pattern, a)
    ij = @index(Global, Linear)
    tend[ij] += a * pattern[ij]
end

"""$(TYPEDSIGNATURES)
Draw the noise patterns for the next time step, transform the spectral one to grid space."""
function wiener_draw!(vars::Variables, noise::AdditiveWienerNoise, model::AbstractModel)
    random_process!(vars.prognostic.wiener_pattern, noise.spectral)
    haskey(vars.prognostic, :wiener_pattern_grid) &&
        random_process!(vars.prognostic.wiener_pattern_grid, noise.grid)
    if haskey(vars.grid, :wiener_pattern)
        transform!(
            vars.grid.wiener_pattern, vars.prognostic.wiener_pattern,
            vars.scratch.transform_memory, model.spectral_transform,
        )
    end
    return nothing
end
