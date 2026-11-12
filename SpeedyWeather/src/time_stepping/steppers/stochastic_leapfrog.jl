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

The amplitudes are given per variable by the name of its tendency, e.g.
`AdditiveWienerNoise(spectral_grid; vorticity = 3e-8, η = 1e-3, soil_temperature = 1e-4)`,
variables not listed get no noise. Units are those of the variable per √s, e.g. 1/s^(3/2) for
vorticity and divergence, K/s^(1/2) for temperature. Which variables there are is inferred from
the tendencies of the model (`tendency_names`, `ocean_tendency_names`, `land_tendency_names`):

- atmospheric variables are forced with one `SpectralWhiteNoiseProcess` (spectrally flat by
  default) that has a pattern slot for every layer of every variable (variables with layers first,
  then single-layer variables, each in the order of `tendency_names`), also for zero amplitudes,
  so that switching one amplitude on or off does not change the noise realisation of the others.
  The noise is added where the dynamical core reads the tendency: in grid space (the grid
  transform of the pattern) for the variables whose grid tendencies the dynamical core transforms
  to spectral space (temperature, humidity, surface pressure in the primitive equation models),
  in spectral space otherwise (vorticity, divergence, η), so the noise covariance is diagonal
  in the spectral coefficients. Variables with layers are tapered with `tapering(σ)`.
- ocean and land variables (e.g. sea surface temperature) with a nonzero amplitude are forced by a
  `GridWhiteNoiseProcess`, independently at every grid point and layer.

The amplitudes multiply a pattern of standard deviation `standard_deviation` (1/3 by default,
as for `SpectralAR1Process`). The defaults are the amplitudes used for T32 L8 PrimitiveWet
models in a 1500-year stochastic control run: climate-neutral, stable, humidity bounded by
clipping rather than by the climate (see the docstring of `StochasticLeapfrog`). They are not
universal, check them for other resolutions, time steps or models. The noise is drawn on
the model's architecture (see `architecture_rng`); `draw_on_cpu = true` draws it on the CPU also
on GPU (slower, but the same numbers as on CPU).

Fields are $(TYPEDFIELDS)"""
struct AdditiveWienerNoise{NF, A <: NamedTuple, F, VectorType, SP, GP} <: AbstractWienerNoise
    "[OPTION] Noise amplitude per variable, by the name of its tendency, zero for all others"
    amplitudes::A

    "[OPTION] Vertical tapering as a function of sigma, for atmospheric variables with layers"
    tapering::F

    "[DERIVED] Vertical tapering per layer"
    taper::VectorType

    "[OPTION] White noise process in spectral space for the atmospheric variables"
    spectral::SP

    "[OPTION] White noise process in grid space for the ocean and land variables"
    grid::GP
end

"""Default amplitudes of `AdditiveWienerNoise`, for T32 L8 PrimitiveWet models."""
const DEFAULT_WIENER_AMPLITUDES = (
    vorticity = 3.0e-8,                 # [1/s^(3/2)]
    divergence = 3.0e-8,                # [1/s^(3/2)]
    temperature = 3.0e-3,               # [K/s^(1/2)]
    humidity = 3.0e-8,                  # [kg/kg/s^(1/2)]
    sea_surface_temperature = 1.0e-3,   # [K/s^(1/2)]
)

"""$(TypedMethodSignatures(defaults = true))
Additive Wiener noise for the prognostic variables, all keyword arguments other than the ones
listed are noise amplitudes by variable name, added to (or replacing) `DEFAULT_WIENER_AMPLITUDES`.
`seed` seeds both white noise processes (they draw from independent streams), `wavenumber` and
`standard_deviation` are passed on to the `SpectralWhiteNoiseProcess`, `draw_on_cpu` to both."""
function AdditiveWienerNoise(
        SG::SpectralGrid;
        tapering = σ -> 1,
        seed = 0,
        wavenumber = 0,
        standard_deviation = 1 / 3,
        draw_on_cpu = false,
        amplitudes...,
    )
    NF = SG.NF
    amplitudes = map(NF, merge(DEFAULT_WIENER_AMPLITUDES, amplitudes))
    taper = on_architecture(SG.architecture, zeros(NF, SG.nlayers))
    spectral = SpectralWhiteNoiseProcess(SG; seed, wavenumber, standard_deviation, draw_on_cpu)
    grid = GridWhiteNoiseProcess(SG; seed, standard_deviation, draw_on_cpu)
    return AdditiveWienerNoise{NF, typeof(amplitudes), typeof(tapering), typeof(taper), typeof(spectral), typeof(grid)}(
        amplitudes, tapering, taper, spectral, grid,
    )
end

"""$(TYPEDSIGNATURES) Noise amplitude of the variable `name`, zero if not given."""
amplitude(noise::AdditiveWienerNoise{NF}, name::Symbol) where {NF} = get(noise.amplitudes, name, zero(NF))

"""$(TYPEDSIGNATURES)
The tendency variables of `model`, i.e. the definitions behind `tendency_names` and co., as they
are needed before the `Variables` exist. Collected from the model and its components except the
time stepper (its variables include the noise patterns, which are defined from these)."""
function tendency_variables(model::AbstractModel)
    components = (getproperty(model, name) for name in propertynames(model) if name != :time_stepping)
    all = (variables(model)..., Iterators.flatten(variables(c, model) for c in components)...)
    return unique(v -> (v.namespace, v.name), filter(v -> v isa TendencyVariable, collect(all)))
end

# number of pattern slots of a tendency (of one time step): its atmosphere, land or ocean layers, else 1
nslots(tendency) = size(tendency.data, 2)
nslots(v::TendencyVariable, model::AbstractModel) = nslots(get_step(allocate(v, model), 1))
has_layers(tendency) = ndims(tendency.data) == 2

function variables(noise::AdditiveWienerNoise, model::AbstractModel)
    tendencies = tendency_variables(model)
    atmosphere = filter(v -> v.namespace == Symbol(), tendencies)
    surface = filter(v -> v.namespace in (:ocean, :land) && !iszero(amplitude(noise, v.name)), tendencies)
    n = sum(v -> nslots(v, model), atmosphere)

    # forced in grid space: variables whose grid tendencies the dynamical core transforms to spectral
    on_grid = any(tendencies) do v
        v.namespace == :grid && v.fuse == :grid_tendencies && any(a -> a.name == v.name, atmosphere)
    end

    spectral = PrognosticVariable(
        :wiener_pattern, Spectral3D(; n), units = "1",
        desc = "Pattern of the Wiener increment of the next time step, one slot per variable and layer",
    )
    grid = on_grid ? (
            GridVariable(
                :wiener_pattern, Grid3D(; n), units = "1",
                desc = "Pattern of the Wiener increment of the next time step on the grid",
            ),
        ) : ()
    surface = isempty(surface) ? () : (
            PrognosticVariable(
                :wiener_pattern_surface, Grid3D(; n = sum(v -> nslots(v, model), surface)), units = "1",
                desc = "Grid-point pattern of the Wiener increment of ocean and land variables of the next time step",
            ),
        )
    return (spectral, grid..., surface...)
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

    # amplitudes for variables this model does not time step are ignored, warn unless a default
    names = [v.name for v in tendency_variables(model) if v.namespace in (Symbol(), :ocean, :land)]
    ignored = filter(name -> !(name in names || haskey(DEFAULT_WIENER_AMPLITUDES, name)), keys(noise.amplitudes))
    isempty(ignored) || @warn "AdditiveWienerNoise: the model has no tendency for $ignored, amplitude ignored."
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
`vars.prognostic.wiener_pattern_surface`), so the first time step of a new simulation is
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
function random_process!(vars::Variables, L::StochasticLeapfrog, model::AbstractModel)
    wiener_tendencies!(vars, L.noise, L, model)
    wiener_draw!(vars, L.noise, model)
    return nothing
end

"""$(TYPEDSIGNATURES)
Add the stored noise patterns times `amplitude/√Δt` to the tendencies, see `AdditiveWienerNoise`
for which variable is forced where."""
function wiener_tendencies!(
        vars::Variables,
        noise::AdditiveWienerNoise,
        L::AbstractTimeStepper,
        model::AbstractModel,
    )
    isdt = 1 / sqrt(L.Δt)               # Δt in [s], the Wiener increment scaling
    R = vars.prognostic.scale[]         # radius scaling of vorticity and divergence

    # ATMOSPHERE, slots of the spectral pattern: variables with layers first, then single-layer ones
    on_grid = haskey(vars.fused, :grid_tendencies) ? keys(vars.fused.grid_tendencies.slot_map) : ()
    offset = 0
    for layers in (true, false), name in tendency_names(vars)
        tendency = get_tendency_step(getfield(vars.tendencies, name), L, noise)
        has_layers(tendency) == layers || continue

        # scale_tendencies! multiplies every tendency by the radius R and the time step is Δt/R,
        # vorticity and divergence are stored scaled by R (scale_prognostic!), so need R once more
        a = amplitude(noise, name) * isdt * (name in (:vorticity, :divergence) ? R : one(R))
        taper = layers ? noise.taper : nothing

        # the grid tendency if the dynamical core transforms it to (and overwrites) the spectral one
        if name in on_grid
            grid_tendency = get_tendency_step(getfield(vars.tendencies.grid, name), L, noise)
            wiener!(grid_tendency, vars.grid.wiener_pattern, a, taper, offset)
        else
            wiener!(tendency, vars.prognostic.wiener_pattern, a, taper, offset)
        end
        offset += nslots(tendency)
    end

    # OCEAN, LAND, slots of the grid-point pattern for variables with a nonzero amplitude, not
    # radius-scaled: ocean and land variables are stepped with the unscaled time step
    offset = 0
    for (namespace, names) in ((:ocean, ocean_tendency_names(vars)), (:land, land_tendency_names(vars))), name in names
        a = amplitude(noise, name)
        iszero(a) && continue
        tendency = get_tendency_step(getfield(getfield(vars.tendencies, namespace), name), L, noise)
        wiener!(tendency, vars.prognostic.wiener_pattern_surface, a * isdt, nothing, offset)
        offset += nslots(tendency)
    end
    return nothing
end

# tendency[i, k] += a * taper[k] * pattern[i, offset + k] for every spectral or grid index i, layer k
function wiener!(tendency, pattern, a, taper, offset)
    iszero(a) && return nothing
    workorder = tendency isa LowerTriangularArray ? SpectralWorkOrder : RingGridWorkOrder
    launch!(architecture(tendency), workorder, size(tendency), wiener_kernel!, tendency, pattern, a, taper, offset)
    return nothing
end

@kernel inbounds = true function wiener_kernel!(tendency, pattern, a, taper, offset)
    I = @index(Global, NTuple)
    i = I[1]
    k = length(I) > 1 ? I[2] : 1
    t = isnothing(taper) ? one(a) : taper[k]
    tendency[I...] += a * t * pattern[i, offset + k]
end

"""$(TYPEDSIGNATURES)
Draw the noise patterns for the next time step, transform the spectral one to grid space."""
function wiener_draw!(vars::Variables, noise::AdditiveWienerNoise, model::AbstractModel)
    update_random!(vars.prognostic.wiener_pattern, noise.spectral)
    haskey(vars.prognostic, :wiener_pattern_surface) &&
        update_random!(vars.prognostic.wiener_pattern_surface, noise.grid)
    if haskey(vars.grid, :wiener_pattern)
        transform!(
            vars.grid.wiener_pattern, vars.prognostic.wiener_pattern,
            vars.scratch.transform_memory, model.spectral_transform,
        )
    end
    return nothing
end
