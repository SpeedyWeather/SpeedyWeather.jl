export AdditiveWienerNoise

"""Supertype of the noise components of stochastic time steppers, see `StochasticLeapfrog`."""
abstract type AbstractWienerNoise <: AbstractModelComponent end

"""Additive Wiener noise for the prognostic variables, used by `StochasticLeapfrog`.
Every prognostic variable `x` and every layer is forced by its own independent white noise
pattern `ξ` (independent in time, between variables and between layers) with a fixed
amplitude `a`, so the noise is additive: its covariance depends neither on the state nor on
any parameter of the model. With time step `Δt` the tendency receives `a/√Δt ⋅ taper ⋅ ξ`,
so the state receives an increment ∝ `a√Δt ⋅ ξ` per step, a Wiener increment whose variance
is independent of the time step. `Δt` is the time step of the respective time stepper at the
current step, half of it for a leapfrog step, see `wiener_scale`.

The amplitudes are given per variable by the name of its tendency, e.g.
`AdditiveWienerNoise(spectral_grid; vorticity = 3e-8, η = 1e-3, soil_temperature = 1e-4)`,
variables not listed get no noise. Units are those of the variable per √s, e.g. 1/s^(3/2) for
vorticity and divergence, K/s^(1/2) for temperature. Which variables there are is inferred from
the tendencies of the model (`tendency_names`, `ocean_tendency_names`, `land_tendency_names`):

- atmospheric variables are forced with one `SpectralWhiteNoiseProcess` (spectrally flat by
  default) that has a pattern slot for every layer of every variable (variables with layers first,
  then single-layer variables, each in the order of `tendency_names`), also for zero amplitudes,
  so that switching one amplitude on or off does not change the noise realisation of the others.
  The noise is added to the spectral tendencies, so its covariance is diagonal in the spectral
  coefficients. The pattern has no global mean, the noise therefore conserves e.g. the mass
  (surface pressure). Variables with layers are tapered with `tapering(σ)`.
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

"""$(TYPEDSIGNATURES)
Additive Wiener noise for the prognostic variables, all keyword arguments other than the ones
listed are noise amplitudes by variable name, added to (or replacing) `DEFAULT_WIENER_AMPLITUDES`.
`seed` seeds both white noise processes (salted: independent streams, also of other random
processes with the same seed, e.g. `model.random_process`), `wavenumber` and
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

    spectral = PrognosticVariable(
        :wiener_pattern, Spectral3D(; n), units = "1",
        desc = "Pattern of the Wiener increment of the next time step, one slot per variable and layer",
    )
    surface = isempty(surface) ? () : (
            PrognosticVariable(
                :wiener_pattern_surface, Grid3D(; n = sum(v -> nslots(v, model), surface)), units = "1",
                desc = "Grid-point pattern of the Wiener increment of ocean and land variables of the next time step",
            ),
        )
    return (spectral, surface...)
end

# seed salts of the two white noise processes, see `salted_seed`, ASCII "WIENERSP" and "WIENERGR"
const WIENER_SPECTRAL_SALT = 0x5749454e45525350
const WIENER_GRID_SALT = 0x5749454e45524752

function initialize!(noise::AdditiveWienerNoise, model::AbstractModel)
    initialize!(noise.spectral, model)
    initialize!(noise.grid, model)

    # salt the seeds: streams independent of each other and of other random processes
    # with the same seed, e.g. model.random_process
    salt_seed!(noise.spectral, WIENER_SPECTRAL_SALT)
    salt_seed!(noise.grid, WIENER_GRID_SALT)

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


"""$(TYPEDSIGNATURES)
Add the stored noise patterns times `amplitude/√Δt` to the tendencies: of the atmospheric
variables in spectral space, radius-scaled as the dynamical core leaves them, of the ocean and
land variables in grid space, see `AdditiveWienerNoise` and `wiener_scale` for `Δt`."""
function wiener_tendencies!(
        vars::Variables,
        noise::AdditiveWienerNoise,
        L::AbstractTimeStepper,
        model::AbstractModel,
    )
    clock = vars.prognostic.clock
    wiener_scaling = wiener_scale(L, time_step(L, clock))   # 1/√Δt [1/√s]
    R = vars.prognostic.scale[]                             # radius scaling
    wiener_atmosphere!(vars, noise, L, wiener_scaling, R)
    wiener_surface!(vars, noise, L, clock)
    return nothing
end

"""$(TYPEDSIGNATURES)
Scaling `1/√Δt` [1/√s] of the Wiener increment `a/√Δt ⋅ ξ` added to the tendencies of the
variables that `stepper` steps with the time step `Δt` [s] (from `time_step`), so that the
increment of the state `Δt ⋅ a/√Δt ⋅ ξ` has the variance `a²Δt` (Euler-Maruyama)."""
wiener_scale(::AbstractTimeStepper, Δt) = 1 / sqrt(Δt)

# The drivers below are unrolled over the variable names, which are compile-time constants (they
# come from the type of the `Variables`), so that every `getfield` infers concretely; a loop over
# `tendency_names(vars)` would dispatch at runtime on every variable. Same pattern as
# `update_prognostic_atmosphere!` and `update_prognostic_namespaces!`.

#TODO: might update the function below to use the fused vairbales systems
"""$(TYPEDSIGNATURES)
Add the noise to the atmospheric tendencies. Slots of the spectral pattern: variables with layers
first, then single-layer variables, each in the order of `tendency_names`. After the dynamical core
the tendencies are scaled with the radius R (`scale_tendencies!`, the time step is Δt/R), vorticity
and divergence are stored scaled by R in addition (`scale_prognostic!`), so need R once more."""
@generated function wiener_atmosphere!(vars::Variables{Po, G, T}, noise, L, wiener_scaling, R) where {Po, G, T}
    calls = Expr[]
    for layers in (true, false), name in _tendency_names(T)
        # R: tendencies are radius-scaled (time step Δt/R); R²: vorticity, divergence are also stored × R
        a = name in (:vorticity, :divergence) ? :(amplitude(noise, $(QuoteNode(name))) * wiener_scaling * R * R) :
            :(amplitude(noise, $(QuoteNode(name))) * wiener_scaling * R)
        push!(
            calls, :(
                offset = wiener_slots!(
                    getfield(vars.tendencies, $(QuoteNode(name))), vars.prognostic.wiener_pattern,
                    $a, noise.taper, offset, Val($layers), L, noise,
                )
            )
        )
    end
    return Expr(:block, :(offset = 0), calls..., :(return nothing))
end

# add the noise to the tendency `var` if its layers match the pass `layers` (with layers first, then
# without), return the next slot of the pattern; `has_layers` of a concrete type is a constant
@inline function wiener_slots!(var, pattern, a, taper, offset, ::Val{layers}, L, noise) where {layers}
    tendency = get_tendency_step(var, L, noise)
    has_layers(tendency) == layers || return offset
    wiener!(tendency, pattern, a, layers ? taper : nothing, offset)
    return offset + nslots(tendency)
end

"""$(TYPEDSIGNATURES)
Add the noise to the ocean and land tendencies. Slots of the grid-point pattern only for variables
with a nonzero amplitude, not radius-scaled: ocean and land variables are stepped with the unscaled
time step. Without any such variable there is no `wiener_pattern_surface` and nothing to do."""
@generated function wiener_surface!(vars::Variables{Po, G, T}, noise, L, clock) where {Po, G, T}
    :wiener_pattern_surface in fieldnames(Po) || return :(return nothing)
    calls = Expr[]
    for namespace in (:ocean, :land), name in _namespace_names(T, namespace)
        push!(
            calls, :(
                offset = wiener_surface_slots!(
                    getfield(getfield(vars.tendencies, $(QuoteNode(namespace))), $(QuoteNode(name))),
                    vars.prognostic.wiener_pattern_surface, amplitude(noise, $(QuoteNode(name))),
                    offset, L, Val($(QuoteNode(namespace))), clock, noise,
                )
            )
        )
    end
    return Expr(:block, :(offset = 0), calls..., :(return nothing))
end

# the ocean or land time stepper (e.g. EulerForward) decides the tendency step and the Wiener scaling,
# with the time step of the ocean or land variables within `L`, e.g. Δt/2 on the first two Leapfrog steps
@inline function wiener_surface_slots!(var, pattern, a, offset, L, component::Val{C}, clock, noise) where {C}
    iszero(a) && return offset
    stepper = time_stepper(L, component)
    tendency = get_tendency_step(var, stepper, noise)
    wiener_scaling = wiener_scale(stepper, time_step(L, C, clock))
    wiener!(tendency, pattern, a * wiener_scaling, nothing, offset)
    return offset + nslots(tendency)
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
Draw the noise patterns for the next time step."""
function wiener_draw!(vars::Variables, noise::AdditiveWienerNoise)
    update_random!(vars.prognostic.wiener_pattern, noise.spectral)
    haskey(vars.prognostic, :wiener_pattern_surface) &&
        update_random!(vars.prognostic.wiener_pattern_surface, noise.grid)
    return nothing
end
