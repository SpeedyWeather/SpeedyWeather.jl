"""$(TYPEDSIGNATURES)
Implementation of random processes for stochastic forcings and parameterizations in SpeedyWeather.
Any random process <: AbstractRandomProcess has to implement the following interface:
 - `variables(process)` defines the variables of `model.random_process`, by convention the
   spectral `vars.prognostic.random_pattern` and its grid counterpart `vars.grid.random_pattern`
 - `initialize!(process, model)` precomputes and (re)seeds the random number generator
 - `update_random!(pattern, process)` advances the process with `pattern` as its state, i.e.
   draws the next pattern into the array `pattern` (for processes with memory like the AR1
   the previous pattern is read from it). Any model component can own a process and its
   own pattern this way, e.g. the `StochasticLeapfrog` time stepper
 - `update_random!(vars, process)` advances `model.random_process` at the end of every time
   step, defaults to the spectral `vars.prognostic.random_pattern`
 - `transform!(vars, process, spectral_transform)` makes `vars.grid.random_pattern` from it.

Random numbers are drawn on the architecture of the model: the processes' `random_number_generator`
is a `Random.Xoshiro` on CPU and the counter-based `GPUArrays.RNG` on GPU, which draws in a kernel
on the device. `draw_on_cpu = true` uses a `Xoshiro` also on GPU: drawn on the CPU and copied over,
slower, but the same numbers on every architecture."""
abstract type AbstractRandomProcess <: AbstractModelComponent end

"""$(TYPEDSIGNATURES)
Advance the random process of the model, `model.random_process`, in its default variable
`vars.prognostic.random_pattern` by one time step."""
update_random!(vars::Variables, process::AbstractRandomProcess) =
    update_random!(vars.prognostic.random_pattern, process)

"""$(TYPEDSIGNATURES)
General transform for `random processes <: AbstractRandomProcess`.
Takes the spectral `random_pattern` in the prognostic variables
and transforms it to spectral space in `diagn.grid.random_pattern`."""
function SpeedyTransforms.transform!(
        vars::Variables,
        random_process::AbstractRandomProcess,
        spectral_transform::AbstractSpectralTransform,
    )
    pattern = vars.prognostic.random_pattern
    pattern_grid = vars.grid.random_pattern
    scratch_memory = vars.scratch.transform_memory
    transform!(pattern_grid, pattern, scratch_memory, spectral_transform)

    if :clamp in fieldnames(typeof(random_process))
        lo, hi = random_process.clamp
        @. pattern_grid = clamp(pattern_grid, lo, hi)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
`random_process=nothing` does not need to transform any random pattern from
spectral to grid space."""
SpeedyTransforms.transform!(::Variables, ::Nothing, ::AbstractSpectralTransform) = nothing
update_random!(::Variables, process::Nothing) = nothing

"""$(TYPEDSIGNATURES)
Random number generator for drawing on `architecture`: `Random.Xoshiro` on CPU, the counter-based
`GPUArrays.RNG` on GPU (random numbers generated in a kernel on the device). With `draw_on_cpu = true`
always a `Xoshiro`, its numbers are drawn on the CPU and copied to the device."""
architecture_rng(architecture, seed; draw_on_cpu::Bool = false) =
    (draw_on_cpu || architecture isa CPU) ? Random.Xoshiro(seed) : GPUArrays.RNG{array_type(architecture)}(seed)

"""$(TYPEDSIGNATURES)
Fill `r` in place with complex normals drawn from `rng`, unit variance in real and imaginary part."""
function complex_normals!(r::LowerTriangularArray, rng::Random.AbstractRNG)
    Random.randn!(rng, r.data)

    # Julia's randn(Complex) has total variance E|z|² = 1, i.e. variance 1/2 in real and imaginary part
    # each. The noise factors of the random processes are normalised for unit variance in each part,
    # so multiply by √2 to get the grid-point standard deviation the process specifies.
    # Written as 2/√2 to keep the random
    # numbers of earlier versions bit for bit; drawing real normals directly into `reinterpret(real, r)`
    # would be the same stream but rounded differently
    r.data .*= convert(real(eltype(r)), 2 / sqrt(2))
    return r
end

export SpectralAR1Process

"""First-order auto-regressive random process (AR1) in spectral space,
evolving `wavenumbers` with respectice `time_scales` and `standard_deviations`
independently. Transformed after every time step to grid space with a
`clamp` applied to limit extrema. For reproducability `seed` can be
provided and an independent `random_number_generator` is used
that is reseeded on every `initialize!`. `n_patterns > 1` evolves that many
independent patterns (`Spectral3D`/`Grid3D` with `n_patterns` slots).

Every time step the coefficients `ξ` of degree l of the pattern are updated as

    ξ ← a ξ + Fₗ r,   a = exp(-Δt / time_scale)

with complex normals `r` (unit variance in real and imaginary part) and the noise factors
(ECMWF Tech Memorandum 598, Appendix 8)

    Fₗ = F₀ exp(-l(l+1) / (2k(k+1))),   F₀ = σ √(1 - a²) norm_sphere / √(2 ∑ₗ (2l+1) exp(-l(l+1) / (k(k+1))))

with `k = wavenumber`, `σ = standard_deviation` and the sum over the degrees that get noise,
1 ≤ l ≤ lmax (only the global mean l = 0 gets none). `F₀` normalises the stationary pattern
to the grid-point standard deviation `σ`, the factor √(1 - a²) compensates for the memory
(the stationary variance is Fₗ²/(1 - a²)). The variance decays like a Gaussian in l: the
large scales keep their variance, at l = k it has dropped to 1/e, so features have a size
of about radius/k (k = 12: about 500 km). This limits the noise to the large scales, e.g. for
the SPPT stochastic physics scheme. `wavenumber` has to be positive: unlike for
`SpectralWhiteNoiseProcess`, 0 does not give a spectrally flat pattern.

Fields are: $(TYPEDFIELDS)"""
@kwdef struct SpectralAR1Process{NF, VectorType, S, RNG, IntType, SpectralVariable} <: AbstractRandomProcess
    truncation::IntType

    "[OPTION] Time scale of the AR1 process"
    time_scale::Second = Hour(6)

    "[OPTION] Wavenumber (0-based) of the AR1 process"
    wavenumber::IntType = 12

    "[OPTION] Number of independent patterns, 1 = a single pattern"
    n_patterns::IntType = 1

    "[OPTION] Standard deviation of the AR1 process"
    standard_deviation::NF = 1 / 3

    "[OPTION] Range to clamp values into after every transform into grid space"
    clamp::NTuple{2, NF} = (-1, 1)

    "[OPTION] Random number generator seed, 0=randomly seed from Julia's GLOBAL_RNG"
    seed::S = 0

    "Independent random number generator for this random process, see `architecture_rng`"
    random_number_generator::RNG = Random.Xoshiro(seed)

    "Precomputed auto-regressive factor [1], function of time scale and model time step"
    autoregressive_factor::Base.RefValue{NF} = Ref(zero(NF))

    "Precomputed noise factors [1] for every total wavenumber l"
    noise_factors::VectorType = zeros(NF, truncation + 1)

    "Buffer for the random numbers drawn every time step, same size as the pattern"
    random_numbers::SpectralVariable
end

# generator function
function SpectralAR1Process(
        SG::SpectralGrid;
        seed = 0,
        n_patterns = 1,
        draw_on_cpu::Bool = false,
        random_number_generator = architecture_rng(SG.architecture, seed; draw_on_cpu),
        kwargs...,
    )
    RNG = typeof(random_number_generator)
    (; truncation, spectrum) = SG
    random_numbers = n_patterns <= 1 ? zeros(SG.SpectralVariable2D, spectrum) :
        zeros(SG.SpectralVariable3D, spectrum, n_patterns)
    return SpectralAR1Process{SG.NF, SG.VectorType, typeof(seed), RNG, typeof(truncation), typeof(random_numbers)}(;
        truncation, n_patterns, seed, random_number_generator, random_numbers, kwargs...,
    )
end

"""$(TYPEDSIGNATURES)
Variables of a spectral random process: the spectral `random_pattern` and its grid counterpart,
2D for a single pattern, 3D with `n` slots for `n > 1` independent patterns."""
function random_pattern_variables(n::Integer)
    desc = "Random pattern for the random process"
    return n <= 1 ? (
            PrognosticVariable(:random_pattern, Spectral2D(); desc, units = "1"),
            GridVariable(:random_pattern, Grid2D(); desc, units = "1"),
        ) : (
            PrognosticVariable(:random_pattern, Spectral3D(; n); desc, units = "1"),
            GridVariable(:random_pattern, Grid3D(; n); desc, units = "1"),
        )
end

variables(process::SpectralAR1Process) = random_pattern_variables(process.n_patterns)

function initialize!(
        process::SpectralAR1Process,
        model::AbstractModel,
    )
    # auto-regressive factor in the AR1 process
    dt = model.time_stepping.Δt             # in seconds
    process.autoregressive_factor[] = exp(-dt / Second(process.time_scale).value)

    # noise factors per total wavenumber in the AR1 process
    k = process.wavenumber
    a = process.autoregressive_factor[]
    σ = process.standard_deviation

    # ECMWF Tech Memorandum 598, Appendix 8, eq. 18
    # TODO *norm_sphere seems to be needed, maybe ECMWF uses another normalization of the harmonics?
    F₀_denominator = 2 * sum([(2l + 1) * exp(-l * (l + 1) / (k * (k + 1))) for l in 1:process.truncation])
    F₀ = sqrt(σ^2 * (1 - a^2) / F₀_denominator) * model.spectral_transform.norm_sphere

    # ECMWF Tech Memorandum 598, Appendix 8, eq. 17
    # eigenvalue = l * (l - 1), 1-based l so l=1 (mean) gets factor 0 instead of F₀
    NF = eltype(process.noise_factors)
    ls = 1:length(process.noise_factors)
    denom = NF(2k * (k + 1))
    F₀_NF = NF(F₀)
    @. process.noise_factors = ifelse(ls > 1, F₀_NF * exp(-ls * (ls - 1) / denom), zero(NF))

    # reseed the random number generator, for seed=0 randomly seed from Julia's global RNG
    seed = process.seed == 0 ? rand(UInt) : process.seed
    Random.seed!(process.random_number_generator, seed)
    return nothing
end

function update_random!(
        random_pattern::LowerTriangularArray,
        process::SpectralAR1Process{NF},
    ) where {NF}
    a = process.autoregressive_factor[]
    r = complex_normals!(process.random_numbers, process.random_number_generator)

    # linear index over all patterns, l_indices are per lm
    nlm = size(random_pattern, 1)
    launch!(
        architecture(random_pattern), SpectralWorkOrder, size(random_pattern), spectral_ar1_kernel!,
        random_pattern, r, process.noise_factors, random_pattern.spectrum.l_indices, a, nlm,
    )
    return nothing
end

@kernel inbounds = true function spectral_ar1_kernel!(
        random_pattern, r, noise_factors, l_indices, a, nlm,
    )
    i = @index(Global, Linear)
    l = l_indices[mod1(i, nlm)]
    ξ = noise_factors[l]
    random_pattern[i] = a * random_pattern[i] + ξ * r[i]
end

export SpectralWhiteNoiseProcess, GridWhiteNoiseProcess

"""Random process that is white in time in spectral space: the pattern is redrawn
from scratch every time step, it has no memory. Same units and `standard_deviation`
semantics as `SpectralAR1Process`.

- `wavenumber = 0` (default) gives a spectrally flat pattern, `wavenumber = k > 0` the
  ECMWF shaping (Tech Memorandum 598, Appendix 8) of the AR1 process: the real and
  imaginary parts of the m > 0 coefficients of degree l have the standard deviation

      Fₗ = F₀ exp(-l(l+1) / (2k(k+1))),   F₀ = σ norm_sphere / √(2 ∑ₗ (2l+1) exp(-l(l+1) / (k(k+1))))

  with `σ = standard_deviation` and the sum over the degrees that get noise, 1 ≤ l < lmax.
  `F₀` normalises the pattern to the grid-point standard deviation `σ`. The variance
  Fₗ² decays like a Gaussian in l: the large scales keep their
  variance, at l = k it has dropped to 1/e, so features have a size of about radius/k. This
  may be used e.g. in the SPPT stochastic physics scheme to limit the noise to the large scales.
  With this option the pattern is white in time, but coloured in space.
  `wavenumber = 0` stands for k → ∞ (exp → 1, the same variance at every degree),
  not for wavenumber zero (used e.g. for the GFDT sensitivity computations). With k > 0
  the pattern is correlated in space but remains white in time (independent draws).
- Neither the global mean (l = 0) nor the last degree l = lmax (not a resolved degree
  in SpeedyWeather, its coefficients do not survive a grid round trip) get noise.
- The m = 0 coefficients are real, with their real part scaled by √2 so that the
  pattern is statistically homogeneous and a band-limited grid field: a spectral to
  grid to spectral round trip returns it to rounding.
- No clamp by default, `clamp = (-Inf, Inf)`: clamping makes the increment
  non-Gaussian and the grid pattern no longer band-limited.

Fields are $(TYPEDFIELDS)"""
@kwdef struct SpectralWhiteNoiseProcess{NF, VectorType, S, RNG, IntType} <: AbstractRandomProcess
    truncation::IntType

    "[OPTION] Wavenumber (0-based) of the spectral shaping, 0 = spectrally flat"
    wavenumber::IntType = 0

    "[OPTION] Number of independent patterns, 1 = a single pattern"
    n_patterns::IntType = 1

    "[OPTION] Standard deviation of the pattern in grid space"
    standard_deviation::NF = 1 / 3

    "[OPTION] Range to clamp values into after every transform into grid space, default no clamping"
    clamp::NTuple{2, NF} = (-Inf, Inf)

    "[OPTION] Random number generator seed, 0 = randomly seed from Julia's GLOBAL_RNG"
    seed::S = 0

    "Independent random number generator for this random process, see `architecture_rng`"
    random_number_generator::RNG = Random.Xoshiro(seed)

    "[DERIVED] Noise factors [1] for every total wavenumber l"
    noise_factors::VectorType = zeros(NF, truncation + 1)
end

# generator function
function SpectralWhiteNoiseProcess(
        SG::SpectralGrid;
        seed = 0,
        draw_on_cpu::Bool = false,
        random_number_generator = architecture_rng(SG.architecture, seed; draw_on_cpu),
        kwargs...,
    )
    RNG = typeof(random_number_generator)
    (; truncation) = SG
    return SpectralWhiteNoiseProcess{SG.NF, SG.VectorType, typeof(seed), RNG, typeof(truncation)}(;
        truncation, seed, random_number_generator, kwargs...,
    )
end

variables(process::SpectralWhiteNoiseProcess) = random_pattern_variables(process.n_patterns)

function initialize!(process::SpectralWhiteNoiseProcess, model::AbstractModel)
    σ = process.standard_deviation
    k = process.wavenumber
    NF = eltype(process.noise_factors)
    ls = 1:length(process.noise_factors)    # 1-based degree l + 1

    # no noise on the global mean (l = 0) and the last degree l = lmax
    live = (ls .> 1) .& (ls .< length(ls))

    # computes the noise factors for the wavenumbers as explained in the docstring of 
    # `SpectralWhiteNoiseProcess`
    if k <= 0
        denom = 2 * sum(2l + 1 for l in 1:(process.truncation - 1))
        F₀ = sqrt(σ^2 / denom) * model.spectral_transform.norm_sphere
        copyto!(process.noise_factors, NF[live[i] ? NF(F₀) : zero(NF) for i in ls])
    else
        denom0 = 2 * sum(
            (2l + 1) * exp(-l * (l + 1) / (k * (k + 1)))
                for l in 1:(process.truncation - 1)
        )
        F₀ = sqrt(σ^2 / denom0) * model.spectral_transform.norm_sphere
        dn = NF(2k * (k + 1))
        copyto!(
            process.noise_factors,
            NF[live[i] ? NF(F₀) * exp(-NF(i) * NF(i - 1) / dn) : zero(NF) for i in ls],
        )
    end

    # reseed the random number generator, for seed=0 randomly seed from Julia's global RNG
    seed = process.seed == 0 ? rand(UInt) : process.seed
    Random.seed!(process.random_number_generator, seed)
    return nothing
end

"""$(TYPEDSIGNATURES)
Draw a new white noise `pattern` (all its slots) from `process`, independently of
the previous pattern."""
function update_random!(
        pattern::LowerTriangularArray,
        process::SpectralWhiteNoiseProcess{NF},
    ) where {NF}
    complex_normals!(pattern, process.random_number_generator)

    # linear index over all slots, l_indices and m_indices are per lm
    nlm = size(pattern, 1)
    launch!(
        architecture(pattern), SpectralWorkOrder, size(pattern), white_noise_kernel!,
        pattern, process.noise_factors, pattern.spectrum.l_indices, pattern.spectrum.m_indices, nlm,
    )
    return nothing
end

@kernel inbounds = true function white_noise_kernel!(
        pattern, noise_factors, l_indices, m_indices, nlm,
    )
    i = @index(Global, Linear)
    lm = mod1(i, nlm)
    l = l_indices[lm]
    x = noise_factors[l] * pattern[i]

    # m = 0 is real (a real grid field has real zonal-mean coefficients), and carries √2 × the
    # real part to restore the variance of the dropped imaginary part, keeping the pattern homogeneous
    rt2 = oftype(real(x), 1.4142135623730951)
    pattern[i] = ifelse(m_indices[lm] == 1, complex(rt2 * real(x)), x)
end

"""Random process that is white in time and in space on the grid: every grid point
of every pattern is an independent normal with standard deviation `standard_deviation`,
redrawn every time step. As a random process it fills `vars.grid.random_pattern` directly
(nothing is transformed). On a reduced grid the boxes differ in area, so the pattern is
uniform in the grid coordinates, not statistically homogeneous on the sphere.
Fields are $(TYPEDFIELDS)"""
@kwdef struct GridWhiteNoiseProcess{NF, S, RNG, IntType} <: AbstractRandomProcess
    "[OPTION] Number of independent patterns"
    n_patterns::IntType = 1

    "[OPTION] Standard deviation of the pattern at every grid point"
    standard_deviation::NF = 1 / 3

    "[OPTION] Random number generator seed, 0 = randomly seed from Julia's GLOBAL_RNG"
    seed::S = 0

    "Independent random number generator for this random process, see `architecture_rng`"
    random_number_generator::RNG = Random.Xoshiro(seed)

end

function GridWhiteNoiseProcess(
        SG::SpectralGrid;
        seed = 0,
        draw_on_cpu::Bool = false,
        random_number_generator = architecture_rng(SG.architecture, seed; draw_on_cpu),
        kwargs...,
    )
    RNG = typeof(random_number_generator)
    return GridWhiteNoiseProcess{SG.NF, typeof(seed), RNG, typeof(SG.truncation)}(; seed, random_number_generator, kwargs...)
end

function variables(process::GridWhiteNoiseProcess)
    n = process.n_patterns
    dims = n <= 1 ? Grid2D() : Grid3D(; n)
    return (GridVariable(:random_pattern, dims, desc = "Random pattern for the random process", units = "1"),)
end

function initialize!(process::GridWhiteNoiseProcess, model::AbstractModel)
    seed = process.seed == 0 ? rand(UInt) : process.seed
    Random.seed!(process.random_number_generator, seed)
    return nothing
end

# the pattern is drawn in grid space, nothing to transform
SpeedyTransforms.transform!(::Variables, ::GridWhiteNoiseProcess, ::AbstractSpectralTransform) = nothing

# as `model.random_process` the pattern lives on the grid
update_random!(vars::Variables, process::GridWhiteNoiseProcess) =
    update_random!(vars.grid.random_pattern, process)

"""$(TYPEDSIGNATURES)
Draw a new white noise `pattern` (all its slots) from `process`."""
function update_random!(pattern::AbstractField, process::GridWhiteNoiseProcess)
    data = parent(pattern)
    Random.randn!(process.random_number_generator, data)
    data .*= process.standard_deviation
    return nothing
end
