abstract type AbstractRandomProcess <: AbstractModelComponent end

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
random_process!(::Variables, process::Nothing) = nothing

export SpectralAR1Process

"""First-order auto-regressive random process (AR1) in spectral space,
evolving `wavenumbers` with respectice `time_scales` and `standard_deviations`
independently. Transformed after every time step to grid space with a
`clamp` applied to limit extrema. For reproducability `seed` can be
provided and an independent `random_number_generator` is used
that is reseeded on every `initialize!`. Fields are: $(TYPEDFIELDS)"""
@kwdef struct SpectralAR1Process{NF, VectorType, S, RNG, IntType, RefV, TS} <: AbstractRandomProcess
    trunc::IntType

    "[OPTION] Time scale of the AR1 process"
    time_scale::TS = Hour(6)

    "[OPTION] Wavenumber of the AR1 process"
    wavenumber::IntType = 12

    "[OPTION] Standard deviation of the AR1 process"
    standard_deviation::NF = 1 / 3

    "[OPTION] Range to clamp values into after every transform into grid space"
    clamp::NTuple{2, NF} = (-1, 1)

    "[OPTION] Random number generator seed, 0=randomly seed from Julia's GLOBAL_RNG"
    seed::S = 0

    "Independent random number generator for this random process"
    random_number_generator::RNG = Random.Xoshiro(seed)

    "Precomputed auto-regressive factor [1], function of time scale and model time step"
    autoregressive_factor::RefV = Ref(zero(NF))

    "Precomputed noise factors [1] for every total wavenumber l"
    noise_factors::VectorType = zeros(NF, trunc + 2)
end

# generator function
function SpectralAR1Process(SG::SpectralGrid; kwargs...)
    RNG = haskey(kwargs, :random_number_generator) ? typeof(kwargs[:random_number_generator]) : typeof(Random.Xoshiro())
    SeedType = haskey(kwargs, :seed) ? typeof(kwargs[:seed]) : Int
    return SpectralAR1Process{SG.NF, SG.VectorType, SeedType, RNG, typeof(SG.trunc), Base.RefValue{SG.NF}, Dates.Second}(trunc = SG.trunc; kwargs...)
end

function variables(::SpectralAR1Process)
    return (
        PrognosticVariable(:random_pattern, Spectral2D(), desc = "Random pattern for the random process", units = "1"),
        GridVariable(:random_pattern, Grid2D(), desc = "Random pattern for the random process", units = "1"),
    )
end

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
    F₀_denominator = 2 * sum([(2l + 1) * exp(-l * (l + 1) / (k * (k + 1))) for l in 1:process.trunc])
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

function random_process!(
        vars::Variables,
        process::SpectralAR1Process{NF},
    ) where {NF}

    (; random_pattern) = vars.prognostic
    a = process.autoregressive_factor[]
    RNG = process.random_number_generator
    s = convert(NF, 2 / sqrt(2))              # to scale: std(real(randn(Complex))) = √2/2 to 1

    arch = architecture(random_pattern)
    n = length(random_pattern)                # total number of harmonics

    # draw all complex normals on CPU through the seeded RNG to keep reproducibility,
    # then transfer to the device in one go (GPU has no compatible Random.Xoshiro RNG)
    r_cpu = Vector{Complex{NF}}(undef, n)
    @inbounds for i in 1:n
        r_cpu[i] = s * randn(RNG, Complex{NF})
    end
    r = on_architecture(arch, r_cpu)

    launch!(
        arch, SpectralWorkOrder, size(random_pattern), spectral_ar1_kernel!,
        random_pattern, r, process.noise_factors, random_pattern.spectrum.l_indices, a
    )
    return nothing
end

@kernel inbounds = true function spectral_ar1_kernel!(
        random_pattern, r, noise_factors, l_indices, a
    )
    lm = @index(Global, Linear)
    l = l_indices[lm]
    ξ = noise_factors[l]
    random_pattern[lm] = a * random_pattern[lm] + ξ * r[lm]
end
