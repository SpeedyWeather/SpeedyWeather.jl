abstract type AbstractRandomProcess <: AbstractModelComponent end

"""$(TYPEDSIGNATURES)
General transform for `random processes <: AbstractRandomProcess`.
Takes the spectral `random_pattern` in the prognostic variables
and transforms it to spectral space in `diagn.grid.random_pattern`."""
function SpeedyTransforms.transform!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Integer,
    random_process::AbstractRandomProcess,
    spectral_transform::SpectralTransform,
)
    grid = diagn.grid.random_pattern
    spec = progn.random_pattern
    transform!(grid, spec, diagn, spectral_transform)

    if :clamp in fieldnames(typeof(random_process))
        clamp!(grid, random_process.clamp...)
    end
end

export NoRandomProcess
"""Dummy type for no random process."""
struct NoRandomProcess <: AbstractRandomProcess end
NoRandomProcess(::SpectralGrid) = NoRandomProcess()

"""$(TYPEDSIGNATURES)
`NoRandomProcess` does not need to transform any random pattern from
spectral to grid space."""
function SpeedyTransforms.transform!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Integer,
    random_process::NoRandomProcess,
    spectral_transform::SpectralTransform,
)
    return nothing
end

initialize!(process::NoRandomProcess, model::AbstractModel) = nothing
random_process!(progn::PrognosticVariables, process::NoRandomProcess) = nothing

export SpectralAR1Process

"""First-order auto-regressive random process (AR1) in spectral space,
evolving `wavenumbers` with respectice `time_scales` and `standard_deviations`
independently. Transformed after every time step to grid space with a
`clamp` applied to limit extrema. For reproducability `seed` can be
provided and an independent `random_number_generator` is used
that is reseeded on every `initialize!`. Fields are $(TYPEDFIELDS)"""
@kwdef struct SpectralAR1Process{NF, VectorType} <: AbstractRandomProcess
    trunc::Int
    
    "[OPTION] Time scale of the AR1 process"
    time_scale::Second = Hour(6)

    "[OPTION] Wavenumber of the AR1 process"
    wavenumber::Int = 12

    "[OPTION] Standard deviation of the AR1 process"
    standard_deviation::NF = 1/3

    "[OPTION] Range to clamp values into after every transform into grid space"
    clamp::NTuple{2, NF} = (-1, 1)

    "[OPTION] Random number generator seed, 0=randomly seed from Julia's GLOBAL_RNG"
    seed::Int = 0

    "Independent random number generator for this random process"
    random_number_generator::Random.Xoshiro = Random.Xoshiro(seed)

    "Precomputed auto-regressive factor [1], function of time scale and model time step"
    autoregressive_factor::Base.RefValue{NF} = Ref(zero(NF))

    "Precomputed noise factors [1] for every total wavenumber l"
    noise_factors::VectorType = zeros(NF, trunc+2)
end

# generator function
SpectralAR1Process(SG::SpectralGrid; kwargs...) = SpectralAR1Process{SG.NF, SG.VectorType}(trunc=SG.trunc; kwargs...)

function initialize!(
    process::SpectralAR1Process,
    model::AbstractModel,
)
    # auto-regressive factor in the AR1 process
    dt = model.time_stepping.Δt_sec         # in seconds
    process.autoregressive_factor[] = exp(-dt/Second(process.time_scale).value)

    # noise factors per total wavenumber in the AR1 process
    k = process.wavenumber
    a = process.autoregressive_factor[]
    σ = process.standard_deviation

    # ECMWF Tech Memorandum 598, Appendix 8, eq. 18
    # TODO *norm_sphere seems to be needed, maybe ECMWF uses another normalization of the harmonics?
    F₀_denominator = 2*sum([(2l + 1)*exp(-l*(l+1)/(k*(k+1))) for l in 1:process.trunc])
    F₀ = sqrt(σ^2 * (1-a^2) / F₀_denominator)*model.spectral_transform.norm_sphere

    for l in eachindex(process.noise_factors)       # total wavenumber, but 1-based
        eigenvalue = l*(l-1)                        # (negative) eigenvalue l*(l+1) but 1-based l->l-1

        # ECMWF Tech Memorandum 598, Appendix 8, eq. 17
        process.noise_factors[l] = F₀*exp(-eigenvalue/(2k*(k+1)))
    end

    # set mean of random pattern to zero
    process.noise_factors[1] = 0

    # reseed the random number generator, for seed=0 randomly seed from Julia's global RNG
    seed = process.seed == 0 ? rand(UInt) : process.seed
    Random.seed!(process.random_number_generator, seed)
    return nothing
end

function random_process!(
    progn::PrognosticVariables,
    process::SpectralAR1Process{NF},
) where NF

    (; random_pattern) = progn
    lmax, mmax = size(random_pattern, OneBased, as=Matrix)  # max degree l, order m of harmonics (1-based)

    a = process.autoregressive_factor[]
    RNG = process.random_number_generator
    s = convert(NF, 2/sqrt(2))              # to scale: std(real(randn(Complex))) = √2/2 to 1

    lm = 0
    @inbounds for m in 1:mmax
        for l in m:lmax
            lm += 1

            # draw from independent N(0,1) in real and imaginary parts
            r = s*randn(RNG, Complex{NF})   # scale to unit variance in real/imaginary

            # ECMWF Tech Memorandum 598, Appendix 8, eq. 14
            ξ = process.noise_factors[l]
            random_pattern[lm] *= a         # auto-regressive term
            random_pattern[lm] += ξ*r       # noise term
        end
    end
end